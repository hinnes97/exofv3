! $Id: fv_control.F90,v 18.0.2.1 2010/05/13 17:36:18 pjp Exp $
!
!----------------
! FV contro panel
!----------------

module fv_control_mod

   use constants_mod,      only: pi, kappa, radius
   use field_manager_mod,  only: MODEL_ATMOS
   use mpp_mod,            only: FATAL, mpp_error, mpp_pe, stdlog, &
                                 mpp_npes, mpp_get_current_pelist, get_unit
   use mpp_domains_mod,    only: mpp_get_data_domain, mpp_get_compute_domain
   use tracer_manager_mod, only: tm_get_number_tracers => get_number_tracers, &
                                 tm_get_tracer_index   => get_tracer_index,   &
                                 tm_get_tracer_indices => get_tracer_indices, &
                                 tm_set_tracer_profile => set_tracer_profile, &
                                 tm_get_tracer_names   => get_tracer_names,   &
                                 tm_check_if_prognostic=> check_if_prognostic,&
                                 tm_register_tracers   => register_tracers

   use fv_io_mod,          only: fv_io_exit
   use fv_restart_mod,     only: fv_restart_init, fv_restart_end
   use fv_arrays_mod,      only: fv_atmos_type
   use fv_grid_utils_mod,  only: grid_utils_init, grid_utils_end, ptop_min, deglat, &
                                 da_min_c, da_min
   use fv_grid_tools_mod,  only: init_grid, cosa, sina, area, area_c, dx, dy, dxa, dya, &
                                 dxc, dyc, grid_type, dx_const, dy_const,                         &
                                 deglon_start, deglon_stop, deglat_start, deglat_stop, &
                                 read_grid, debug_message_size, write_grid_char_file
   use fv_mp_mod,          only: mp_start, domain_decomp, domain, &
                                 ng, tile, npes_x, npes_y, gid, io_domain_layout
   use test_cases_mod,     only: test_case, alpha
   use fv_timing_mod,      only: timing_on, timing_off, timing_init, timing_prt

   implicit none
   private

!-----------------------------------------------------------------------
! Grid descriptor file setup
!-----------------------------------------------------------------------
   character(len=80) :: grid_name = 'Gnomonic'
   character(len=120):: grid_file = 'Inline'
!   integer      :: grid_type = 0    ! -1: read from file; 0: ED Gnomonic
!                                    !  0: the "true" equal-distance Gnomonic grid
!                                    !  1: the traditional equal-distance Gnomonic grid
!                                    !  2: the equal-angular Gnomonic grid
!                                    !  3: the lat-lon grid -- to be implemented
!                                    !  4: double periodic boundary condition on Cartesian grid
!                                    !  5: channel flow on Cartesian grid
!  -> moved to grid_tools

! Momentum (or KE) options:
   integer :: hord_mt = 9    ! the best option for Gnomonic grids  
   integer :: kord_mt = 8    ! vertical mapping option

! Vorticity & w transport options:
   integer :: hord_vt = 9    ! 10 not recommended (noisy case-5) 

! Heat & air mass (delp) transport options:
   integer :: hord_tm = 9    ! virtual potential temperature
   integer :: hord_dp = 9    ! delp (positive definite)
   integer :: hord_ze = 9    ! edge height (can be negative near the surface)
   integer :: kord_tm =-8    !

! Tracer transport options:
   integer :: hord_tr = 12   !11: PPM mono constraint (Lin 2004); fast 
                             !12: Huynh 2nd constraint (Lin 2004) +
                             !    positive definite (Lin & Rood 1996); slower
                             !>12: positive definite only (Lin & Rood 1996); fastest
   integer :: kord_tr = 8    ! 

   integer :: nord=3         ! 0: del-2, 1: del-4, 2: del-6, 3: del-8 divergence damping
                             ! Alternative setting for high-res: nord=1; d4_bg = 0.075
   real    :: dddmp = 0.0    ! coefficient for del-2 divergence damping (0.2)
                             ! for C90 or lower: 0.2
   real    :: d2_bg = 0.0    ! coefficient for background del-2 divergence damping
   real    :: d4_bg = 0.16   ! coefficient for background del-4(6) divergence damping
                             ! for stability, d4_bg must be <=0.16 if nord=3
   real    :: vtdm4 = 0.0    ! coefficient for del-4 vorticity damping
   real    :: d2_bg_k1 = 4.         ! factor for d2_bg (k=1)
   real    :: d2_bg_k2 = 2.         ! factor for d2_bg (k=2)
   real    :: d2_divg_max_k1 = 0.05 ! d2_divg max value (k=1)
   real    :: d2_divg_max_k2 = 0.02 ! d2_divg max value (k=2)
   real    :: damp_k_k1 = 0.05      ! damp_k value (k=1)
   real    :: damp_k_k2 = 0.025     ! damp_k value (k=2)

! PG off centering:
   real    :: beta1 = 0.0    !
   real    :: beta  = 0.25   ! 0.5 is "neutral" but it may not be stable
#ifdef SW_DYNAMICS
   integer :: n_sponge = 0   ! Number of sponge layers at the top of the atmosphere
   real    :: d_ext = 0.    
   integer :: nwat  = 0      ! Number of water species
   logical :: warm_start = .false. 
#else
   integer :: n_sponge = 1   ! Number of sponge layers at the top of the atmosphere
   real    :: d_ext = 0.02   ! External model damping (was 0.02)
   integer :: nwat  = 3      ! Number of water species
                             ! Set to .F. if cold_start is desired (including terrain generation)
   logical :: warm_start = .true. 
#endif
   integer :: m_riem  = 0    ! Time scheme for Riem solver subcycling
   integer :: k_top   = 1    ! Starting layer for non-hydrostatic dynamics
   integer :: n_split = 0    ! Number of time splits for the lagrangian dynamics
                             ! Default = 0 (automatic computation of best value)
   integer :: m_split = 0    ! Number of time splits for Riemann solver
   integer :: k_split = 1    ! Number of time splits for Remapping

!            For doubly periodic domain with sim_phys
!                     5km        150         20 (7.5 s)  2
!
!                     Estimates for Gnomonic grids:
            !===================================================
            !        dx (km)    dt (sc)    n_split    m_split
            !===================================================
            ! C1000:  ~10        150         16          3
            ! C2000:   ~5         90         18 (5 s)    2
            !===================================================
! The nonhydrostatic algorithm is described in Lin 2006, QJ, (submitted)
! C2000 should easily scale to at least 6 * 100 * 100 = 60,000 CPUs  
! For a 1024 system: try 6 x 13 * 13 = 1014 CPUs
  
   integer :: q_split = 0    ! Number of time splits for tracer transport
   integer :: print_freq = 0 ! Print max/min of selected fields
                             ! 0: off
                             ! positive n: every n hours
                             ! negative n: every time step

!------------------------------------------
! Model Domain parameters
!------------------------------------------
   integer :: npx                     ! Number of Grid Points in X- dir
   integer :: npy                     ! Number of Grid Points in Y- dir
   integer :: npz                     ! Number of Vertical Levels
   integer :: npz_rst = 0             ! Original Vertical Levels (in the restart)
                                      ! 0: no change (default)
   integer :: layout(2)=(/2,2/)       ! Processor layout
   integer :: io_layout(2)=(/0,0/)    ! IO domain processor layout
   integer :: ncnst = 0               ! Number of advected consituents
   integer :: pnats = 0               ! Number of non-advected consituents
   integer :: ntiles                  ! Number or tiles that make up the Grid 
   integer :: ntilesMe                ! Number of tiles on this process =1 for now
   integer, parameter:: ndims = 2     ! Lat-Lon Dims for Grid in Radians
   integer :: nf_omega  = 1           ! Filter omega "nf_omega" times
   integer :: fv_sg_adj = -1          ! Perform grid-scale dry adjustment if > 0
                                      ! Relaxzation time  scale (sec) if positive
#ifdef MARS_GCM
   real    :: p_ref = 600.
   real    :: reference_sfc_pres = 7.7E2
   real    :: sponge_damp=   1.0
   real    :: dry_mass = 7.7E2
#else
   real    :: p_ref = 1.E5
   real    :: dry_mass = 98290.
#endif
   integer :: nt_prog = 0
   integer :: nt_phys = 0
   real    :: tau_h2o = 0.            ! Time scale (days) for ch4_chem

   real    :: too_big  = 1.E35
   real    :: d_con = 0.
   real    :: consv_te = 0.
   real    :: tau = 0.                ! Time scale (days) for Rayleigh friction
   real    :: rf_center = 0.          ! Center position of the hyper-tan profile
                                      ! 0: use the top layer center
                                      ! > 0, [Pascal]
   real    :: rf_lowest = 1.e1                                      
   logical :: tq_filter = .false.
   logical :: filter_phys = .false.
   logical :: dwind_2d = .false.
   logical :: inline_q = .false.
   logical :: breed_vortex_inline = .false.
   logical :: no_cgrid = .false.
   logical :: init_wind_m = .false.
   logical :: range_warn = .false.
   logical :: fill = .false.
   logical :: fill_dp = .false.
   logical :: non_ortho = .true.
   logical :: adiabatic = .false.     ! Run without physics (full or idealized).
   logical :: moist_phys = .true.     ! Run with moist physics
   logical :: do_Held_Suarez = .false.
   logical :: uniform_vert_spacing=.false. ! Uniform vertical spacing for any npz
   logical :: reproduce_sum = .true.  ! Make global sum for consv_te reproduce
   logical :: adjust_dry_mass = .false.
   logical :: fv_debug  = .false.
   logical :: srf_init  = .false.
   logical :: mountain  = .true.
   logical :: uniform_ppm = .true.
   logical :: remap_t  = .true.
   logical :: z_tracer = .true.       ! transport tracers layer by layer with independent
                                      ! time split; use this if tracer number is huge and/or
                                      ! high resolution (nsplt > 1)

   logical :: non_dilute = .false.
   logical :: old_divg_damp = .false. ! parameter to revert damping parameters back to values
                                      ! defined in a previous revision
                                      ! old_values:
                                      !    d2_bg_k1 = 6.           d2_bg_k2 = 4.
                                      !    d2_divg_max_k1 = 0.02   d2_divg_max_k2 = 0.01
                                      !    damp_k_k1 = 0.          damp_k_k2 = 0.
                                      ! current_values:
                                      !    d2_bg_k1 = 4.           d2_bg_k2 = 2.
                                      !    d2_divg_max_k1 = 0.05   d2_divg_max_k2 = 0.02
                                      !    damp_k_k1 = 0.05        damp_k_k2 = 0.025
   logical :: master

   logical :: fv_land = .false.       ! To cold starting the model with USGS terrain
!--------------------------------------------------------------------------------------
! The following options are useful for NWP experiments using datasets on the lat-lon grid
!--------------------------------------------------------------------------------------
   logical :: nudge = .false.         ! Perform nudging
   logical :: ncep_ic = .false.       ! use NCEP ICs 
   logical :: fv_diag_ic = .false.    ! reconstruct IC from fv_diagnostics on lat-lon grid
   logical :: external_ic = .false.   ! use ICs from external sources; e.g. lat-lon FV core
                                      ! or NCEP re-analysis; both vertical remapping & horizontal
                                      ! (lat-lon to cubed sphere) interpolation will be done
! Default restart files from the "Memphis" latlon FV core:
   character(len=128) :: res_latlon_dynamics = 'INPUT/fv_rst.res.nc'
   character(len=128) :: res_latlon_tracers  = 'INPUT/atmos_tracers.res.nc'
! The user also needs to copy the "cold start" cubed sphere restart files (fv_core.res.tile1-6)
! to the INPUT dir during runtime
!------------------------------------------------
! Parameters related to non-hydrostatic dynamics:
!------------------------------------------------
   logical :: hydrostatic = .true.
   logical :: phys_hydrostatic = .true.    ! heating/cooling term from the physics is hydrostatic
   logical :: hybrid_z    = .false. ! use hybrid_z for remapping
   logical :: quick_p_c   = .false. ! Use quick (approximated) algorithm for Riemann Solver (C grid)
   logical :: quick_p_d   = .false. ! Use quick (approximated) algorithm for Riemann Solver (D grid)
                                    ! The above two options run much faster; but it may be unstable
   logical :: Make_NH     = .false. ! Initialize (w, delz) from hydro restart file 
   logical :: make_hybrid_z  = .false. ! transform hydrostatic eta-coord IC into non-hydrostatic hybrid_z
   integer :: m_grad_p = 0   ! method for non-hydrostatic grad-p
                             ! m_grad_p=1:  one-stage full pressure for grad_p; this option is faster
                             !              but it is not suitable for low horizontal resolution
                             ! m_grad_p=0:  two-stage grad computation (best for low resolution runs)
   integer :: a2b_ord = 4    ! order for interpolation from A to B Grid (corners)
   integer :: c2l_ord = 4    ! order for interpolation from D to lat-lon A winds for phys & output
   real    :: ppm_limiter = 2.
   public :: npx,npy,npz, npz_rst, ntiles, ncnst, pnats, nwat
   public :: hord_mt, hord_vt, kord_mt, hord_tm, hord_dp, hord_ze, kord_tm, hord_tr, kord_tr
   public :: nord, no_cgrid, fill_dp, inline_q, breed_vortex_inline, dwind_2d, filter_phys, tq_filter 
   public :: k_split, n_split, m_split, q_split, master
   public :: dddmp, d2_bg, d4_bg, d_ext, vtdm4, beta1, beta, init_wind_m, ppm_limiter
   public :: k_top, m_riem, n_sponge, p_ref, mountain
   public :: uniform_ppm, remap_t,  z_tracer, fv_debug
   public :: external_ic, ncep_ic, fv_diag_ic, res_latlon_dynamics, res_latlon_tracers, fv_land
   public :: fv_sg_adj, tau, tau_h2o, rf_center, rf_lowest, d_con
   public :: fv_init, fv_end
   public :: domain
   public :: adiabatic, nf_omega, moist_phys, range_warn, non_dilute
   public :: hydrostatic, phys_hydrostatic,  hybrid_z, quick_p_c, quick_p_d, m_grad_p, a2b_ord
   public :: nt_prog, nt_phys
   public :: d2_bg_k1, d2_bg_k2, d2_divg_max_k1, d2_divg_max_k2, damp_k_k1, damp_k_k2

#ifdef MARS_GCM
   public :: reference_sfc_pres, sponge_damp
#endif MARS

   integer, allocatable :: pelist(:)
   integer :: commID

 contains

!-------------------------------------------------------------------------------
         
 subroutine fv_init(Atm, dt_atmos)

   type(fv_atmos_type), intent(inout) :: Atm(:)
   real,                intent(in)    :: dt_atmos

   integer :: i, j, k, n
   integer :: isc, iec, jsc, jec
   integer :: isd, ied, jsd, jed
   real :: sdt


! tracers
   integer :: num_family          ! output of register_tracers

! Start up MPI
      allocate( pelist(mpp_npes()) )
      call mpp_get_current_pelist( pelist, commID=commID )
      call mp_start(commID)  ! fv_mp_mod will eventually be eliminated

      master = gid==0

    ! Initialize timing routines
      call timing_init
      call timing_on('TOTAL')

    ! Setup the run from namelist 
      ntilesMe = size(Atm(:))
      if(ntilesMe > 1)call mpp_error(FATAL,'More than 1 tile per process not implemented')

      call run_setup(Atm(1),dt_atmos)   ! initializes domain_decomp
                                        ! needs modification for multiple tiles
      k_top = max(1, k_top)   ! to idiot proof


!--------------------------------------------------
! override number of tracers by reading field_table
!--------------------------------------------------

      call tm_register_tracers (MODEL_ATMOS, ncnst, nt_prog, pnats, num_family)
      if(master) write(*,*) 'ncnst=', ncnst,' num_prog=',nt_prog,' pnats=',pnats,' num_family=',num_family

      Atm(1)%npx=npx; Atm(1)%npy=npy; Atm(1)%npz=npz; Atm(1)%ng=ng
      Atm(1)%npz_rst = npz_rst
      Atm(1)%k_split = k_split
      Atm(1)%n_split = n_split
      Atm(1)%m_split = m_split
      Atm(1)%q_split = q_split
      Atm(1)%print_freq = print_freq
      Atm(1)%consv_te = consv_te
      Atm(1)%ncnst = ncnst
      Atm(1)%pnats = pnats
      Atm(1)%nwat  = nwat 
      Atm(1)%range_warn = range_warn
      Atm(1)%fill = fill
      Atm(1)%warm_start = warm_start
      Atm(1)%tq_filter = tq_filter
      Atm(1)%no_cgrid = no_cgrid
      Atm(1)%init_wind_m = init_wind_m
      Atm(1)%z_tracer = z_tracer
      Atm(1)%do_Held_Suarez = do_Held_Suarez
      Atm(1)%reproduce_sum = reproduce_sum
      Atm(1)%moist_phys = moist_phys
      Atm(1)%srf_init = srf_init
      Atm(1)%mountain = mountain
      Atm(1)%non_ortho = non_ortho
      Atm(1)%adjust_dry_mass = adjust_dry_mass
      Atm(1)%fv_sg_adj = fv_sg_adj
      Atm(1)%dry_mass = dry_mass

      Atm(1)%fv_diag_ic = fv_diag_ic
      Atm(1)%ncep_ic = ncep_ic

      if ( ncep_ic .or. fv_diag_ic)  external_ic = .true.
      Atm(1)%external_ic = external_ic

      Atm(1)%nudge = nudge

      Atm(1)%fv_land     = fv_land
      Atm(1)%res_latlon_dynamics =  res_latlon_dynamics
      Atm(1)%res_latlon_tracers  =  res_latlon_tracers

      Atm(1)%hydrostatic = hydrostatic
      Atm(1)%phys_hydrostatic = phys_hydrostatic
      Atm(1)%hybrid_z    = hybrid_z
      Atm(1)%Make_NH     = Make_NH
      Atm(1)%make_hybrid_z  = make_hybrid_z
      Atm(1)%k_top       = k_top

    ! Read Grid from GRID_FILE and setup grid descriptors
    ! needs modification for multiple tiles
      
      if(grid_type <0 .AND. trim(grid_file) == 'INPUT/grid_spec.nc') then
         call read_grid(Atm(1), grid_name, grid_file, npx, npy, npz, ndims, ntiles, ng)
      else
         call init_grid(Atm(1), grid_name, grid_file, npx, npy, npz, ndims, ntiles, ng)
      endif
      Atm(1)%ndims = ndims
      Atm(1)%ntiles = ntiles

    ! Initialize the SW (2D) part of the model
      call grid_utils_init(Atm(1), Atm(1)%npx, Atm(1)%npy, Atm(1)%npz, uniform_vert_spacing, Atm(1)%grid, Atm(1)%agrid, &
                           area, area_c, cosa, sina, dx, dy, dxa, dya, dxc, dyc, non_ortho,   &
                           uniform_ppm, grid_type, c2l_ord)

      if ( master ) then
           sdt =  dt_atmos/real(n_split*k_split)
           write(*,*) ' '
           write(*,*) 'Divergence damping Coefficients * 1.E6:'
           write(*,*) 'For small dt=', sdt
           write(*,*) 'External mode del-2 (m**2/s)=',  d_ext*da_min_c     /sdt*1.E-6
           write(*,*) 'Internal mode del-2 SMAG dimensionless coeff=',  dddmp
           write(*,*) 'Internal mode del-2 background diff=', d2_bg

           if (nord==1) write(*,*) 'Internal mode del-4 background diff=', d4_bg
           if (nord==2) write(*,*) 'Internal mode del-6 background diff=', d4_bg
           if (nord==3) write(*,*) 'Internal mode del-8 background diff=', d4_bg

           write(*,*) 'Vorticity del-4 (m**4/s)=', (vtdm4*da_min)**2/sdt*1.E-6
           write(*,*) ' '
      endif

      Atm(1)%domain =>domain
      do n = 1, ntilesMe
        call mpp_get_compute_domain(domain,isc,iec,jsc,jec,tile_count=n)
        call mpp_get_data_domain(domain,isd,ied,jsd,jed,tile_count=n)
        if ( (iec-isc+1).lt.4 .or. (jec-jsc+1).lt.4 ) &
           call mpp_error(FATAL,'Domain Decomposition:  Cubed Sphere compute domain has a &
                                &minium requirement of 4 points in X and Y, respectively')
        Atm(n)%isc = isc; Atm(n)%iec = iec
        Atm(n)%jsc = jsc; Atm(n)%jec = jec
        Atm(n)%isd = isd; Atm(n)%ied = ied
        Atm(n)%jsd = jsd; Atm(n)%jed = jed

      ! Allocate State Variables
        allocate (    Atm(n)%u(isd:ied  ,jsd:jed+1,npz) )
        allocate (    Atm(n)%v(isd:ied+1,jsd:jed  ,npz) )

        if ( no_cgrid ) then
             allocate ( Atm(n)%um(isd:ied  ,jsd:jed+1,npz) )
             allocate ( Atm(n)%vm(isd:ied+1,jsd:jed  ,npz) )
        else
             allocate ( Atm(n)%um(1,1,1) )
             allocate ( Atm(n)%vm(1,1,1) )
        endif

        allocate (   Atm(n)%pt(isd:ied  ,jsd:jed  ,npz) )
        allocate ( Atm(n)%delp(isd:ied  ,jsd:jed  ,npz) )
        allocate (    Atm(n)%q(isd:ied  ,jsd:jed  ,npz, ncnst) )

      ! Allocate Auxilliary pressure arrays
        allocate (   Atm(n)%ps(isd:ied  ,jsd:jed) )
        allocate (   Atm(n)%tsurf(isd:ied  ,jsd:jed) )
        allocate (   Atm(n)%pe(isc-1:iec+1, npz+1,jsc-1:jec+1) )
        allocate (   Atm(n)%pk(isc:iec    ,jsc:jec  , npz+1) )
        allocate (   Atm(n)%pk_nd(isc:iec, jsc:jec , npz+1))
        allocate ( Atm(n)%peln(isc:iec,npz+1,jsc:jec) )
        allocate (  Atm(n)%pkz(isc:iec,jsc:jec,npz) )

        allocate ( Atm(n)%u_srf(isc:iec,jsc:jec) )
        allocate ( Atm(n)%v_srf(isc:iec,jsc:jec) )

        if ( fv_land ) then
             allocate ( Atm(n)%sgh(isc:iec,jsc:jec) )
             allocate ( Atm(n)%oro(isc:iec,jsc:jec) )
        else
             allocate ( Atm(n)%oro(1,1) )
        endif

      ! Allocate others
        allocate ( Atm(n)%ts(isc:iec,jsc:jec) )
        allocate ( Atm(n)%phis(isd:ied  ,jsd:jed  ) )
        allocate ( Atm(n)%omga(isd:ied  ,jsd:jed  ,npz) ); Atm(n)%omga=0.
        allocate (   Atm(n)%ua(isd:ied  ,jsd:jed  ,npz) )
        allocate (   Atm(n)%va(isd:ied  ,jsd:jed  ,npz) )
        allocate (   Atm(n)%uc(isd:ied+1,jsd:jed  ,npz) )
        allocate (   Atm(n)%vc(isd:ied  ,jsd:jed+1,npz) )
      ! For tracer transport:
        allocate ( Atm(n)%mfx(isc:iec+1, jsc:jec,  npz) )
        allocate ( Atm(n)%mfy(isc:iec  , jsc:jec+1,npz) )
        allocate (  Atm(n)%cx(isc:iec+1, jsd:jed, npz) )
        allocate (  Atm(n)%cy(isd:ied ,jsc:jec+1, npz) )

!--------------------------
! Non-hydrostatic dynamics:
!--------------------------
      if ( hydrostatic ) then
          allocate (    Atm(n)%w(1, 1  ,1) )
          allocate ( Atm(n)%delz(1, 1  ,1) )
          allocate (  Atm(n)%ze0(1, 1  ,1) )
      else
          allocate (    Atm(n)%w(isd:ied, jsd:jed  ,npz  ) )
          allocate ( Atm(n)%delz(isc:iec, jsc:jec  ,npz) )
          if( hybrid_z ) then
             allocate (  Atm(n)%ze0(isc:iec, jsc:jec ,npz+1) )
          else
             allocate (  Atm(n)%ze0(1, 1  ,1) )
          endif
!         allocate ( mono(isd:ied, jsd:jed, npz))
      endif
 !MDH changed
        Atm(n)%ts   = 300.
        Atm(n)%phis = too_big
! The following statements are to prevent the phatom corner regions from
! growing instability
        Atm(n)%u  = 0.
        Atm(n)%v  = 0.
        Atm(n)%ua = too_big
        Atm(n)%va = too_big

        if ( no_cgrid ) then
             Atm(n)%um = 0.
             Atm(n)%vm = 0.
        endif

      end do
      
    ! Initialize restart functions
      call fv_restart_init()

 end subroutine fv_init
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
         
 subroutine fv_end(Atm)

    type(fv_atmos_type), intent(inout) :: Atm(:)

    integer :: n

    call timing_off('TOTAL')
    call timing_prt( mpp_pe() )

    call fv_restart_end(Atm)
    call fv_io_exit()

  ! Free temporary memory from sw_core routines

  ! Deallocate
    call grid_utils_end( uniform_ppm )

    do n = 1, ntilesMe
      deallocate (    Atm(n)%u )
      deallocate (    Atm(n)%v )
      deallocate (    Atm(n)%um )
      deallocate (    Atm(n)%vm )
      deallocate (   Atm(n)%pt )
      deallocate ( Atm(n)%delp )
      deallocate (    Atm(n)%q )
      deallocate (   Atm(n)%ps )
      deallocate (   Atm(n)%tsurf )
      deallocate (   Atm(n)%pe )
      deallocate (   Atm(n)%pk )
      deallocate ( Atm(n)%peln )
      deallocate (  Atm(n)%pkz )
      deallocate ( Atm(n)%phis )
      deallocate ( Atm(n)%omga )
      deallocate (   Atm(n)%ua )
      deallocate (   Atm(n)%va )
      deallocate (   Atm(n)%uc )
      deallocate (   Atm(n)%vc )
      deallocate ( Atm(n)%mfx )
      deallocate ( Atm(n)%mfy )
      deallocate (  Atm(n)%cx )
      deallocate (  Atm(n)%cy )
      deallocate (  Atm(n)%ak )
      deallocate (  Atm(n)%bk )

      deallocate ( Atm(n)%u_srf )
      deallocate ( Atm(n)%v_srf )
      if( fv_land ) deallocate ( Atm(n)%sgh )
      deallocate ( Atm(n)%oro )

! Non-hydrostatic:
      deallocate ( Atm(n)%w )
      deallocate ( Atm(n)%delz  )
      deallocate ( Atm(n)%ze0   )
    end do


 end subroutine fv_end
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
!     run_setup :: initialize run from namelist
!
      subroutine run_setup(Atm, dt_atmos)
      type(fv_atmos_type), intent(inout) :: Atm
      real, intent(in)                   :: dt_atmos

      character(len=80) :: filename, tracerName
      integer :: ios, f_unit
      logical :: exists

      real :: dim0 = 180.           ! base dimension
      real :: dt0  = 1800.          ! base time step
      real :: ns0  = 5.             ! base nsplit for base dimension 
                                    ! For cubed sphere 5 is better
      real :: umax = 350.           ! max wave speed for grid_type>3
      real :: dimx, dl, dp, dxmin, dymin, d_fac

      integer :: n0split
      integer :: unit

      namelist /mpi_nml/npes_x,npes_y  ! Use of this namelist is deprecated
      namelist /fv_grid_nml/grid_name,grid_file
      namelist /fv_core_nml/npx, npy, ntiles, npz, uniform_vert_spacing, npz_rst, layout, io_layout, ncnst, nwat, &
                            k_split, n_split, m_split, q_split, print_freq,   &
                            hord_mt, hord_vt, hord_tm, hord_dp, hord_ze, hord_tr, &
                            kord_mt, kord_tm, kord_tr, fv_debug, fv_land, nudge,  &
                            external_ic, ncep_ic, fv_diag_ic, res_latlon_dynamics, res_latlon_tracers, &
                            dddmp, d2_bg, d4_bg, vtdm4, d_ext, beta1, beta, non_ortho, n_sponge, &
                            warm_start, adjust_dry_mass, mountain, d_con, nord, no_cgrid, init_wind_m, &
                            dry_mass, grid_type, do_Held_Suarez, consv_te, fill, tq_filter, filter_phys, fill_dp, &
                            range_warn, dwind_2d, inline_q, z_tracer, reproduce_sum, adiabatic,        &
                            tau, tau_h2o, rf_center, rf_lowest, nf_omega, hydrostatic, fv_sg_adj, breed_vortex_inline,  &
                            hybrid_z, quick_p_c, quick_p_d, Make_NH, m_grad_p,   &
                            a2b_ord, uniform_ppm, remap_t, k_top, m_riem, p_ref, &
#ifdef MARS_GCM
                            sponge_damp, reference_sfc_pres,                     &
#endif
                            c2l_ord, dx_const, dy_const, umax, deglat,     &
                            deglon_start, deglon_stop, deglat_start, deglat_stop, &
                            phys_hydrostatic, make_hybrid_z, ppm_limiter, old_divg_damp, &
                            debug_message_size, write_grid_char_file, non_dilute


      namelist /test_case_nml/test_case,alpha

! Make alpha = 0 the default:
      alpha = 0.
      test_case = 11   ! (USGS terrain)

      filename = "input.nml"
      inquire(file=filename,exist=exists)

      if (.not. exists) then  ! This will be replaced with fv_error wrapper
        if(master) write(6,*) "file ",trim(filename)," doesn't exist" 
        call mpp_error(FATAL,'FV core terminating')
      else

        f_unit = get_unit()
        open (f_unit,file=filename)
 ! Read Main namelist
        rewind (f_unit)
        read (f_unit,fv_grid_nml,iostat=ios)
        if (ios .gt. 0) then
           if(master) write(6,*) 'fv_grid_nml ERROR: reading ',trim(filename),', iostat=',ios
           call mpp_error(FATAL,'FV core terminating')
        endif
        unit = stdlog()
        write(unit, nml=fv_grid_nml)

 ! Read FVCORE namelist 
        rewind (f_unit)
        read (f_unit,fv_core_nml,iostat=ios)
        if (ios .ne. 0) then
           if(master) write(6,*) 'fv_core_nml ERROR: reading ',trim(filename),', iostat=',ios
           call mpp_error(FATAL,'FV core terminating')
        endif
        unit = stdlog()
        write(unit, nml=fv_core_nml)
        
!*** single tile for Cartesian grids
        if (grid_type>3) then
           ntiles=1
           non_ortho = .false.
           uniform_ppm = .true.
           nf_omega = 0
        endif

!*** remap_t is NOT yet supported for non-hydrostatic option *****
        if ( .not. hydrostatic ) remap_t = .false.

        npes_x = layout(1)
        npes_y = layout(2)
        io_domain_layout = io_layout

 ! Read Test_Case namelist
        rewind (f_unit)
        read (f_unit,test_case_nml,iostat=ios)
        if (ios .gt. 0) then
         if(master) write(6,*) 'test_case_nml ERROR: reading ',trim(filename),', iostat=',ios
            call mpp_error(FATAL,'FV core terminating')
        endif
        unit = stdlog()
        write(unit, nml=test_case_nml)

 ! Look for deprecated mpi_nml
        rewind (f_unit)
        read (f_unit,mpi_nml,iostat=ios)
        if (ios == 0) then
           call mpp_error(FATAL,'mpi_nml is deprecated. Use layout in fv_core_nml')
        endif

        close (f_unit)
      endif

! Define n_split if not in namelist
          if (ntiles==6) then
             dimx = 4.0*(npx-1)
#ifdef MARS_GCM
             ns0 = 8
#else
             if ( hydrostatic ) then
                  if ( npx >= 120 ) ns0 = 6
             else
                  if ( npx <= 45 ) then
                       ns0 = 6
                  elseif ( npx <=90 ) then
                       ns0 = 7
                  else
                       ns0 = 8
                  endif
             endif
#endif
          else
             dimx = max ( npx, 2*(npy-1) )
          endif
          
          if (grid_type < 4) then
             n0split = nint ( ns0*abs(dt_atmos)*dimx/(dt0*dim0) + 0.49 )
          elseif (grid_type == 4 .or. grid_type == 7) then
             n0split = nint ( 2.*umax*dt_atmos/sqrt(dx_const**2 + dy_const**2) + 0.49 )
          elseif (grid_type == 5 .or. grid_type == 6) then
             if (grid_type == 6) then
                deglon_start = 0.; deglon_stop  = 360.
             endif
             dl = (deglon_stop-deglon_start)*pi/(180.*(npx-1))
             dp = (deglat_stop-deglat_start)*pi/(180.*(npy-1))

             dxmin=dl*radius*min(cos(deglat_start*pi/180.-ng*dp),   &
                                 cos(deglat_stop *pi/180.+ng*dp))
             dymin=dp*radius
             n0split = nint ( 2.*umax*dt_atmos/sqrt(dxmin**2 + dymin**2) + 0.49 )
          endif
          n0split = max ( 1, n0split )

      if ( n_split == 0 ) then
           n_split = n0split/k_split
           if(master) write(6,*) 'For k_split (remapping)=', k_split
           if(master) write(6,198) 'n_split is set to ', n_split, ' for resolution-dt=',npx,npy,ntiles,dt_atmos
      else
          if(master) write(6,199) 'Using n_split from the namelist: ', n_split
      endif

!----------------------------------------
! Adjust divergence damping coefficients:
!----------------------------------------
!      d_fac = real(n0split)/real(n_split)
!      dddmp = dddmp * d_fac
!      d2_bg = d2_bg * d_fac
!      d4_bg = d4_bg * d_fac
!      d_ext = d_ext * d_fac
!      vtdm4 = vtdm4 * d_fac
      if (old_divg_damp) then
        if (master) write(6,*) " fv_control: using original values for divergence damping "
        d2_bg_k1 = 6.         ! factor for d2_bg (k=1)  - default(4.)
        d2_bg_k2 = 4.         ! factor for d2_bg (k=2)  - default(2.)
        d2_divg_max_k1 = 0.02 ! d2_divg max value (k=1) - default(0.05)
        d2_divg_max_k2 = 0.01 ! d2_divg max value (k=2) - default(0.02)
        damp_k_k1 = 0.        ! damp_k value (k=1)      - default(0.05)
        damp_k_k2 = 0.        ! damp_k value (k=2)      - default(0.025)
      endif

      if ( (.not.hydrostatic) .and. (m_split==0) ) then
           m_split = max(1., 0.5 + abs(dt_atmos)/(n_split*6.) )
           if(master) write(*,198) 'm_split is set to ', m_split
      endif

      if(master) then
         write(6,199) 'Using n_sponge : ', n_sponge
         write(6,197) 'Using non_ortho : ', non_ortho
      endif

      if ( hydrostatic ) m_grad_p = 1

 197  format(A,l7)
 198  format(A,i2.2,A,i4.4,'x',i4.4,'x',i1.1,'-',f9.3)
 199  format(A,i2.2)
! 200  format(A,A,i4.4,A,i4.4,A)
! 201  format(A,A,f5.3,A,i4.4,A,i4.4,A)
! 202  format(A,A,A,i4.4,A,i4.4,A)
! 210  format(A,A,f5.3,A,i4.4,A,i4.4,A,i2.2,A)

      alpha = alpha*pi

      call domain_decomp(npx,npy,ntiles,ng,grid_type)

  end subroutine run_setup


end module fv_control_mod
