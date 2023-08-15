module exo_phys_mod

   !-----------------------------------------------------------------------
   ! Paths:
   ! constants_mod:      src/shared/constants/constants.F90
   ! mpp_update_domains: src/shared/mpp/mpp_domains.F90
   ! time_manager_mod:   src/shared/time_manager/time_manager.F90
   ! diag_manager_mod:   src/shared/diag_manager/diag_manager.F90
   ! fv_grid_tools_mod:  src/atmos_cubed_sphere/tools/fv_grid_tools.F90
   ! fv_grid_utils_mod:  src/atmos_cubed_sphere/model/fv_grid_utils.F90
   ! fv_mp_mod:          src/atmos_cubed_sphere/tools/fv_mp_mod.F90
   ! fv_diagnostics_mod: src/atmos_cubed_sphere/tools/fv_diagnostics.F90
   ! fv_timing_mod:      src/atmos_cubed_sphere/tools/fv_timing.F90
   !
   ! Exo_Tend arguments:
   ! (npx,npy,npz): number of gridpoints in the (latitudes,longitudes,altitudes) grid
   ! is, ie: first and last indices of latitudes over which the model runs
   ! js, je: first and last indices of longitudes over which the model runs
   ! ng: ?
   ! nq: ?
   ! (u,v,w): wind speeds in the (latitudes,longitudes,altitudes) grid [m/s]
   ! pt: temperature array on midpoint levels [K]
   ! q: Moisture
   ! ts: surface temperature [K]
   ! pe: pressure array on interface levels (a.k.a. edges) [Pa]
   ! delp: pressure difference array (delta P) [Pa]
   ! peln: ln(pe) 
   ! pkz: Finite volume average of p**kappa
   ! pdt: pressure tendency? [Pa/s ?]
   ! (ua,va): ?
   ! (u_dt,v_dt): wind speed tendencies in the (latitudes, longitudes) grid [m/s-2 ?]
   ! t_dt: temperature tendency [s-1]
   ! t_dt_conv: temperature tendency due to radiation
   ! ts_dt: surface temperature tendency
   ! agrid: latitudes and longitudes storage array
   ! delz: height difference array [m]
   ! hydrostatic: switch for hydrostatic approximation
   ! (ak,bk): pressure array coefficients 
   ! ks
   ! strat
   ! rayf
   ! master
   ! Time
   !-----------------------------------------------------------------------

  use constants_mod,       only: grav, rdgas, rvgas, RADIAN, kappa, radius, cp_air, rho_cp, &
                                 cp_vapor
   use mpp_domains_mod,     only: mpp_update_domains
   use mpp_mod,             only: FATAL, mpp_error, mpp_pe, stdlog, &
      mpp_npes, mpp_get_current_pelist, get_unit
   use time_manager_mod,    only: time_type, get_date, get_time
   use diag_manager_mod,    only: send_data, register_diag_field
   use tracer_manager_mod,  only: get_tracer_index
   use field_manager_mod,   only: MODEL_ATMOS
   use fv_grid_utils_mod,   only: g_sum
   use fv_diagnostics_mod,  only: prt_maxmin, get_height_field
   use fv_timing_mod,       only: timing_on, timing_off
   use relax_mod,           only: relax_tend
   use fv_mp_mod,           only : is_master
   use simple_boundary_mod, only: simple_boundary_tend
   use dry_adj_mod,         only: dry_adj
   use dry_convection_mod,         only: dry_convection
   use moist_adj_mod,       only: moist_adj
   use rad_coupler_mod,     only: rad_coupler
   use condense_mod,        only: cond_tend
   use netcdf_send_mod,     only: send_data_to_netcdf, send_2D_data_to_netcdf
   use netcdf_reg_mod,      only: id_net_F, id_net_Fs, id_olr, id_opr_IR, id_opr_W1, id_opr_W2, id_opr_UV,&
                                  id_opr_VIS1, id_opr_VIS2, id_cff, id_scff, id_direct_down, id_t_dt_rad,&
                                  id_t_dt_conv, id_t_dt_conv_moist, id_height, id_isr, id_surf_sw_down,&
                                  id_surf_lw_down, id_flux_t


   use         surface_flux_mod, only: surface_flux
   use     vert_turb_driver_mod, only: vert_turb_driver_init,                    &
      vert_turb_driver, vert_turb_driver_end
   use            vert_diff_mod, only: gcm_vert_diff_init, gcm_vert_diff_down, &
      gcm_vert_diff_up, gcm_vert_diff_end, &
      surf_diff_type, gcm_vert_diff
   use          mixed_layer_mod, only: mixed_layer, mixed_layer_init
   use FMS_dry_adj_Ray, only: Ray_dry_adj

   use ding_convection, only: ding_adjust, rain_out_simple, large_scale_cond,&
                              ding_adjust_init, rain_out_revap, calc_enthalpy

   implicit none

   logical :: rf_initialized = .false.

   private
   public :: Exo_Tend, Exo_Tend_init, do_ding_convection

   type(surf_diff_type) :: Tri_surf

   namelist /exo_phys_nml/ relax_module, surface_on, tidally_locked, cp_surf, do_dry_adjustment, do_moist_H2O_adjustment, &
           do_virtual, do_simple_bl, do_vert_diff_bl, do_dry_convection, do_condensation, do_dry_adj_el, conv_timescale, do_ding_convection


   logical :: relax_module = .false.
   logical :: do_dry_adjustment = .true.
   logical :: do_dry_adj_el = .false.
   logical :: do_dry_convection = .false.
   logical :: do_moist_H2O_adjustment = .false.
   logical :: surface_on = .true.
   logical :: do_condensation = .false.
   logical :: tidally_locked = .true.
   real    :: cp_surf = 1.e6
   real    :: conv_timescale = 1.e3
   logical :: do_virtual = .false. ! whether virtual temp used in gcm_vert_diff
   logical :: do_simple_bl = .true.
   logical :: do_vert_diff_bl = .false.
   logical :: do_ding_convection = .true.





contains

   subroutine Exo_Tend_init(is,ie,js,je,npz,nq,Time,axes)
      integer, intent(in) :: is,ie,js,je,npz,nq
      integer :: f_unit, unit, ios ! Namelist reading
      character(80) :: filename
      type(time_type), intent(in) :: Time
      integer :: axes(4)

      filename = "input.nml"
      f_unit = get_unit()
      open (f_unit,file=filename)
      rewind (f_unit)
      read (f_unit,exo_phys_nml ,iostat=ios)
      unit = stdlog()
      write(unit, nml=exo_phys_nml)
      
      close(f_unit)
      if (do_vert_diff_bl) then
         call vert_turb_driver_init (ie-is+1,je-js+1,npz,Time)
         call gcm_vert_diff_init (Tri_surf, is,ie, js,je, npz, .true., do_virtual, .true.)
         call mixed_layer_init(is, ie, js, je, npz, axes, Time)
      end if

      if (do_ding_convection) call ding_adjust_init

   end subroutine Exo_Tend_init


   subroutine Exo_Tend(npx, npy, npz, is, ie, js, je, ng, nq, u, v, w, pt, q, ts, pe, delp, peln, pkz, dt_atmos, ua, va, u_dt, v_dt,&
      t_dt, ts_dt, q_dt, agrid, delz, hydrostatic, ak, bk, ks, strat, rayf, master, non_dilute, Time)

      integer, intent(IN   ) :: npx, npy, npz
      integer, intent(IN   ) :: is, ie, js, je, ng, nq
      logical, intent(IN)    :: hydrostatic, non_dilute
      real   , intent(IN) ::  delz(is:ie,js:je,npz)
      real   , intent(IN) ::  dt_atmos

      real   , intent(INOUT) ::    u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
      real   , intent(INOUT) ::    v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
      real   , intent(INOUT) ::    w(is-ng:ie+1+ng,js-ng:je+  ng,npz)
      real   , intent(INOUT) ::   pt(is-ng:ie+  ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) ::   ts(is-ng:ie+  ng,js-ng:je+  ng)
      real   , intent(INOUT) :: delp(is-ng:ie+  ng,js-ng:je+  ng,npz)
      real   , intent(INOUT) ::    q(is-ng:ie+  ng,js-ng:je+  ng,npz, nq)
      real   , intent(INOUT) ::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
      real   , intent(INOUT) :: peln(is  :ie     ,1:npz+1,js   :je     )
      real   , intent(INOUT) ::  pkz(is  :ie     ,js   :je     ,1:npz)

      real    ::   p_mid(is-1:ie+1 ,js-1:je+1,1:npz)
      real    ::   p_edge(is-1:ie+1 ,js-1:je+1,1:npz+1)
      real    ::   log_pe(is-1:ie+1 ,1:npz+1,js-1:je+1)

      real   , intent(INOUT) ::   ua(is-ng:ie+ng,js-ng:je+ng,npz)
      real   , intent(INOUT) ::   va(is-ng:ie+ng,js-ng:je+ng,npz)

      real    ::   dtau_du(is:ie,js:je)
      real    ::   dtau_dv(is:ie,js:je)
      real    ::   dsens_datmos(is:ie,js:je)
      real    ::   devap_datmos(is:ie,js:je)
      real,    dimension(is:ie,js:je)       :: tau_u, tau_v, sens, evap
      real,    dimension(is:ie,js:je,1:npz)     :: dt_u, dt_v, dt_t, dt_q
      real,    dimension(is:ie,js:je,1:npz,nq)   :: dt_tr
      real,    dimension(is:ie,js:je,1:npz)     :: dissipative_heat


      ! Tendencies:
      real, intent(INOUT):: u_dt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, intent(INOUT):: v_dt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, intent(INOUT):: t_dt(is:ie,js:je,npz)
!    real :: t_dt_rad(is:ie,js:je,npz)
      real :: u_dt_diff(is:ie,js:je,npz)
      real :: v_dt_diff(is:ie,js:je,npz)
      real :: t_dt_diff(is:ie,js:je,npz)
      real :: q_dt_diff(is:ie,js:je,npz,nq)
      real :: tr_dt_diff(is:ie,js:je,npz,nq)

      real :: t_dt_conv_1D(npz)
      real :: t_dt_conv(is:ie,js:je,npz)
      real :: t_dt_conv_moist(is:ie,js:je,npz)
      real :: t_dt_conv_ding(is:ie,js:je,npz)
      real :: q_dt_ding(is:ie,js:je,npz,nq), q_dt_ding_test(is:ie,js:je,npz,nq)
      real :: t_dt_lsc(is:ie,js:je,npz), t_dt_rainout(is:ie,js:je,npz)
      
      ! Create these arrays locally because will be rained out immediately
      real :: q_ice(is:ie,js:je,npz)
      real :: q_liq(is:ie,js:je,npz)
      
      real, intent(INOUT):: ts_dt(is:ie,js:je)
      real, INTENT(INOUT):: q_dt(is:ie,js:je,npz,nq)

      real   , intent(IN   ) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
      real   , intent(IN   ) :: ak(npz+1), bk(npz+1)
      integer, intent(IN   ) :: ks

      !real   , intent(IN   ) :: pdt
      logical, intent(IN   ) :: strat, rayf, master

      type(time_type), intent(in) :: Time
      
      ! Radiation
      real   :: net_F(is:ie,js:je,1:npz+1)
      real   :: olr(is:ie,js:je), opr_IR(is:ie,js:je), opr_W1(is:ie,js:je), opr_W2(is:ie,js:je), opr_UV(is:ie,js:je),&
                opr_VIS1(is:ie,js:je), opr_VIS2(is:ie,js:je)
      real   :: net_Fs(is:ie,js:je)
      real   :: cff(is:ie,js:je,1:npz)
      real   :: scff(is:ie,js:je)
      real   :: direct_down(is:ie,js:je,1:npz+1)
      real   :: surf_lw_down(is:ie,js:je)
      real   :: surf_sw_down(is:ie,js:je)
      integer :: i,j,k
      real :: PI=4.D0*DATAN(1.D0)

      real, allocatable, dimension(:,:,:)   :: t_dt_rad

      real    ::   height(is:ie,js:je,npz+1)
      real    ::   height_mid(is:ie,js:je,npz)

      real  :: h_col(is:ie, js:je), h_col_new(is:ie,js:je)
      real  :: m_col(is:ie, js:je), m_col_new(is:ie,js:je)
      real, dimension(is:ie, js:je, 1:npz) :: q_tmp, t_tmp
      real :: htmp, mtmp
! ---
      logical, dimension(is:ie,js:je) ::                                       &
         avail,                &   ! generate surf. flux (all true)
         land,                 &   ! land points (all false)
         coldT                     ! should precipitation be snow at this point

      real, dimension(is:ie,js:je,1:npz) ::                                        &
         diff_m,               &   ! momentum diffusion coeff.
         diff_t,               &   ! temperature diffusion coeff.
         diss_heat,            &   !  heat dissipated by vertical diffusion
         non_diff_dt_ug,       &   ! zonal wind tendency except from vertical diffusion
         non_diff_dt_vg,       &   ! merid. wind tendency except from vertical diffusion
         non_diff_dt_tg,       &   ! temperature tendency except from vertical diffusion
         non_diff_dt_qg,       &   ! moisture tendency except from vertical diffusion
         conv_dt_tg,           &   ! temperature tendency from convection
         conv_dt_qg,           &   ! moisture tendency from convection
         cond_dt_tg,           &   ! temperature tendency from condensation
         cond_dt_qg                ! moisture tendency from condensation

      real, dimension(is:ie,js:je,1:nq) ::                                        &
         flux_tr                 ! surface tracer flux

      real, dimension(is:ie,js:je)   ::                       &
         z_surf,               &   ! surface height
         t_surf,               &   ! surface temperature
         q_surf,               &   ! surface moisture
         u_surf,               &   ! surface U wind
         v_surf,               &   ! surface V wind
         rough_mom,            &   ! momentum roughness length for surface_flux
         rough_heat,           &   ! heat roughness length for surface_flux
         rough_moist,          &   ! moisture roughness length for surface_flux
         depth_change_lh,      &   ! tendency in bucket depth due to latent heat transfer     ! added by LJJ
         depth_change_cond,    &   ! tendency in bucket depth due to condensation rain        ! added by LJJ
         depth_change_conv,    &   ! tendency in bucket depth due to convection rain          ! added by LJJ
         gust,                 &   ! gustiness constant
         flux_t,               &   ! surface sensible heat flux
         flux_q,               &   ! surface moisture flux
         flux_r,               &   ! surface radiation flux
         flux_u,               &   ! surface flux of zonal mom.
         flux_v,               &   ! surface flux of meridional mom.
         drag_m,               &   ! momentum drag coefficient
         drag_t,               &   ! heat drag coefficient
         drag_q,               &   ! moisture drag coefficient
         w_atm,                &   ! wind speed
         ustar,                &   ! friction velocity
         bstar,                &   ! buoyancy scale
         qstar,                &   ! moisture scale
         dhdt_surf,            &   ! d(sensible heat flux)/d(surface temp)
         dedt_surf,            &   ! d(latent heat flux)/d(surface temp)???
         dedq_surf,            &   ! d(latent heat flux)/d(surface moisture)???
         drdt_surf,            &   ! d(upward longwave)/d(surface temp)
         dhdt_atm,             &   ! d(sensible heat flux)/d(atmos.temp)
         dedq_atm,             &   ! d(latent heat flux)/d(atmospheric mixing rat.)
         dtaudv_atm,           &   ! d(stress component)/d(atmos wind)
         fracland,             &   ! fraction of land in gridbox
         rough                     ! roughness for vert_turb_driver

      real, dimension(is:ie,js:je,2)   :: bucket_depth                ! added by LJJ
      real, dimension(is:ie,js:je) :: dt_psg, dt_bucket, bucket_diffusion, filt   ! added by LJJ

      ! Nondiluteness
      real, dimension(is:ie,js:je,npz) :: cp_local
      logical :: nondilute = .true.
      logical :: stop_switch = .false.
      integer :: vap, cond

      !-----------------------------------------------------------------------------------------------
      !                                         ALLOCATIONS
      !-----------------------------------------------------------------------------------------------
      allocate(t_dt_rad   (is:ie,js:je,npz))

      t_dt(:,:,:) = 0.0
      ts_dt(:,:) = 0.0

      t_dt_diff(:,:,:) = 0.0
      u_dt_diff(:,:,:) = 0.0
      v_dt_diff(:,:,:) = 0.0


      ! CHECK FOR NANS
      do k=1,npz
         do j=js,je
            do i=is,ie
               if (pt(i,j,k) /= pt(i,j,k)) then
                  stop_switch=.true.
                  write(*,*) 'NaN in temperature array @ i,j,k', i,j,k
                  stop
               endif
               if (u(i,j,k) /= u(i,j,k)) then
                  stop_switch=.true.
                  write(*,*) 'NaN in u array @ i,j,k ', i,j,k
                  stop
               endif

               if (v(i,j,k) /= v(i,j,k)) then
                  stop_switch=.true.
                  write(*,*) 'NaN in v array @ i,j,k ', i,j,k
                  stop
               endif
            enddo
         enddo
      enddo

      
      !-----------------------------------------------------------------------------------------------
      !                                   LOCAL ARRAYS DEFINITIONS
      !-----------------------------------------------------------------------------------------------

      !-----------------------------------------------------------------------------------------------
      !                                          RELAXATION
      !-----------------------------------------------------------------------------------------------

      if (relax_module) then
         call relax_tend(npx, npy, npz, is, ie, js, je, ng, nq,   &
            pt, pe, delp, peln, pkz,   &
            t_dt(is:ie,js:je,1:npz), agrid,  &
            delz, hydrostatic,   &
            master, Time)
      end if


      !-----------------------------------------------------------------------------------------------
      !                                          RADIATION
      !-----------------------------------------------------------------------------------------------


      t_dt_rad(is:ie,js:je,1:npz) = 0. ! Ensures t_dt_rad = 0 if radiation is turned off

      do i = is,ie
         do j = js,je
            do k = 1, npz
               p_mid(i,j,k) = delp(i,j,k) / log(pe(i,k+1,j)/pe(i,k,j))
            end do

            do k = 1, npz+1
               p_edge(i,j,k) = pe(i,k,j)
            end do
         end do
      end do

         call rad_coupler(surface_on, tidally_locked, is, ie, js, je, npz, ng, ts(is:ie,js:je), pt(is:ie,js:je,1:npz), p_mid(is:ie,js:je,1:npz), &
            p_edge(is:ie,js:je,1:npz+1), agrid(is:ie,js:je,1:2), net_F, olr, opr_IR, opr_W1, opr_W2, opr_UV, opr_VIS1, opr_VIS2, net_Fs,&
            surf_lw_down(is:ie,js:je), surf_sw_down(is:ie,js:je), cff, scff, direct_down)


      if (non_dilute) then
         vap = get_tracer_index(MODEL_ATMOS, 'vapour')
         do j=js,je
            do i=is,ie
               do k=1,npz
                  cp_local(i,j,k) = cp_air*(1-q(i,j,k,vap)) + cp_vapor*q(i,j,k,vap)
                  t_dt_rad(i,j,k) = grav/cp_local(i,j,k) * (net_F(i,j,k+1) - net_F(i,j,k))/delp(i,j,k)
               enddo
            enddo
         enddo
      else
         do i = is,ie
            do j = js,je
               do k = 1, npz
                  t_dt_rad(i,j,k) = (grav/cp_air) * (net_F(i,j,k+1)-net_F(i,j,k))/delp(i,j,k)
               end do
            end do
         end do
      endif
      
      !-----------------------------------------------------------------------------------------------
      !                                    BOTTOM BOUNDARY EFFECTS
      !-----------------------------------------------------------------------------------------------

      ! Choice between a surface and no surface as a logical input in a namelist?
      ! Physical constants used:
      !    - "stefan" as the Stefan-Boltzman constant
      !    - "cp_air" as specific heat capacity of dry air at constant pressure
      !    - "rdgas" as the specific gas constant of dry air
      ! Previously, the undefined variable "mld" was hardcoded here. Define it and put it in a namelist
      ! if appropriate.

         ! calculate log edge pressure
         do i = is,ie
            do j = js,je
               do k = 1, npz+1
                  log_pe(i,k,j) =  log(pe(i,k,j))
               end do
            end do
         end do

         ! calculate height field
         call get_height_field(is, ie, js, je, ng, npz, hydrostatic, delz, height(is:ie,js:je,1:npz+1),pt(is-ng:ie+ng,js-ng:je+ng,1:npz), q,log_pe(is:ie,1:npz+1,js:je), 0.0)

      if (do_simple_bl) then
         do i = is,ie
            do j = js,je
               ts_dt(i,j) = ts_dt(i,j) + net_Fs(i,j) / (cp_surf)
            end do
         end do

         call simple_boundary_tend(npx, npy, npz, is, ie, js, je, ng, nq,   &
            pt, ts, ua, va, pe, delp, peln, pkz,   &
            t_dt(is:ie,js:je,1:npz), ts_dt(is:ie,js:je), u_dt(is:ie,js:je,1:npz), v_dt(is:ie,js:je,1:npz), agrid,  &
            delz, hydrostatic,   &
            master, Time)


      else if (do_vert_diff_bl) then
         
         ! calculate height of midpoints
         do i = is,ie
            do j = js,je
               do k = 1, npz
                  height_mid(i,j,k) = 0.5*(height(i,j,k+1)+height(i,j,k))
               end do
            end do
         end do

         ! zero solid surface winds and q
         u_surf(is:ie,js:je) = 0.0
         v_surf(is:ie,js:je) = 0.0
         q_surf(is:ie,js:je) = 0.0
         ! don't use bucket
         bucket_depth = 0.0

         ! momentum roughness
         rough_mom(is:ie,js:je) = 3.21e-05
         ! temperature roughness
         rough_heat(is:ie,js:je) = 3.21e-05
         ! q roughness
         rough_moist(is:ie,js:je) = 3.21e-05
         ! gustiness coefficient
         gust(is:ie,js:je) = 0.0

         ! zero "land"
         land(is:ie,js:je) = .false.
         ! calculate surface flux everywhere
         avail(is:ie,js:je) = .true.
         ! zero "land" fraction
         fracland(is:ie,js:je) = 0.0

         ! surface flux routine
         call surface_flux(                        &
            pt(is:ie,js:je,npz),                   &    ! lowest level T
            q(is:ie,js:je,npz,nq),                 &    ! lowest level tracers
            ua(is:ie,js:je,npz),                   &    ! lowest level ua
            va(is:ie,js:je,npz),                   &    ! lowest level va
            p_mid(is:ie,js:je,npz),                &
            height(is:ie,js:je,npz),               &
            p_edge(is:ie,js:je,npz+1),             &
            ts(is:ie,js:je),                       &    ! not land surface T
            ts(is:ie,js:je),                       &    ! land surface T
            q_surf(is:ie,js:je),                   &
            bucket_depth(is:ie,js:je,1),           &    ! not used by default
            depth_change_lh(is:ie,js:je),          &    ! not used by default
            depth_change_conv(is:ie,js:je),        &    ! not used by default
            depth_change_cond(is:ie,js:je),        &    ! not used by default
            u_surf(is:ie,js:je),                   &
            v_surf(is:ie,js:je),                   &
            rough_mom(is:ie,js:je),                &
            rough_heat(is:ie,js:je),               &
            rough_moist(is:ie,js:je),              &
            gust(is:ie,js:je),                     &
            flux_t(is:ie,js:je),                   &
            flux_q(is:ie,js:je),                   &
            flux_r(is:ie,js:je),                   &
            flux_u(is:ie,js:je),                   &
            flux_v(is:ie,js:je),                   &
            drag_m(is:ie,js:je),                   &
            drag_t(is:ie,js:je),                   &
            drag_q(is:ie,js:je),                   &
            w_atm(is:ie,js:je),                    &
            ustar(is:ie,js:je),                    &
            bstar(is:ie,js:je),                    &
            qstar(is:ie,js:je),                    &
            dhdt_surf(is:ie,js:je),                &
            dedt_surf(is:ie,js:je),                &
            dedq_surf(is:ie,js:je),                &
            drdt_surf(is:ie,js:je),                &
            dhdt_atm(is:ie,js:je),                 &
            dedq_atm(is:ie,js:je),                 &
            dtaudv_atm(is:ie,js:je),               &
            dt_atmos,                              &
            land(is:ie,js:je),                     &
            avail(is:ie,js:je)  )


         ! vertical turbulence routine to calculate diffusivity coefficients
         rough(is:ie,js:je) = 0.0 ! M-Y roughness (not used)
         call vert_turb_driver(is, js, Time, Time, & ! only need correct times for M-Y (not used)
            dt_atmos, fracland(is:ie,js:je), p_edge(is:ie,js:je,1:npz+1), p_mid(is:ie,js:je,1:npz),            &
            height(is:ie,js:je,1:npz+1), height_mid(is:ie,js:je,1:npz), ustar(is:ie,js:je),bstar(is:ie,js:je), &
            rough(is:ie,js:je), ua(is:ie,js:je,1:npz), va(is:ie,js:je,1:npz), pt(is:ie,js:je,1:npz),           &
            q(is:ie,js:je,1:npz,nq), ua(is:ie,js:je,1:npz),                                                    & 
            va(is:ie,js:je,1:npz), pt(is:ie,js:je,1:npz),                                                      &
            q(is:ie,js:je,1:npz,nq), u_dt(is:ie,js:je,1:npz),                                                  &
            v_dt(is:ie,js:je,1:npz), t_dt(is:ie,js:je,1:npz),                                                  &
            q_dt(is:ie,js:je,1:npz,nq), diff_t(is:ie,js:je,1:npz),                                             &
            diff_m(is:ie,js:je,1:npz), gust(is:ie,js:je)  )


         ! vertical diffusion downward sweep
         call gcm_vert_diff_down (is, js,                                          &
            dt_atmos,             ua(is:ie,js:je,1:npz),       &
            va(is:ie,js:je,1:npz),  pt(is:ie,js:je,1:npz),       &
            q(is:ie,js:je,1:npz,1),             &
            q(is:ie,js:je,1:npz,1:nq), diff_m(is:ie,js:je,1:npz), &
            diff_t(is:ie,js:je,1:npz),                  p_edge(is:ie,js:je,1:npz+1), &
            p_mid(is:ie,js:je,1:npz),                  height_mid(is:ie,js:je,1:npz), &
            flux_u(is:ie,js:je),                    flux_v(is:ie,js:je),   &
            dtaudv_atm(is:ie,js:je),                               &
            flux_tr(is:ie,js:je,1:nq),                                &
            u_dt_diff(is:ie,js:je,1:npz),                    v_dt_diff(is:ie,js:je,1:npz), &
            t_dt_diff(is:ie,js:je,1:npz),          q_dt_diff(is:ie,js:je,1:npz,1), &
            tr_dt_diff(is:ie,js:je,1:npz,1:nq),         diss_heat(is:ie,js:je,1:npz), &
            Tri_surf)


         ! update surface temperature
         call mixed_layer(                                                 &
            is,                                          &
            ie,                                          &
            js,                                          &
            je,                                          &
            Time,                                        &
            ts(is:ie,js:je),                                 &
            flux_t(is:ie,js:je),                                 &
            flux_q(is:ie,js:je),                                 &
            flux_r(is:ie,js:je),                                 &
            flux_u(is:ie,js:je),                                 &
            dt_atmos,                                 &
            surf_sw_down(is:ie,js:je),                                 &
            surf_lw_down(is:ie,js:je),                                 &
            Tri_surf,                                      &
            dhdt_surf(is:ie,js:je),                                 &
            dedt_surf(is:ie,js:je),                                 &
            dedq_surf(is:ie,js:je),                                 &
            drdt_surf(is:ie,js:je),                                 &
            dhdt_atm(is:ie,js:je),                                 &
            dedq_atm(is:ie,js:je))

         ! vertical diffusion upward sweep
         call gcm_vert_diff_up (is,js, dt_atmos, Tri_surf,  &
            t_dt_diff(is:ie,js:je,1:npz), q_dt_diff(is:ie,js:je,1:npz,1) )

      end if
      !-----------------------------------------------------------------------------------------------
      !                                   DRY CONVECTIVE ADJUSTMENT
      !-----------------------------------------------------------------------------------------------

      t_dt_conv(is:ie,js:je,:) = 0. ! Ensures t_dt_conv = 0 if do_dry_adjustment = .false.
      t_dt_conv_ding(is:ie,js:je,:) = 0.0
      t_dt_lsc(is:ie,js:je,:) = 0.0
      t_dt_rainout(is:ie,js:je,:) = 0.0
      q_dt_ding(is:ie,js:je,:,:)  = 0.0
      
      if (do_dry_adjustment) then

         call dry_adj(pt(is:ie,js:je,:),p_mid(is:ie,js:je,:),p_edge(is:ie,js:je,:),t_dt_conv(is:ie,js:je,:))
         t_dt_conv(is:ie,js:je,:) = t_dt_conv(is:ie,js:je,:) / conv_timescale

      else if (do_dry_convection) then

         call dry_convection(Time,pt(is:ie,js:je,:),p_mid(is:ie,js:je,:),&
                 p_edge(is:ie,js:je,:),t_dt_conv(is:ie,js:je,:))

      else if (do_dry_adj_el) then

      do i = is, ie
         do j = js, je
            call Ray_dry_adj(npz, npz+1, conv_timescale, kappa, pt(i,j,1:npz), p_mid(i,j,1:npz),&
                  p_edge(i,j,1:npz+1), t_dt_conv_1D(1:npz))
            t_dt_conv(i,j,1:npz) = t_dt_conv_1D(1:npz)
         end do
      end do

      end if



      !-----------------------------------------------------------------------------------------------
      !                                  MOIST CONVECTIVE ADJUSTMENT
      !-----------------------------------------------------------------------------------------------

      ! Returns the temperature tendency due to moist convection t_dt_conv_moist(is:ie,js:je,:).
      ! Todo: add in evaporative cooling

      !write(*,*) 'BEFORE MOIST CONV ADJUSTMENT'
      !write(*,*) maxval(q(is:ie,js:je,:,1)), maxval(q(is:ie,js:je,:,2)), maxval(q(is:ie,js:je,:,3))
      
      t_dt_conv_moist(is:ie,js:je,:) = 0. ! Ensures t_dt_conv_moist = 0 if do_moist_H2O_adjustment = .false.
      if (do_moist_H2O_adjustment) then
         ! Moist adjustment loop, returns the temperature tendency due to moist convection t_dt_conv_moist(is:ie,js:je,:).
         ! Valid only for a pure steam atmosphere.
         do i = is, ie
            do j = js, je
               call moist_adj(pt(i,j,:),p_mid(i,j,:),p_edge(i,j,:),t_dt_conv_moist(i,j,:))
            end do
         end do

      else if (do_ding_convection) then
         
         q_ice(is:ie,js:je,:) = 0.0
         q_liq(is:ie,js:je,:) = 0.0
         vap = get_tracer_index(MODEL_ATMOS, 'vapour')

!         ! Calculate initial h and m
!         h_col(is:ie,js:je) = 0.0
!         m_col(is:ie,js:je) = 0.0
!         Do i = is,ie
!           do j=js,je
!               do k=1,npz
!                  call calc_enthalpy(pt(i,j,k), q(i,j,k,vap), 0.0,0.0, htmp)
!                  h_col(i,j) = h_col(i,j) + htmp*delp(i,j,k)
!
!                  m_col(i,j) = m_col(i,j) + q(i,j,k,vap)*delp(i,j,k)
!               enddo
!            enddo
!         enddo
!         
         do i=is,ie
            do j=js,je
              call ding_adjust(p_mid(i,j,:), delp(i,j,:), pt(i,j,:), q(i,j,:,:), &
                                t_dt_conv_ding(i,j,:), q_dt_ding(i,j,:,:), q_liq(i,j,:), q_ice(i,j,:))
            enddo
         enddo
!         !write(*,*) 'POST DING CONVECTION'
!         ! Calculate initial h and m
!
!         ! Add these tendencies to moisture and reset, because we don't want dynamical core to change delp
!         ! as a result of these dq
         q(is:ie,js:je,1:npz,vap) = q(is:ie,js:je,1:npz,vap) + q_dt_ding(is:ie,js:je,1:npz, vap)
         q_dt_ding(is:ie,js:je,1:npz,vap) = 0.0

         ! Add negative here in anticipation of this being rained out
         cond = get_tracer_index(MODEL_ATMOS, 'condensate')
         q_dt_ding(is:ie,js:je,:,cond) = q_dt_ding(is:ie,js:je,:,cond) - q_liq(is:ie,js:je,:) -  q_ice(is:ie,js:je,:)
         q(is:ie,js:je,1:npz, cond) = -q_dt_ding(is:ie,js:je,:,cond)/dt_atmos + 1.e-10
!
!         !! =========================================================================================
!         !! Uncomment for mass and enthalpy diagnostics
         h_col_new(is:ie,js:je) = 0.0
         m_col_new(is:ie,js:je) = 0.0
         Do i = is,ie
            do j=js,je
               do k=1,npz
                  call calc_enthalpy(pt(i,j,k) + t_dt_conv_ding(i,j,k), q(i,j,k,vap) + q_dt_ding(i,j,k,vap), q_liq(i,j,k),&
                       q_ice(i,j,k), htmp)
                  h_col_new(i,j) = h_col_new(i,j) + htmp*delp(i,j,k)

                  m_col_new(i,j) = m_col_new(i,j) + &
                       (q(i,j,k,vap) + q_dt_ding(i,j,k,vap) + q_liq(i,j,k) + q_ice(i,j,k))*delp(i,j,k)
               enddo
            enddo
         enddo

         write(*,*) 'Ding convection mass and enthalpy change'
         write(*,*) 'Enthalpy: ', maxval(abs((h_col_new(is:ie,js:je) - h_col(is:ie,js:je))/h_col(is:ie,js:je)))
         write(*,*) 'H2O mass: ', maxval(abs((m_col_new(is:ie,js:je)-m_col(is:ie,js:je))/m_col(is:ie,js:je)))

         h_col(is:ie,js:je) = h_col_new(is:ie,js:je)
         m_col(is:ie,js:je) = m_col_new(is:ie,js:je)
!         !! =========================================================================================
!                     
!      !-----------------------------------------------------------------------------------------------
!      !                                    LARGE-SCALE CONDENSATION
!      !-----------------------------------------------------------------------------------------------
!
         do i=is,ie
            do j=js,je
               call large_scale_cond(p_mid(i,j,:), pt(i,j,:) + t_dt_conv_ding(i,j,:), &
                    q(i,j,:,:), q_liq(i,j,:), q_ice(i,j,:), &
                    q_dt_ding(i,j,:,:), t_dt_lsc(i,j,:))

         
            enddo
         enddo
!         !! =========================================================================================
!         !! Uncomment for enthalpy and mass diagnostics
!         h_col_new(is:ie,js:je) = 0.0
!         m_col_new(is:ie,js:je) = 0.0
!         Do i = is,ie
!            do j=js,je
!               do k=1,npz
!                  call calc_enthalpy(pt(i,j,k) + t_dt_conv_ding(i,j,k)+ t_dt_lsc(i,j,k), &
!                        q(i,j,k,vap) + q_dt_ding(i,j,k,vap), q_liq(i,j,k),&
!                       q_ice(i,j,k), htmp)
!                 h_col_new(i,j) = h_col_new(i,j) + htmp*delp(i,j,k)
!
!                  m_col_new(i,j) = m_col_new(i,j) + &
!                       (q(i,j,k,vap) + q_dt_ding(i,j,k,vap) + q_liq(i,j,k) + q_ice(i,j,k))*delp(i,j,k)
!               enddo
!            enddo
!        enddo
!
!         write(*,*) 'LSC mass and enthalpy change'
!         write(*,*) 'Enthalpy: ', maxval(abs((h_col_new(is:ie,js:je) - h_col(is:ie,js:je))/h_col(is:ie,js:je)))
!         write(*,*) 'H2O mass: ', maxval(abs((m_col_new(is:ie,js:je)-m_col(is:ie,js:je))/m_col(is:ie,js:je)))
!
!         h_col(is:ie,js:je) = h_col_new(is:ie,js:je)
!         m_col(is:ie,js:je) = m_col_new(is:ie,js:je)
! ====================================================================================================
!
         do i=is,ie
            do j=js,je
               do k=1,npz
                  call rain_out_revap(pt(i,j,:) + t_dt_conv_ding(i,j,:) + t_dt_lsc(i,j,:), &
                       p_mid(i,j,:), delp(i,j,:), q(i,j,:,vap), &
                       q_liq(i,j,:), q_ice(i,j,:), q_dt_ding(i,j,:,:), t_dt_rainout(i,j,:))

               enddo
            enddo
         enddo
!
!         !! =========================================================================================
!         !! Uncomment for enthalpy diagnostics
!         
!         h_col_new(is:ie,js:je) = 0.0
!         m_col_new(is:ie,js:je) = 0.0
!         Do i = is,ie
!            do j=js,je
!              do k=1,npz
!
!                  ! Account for the change in layer masses, especially for enthalpy !
!                  call calc_enthalpy(pt(i,j,k) + t_dt_conv_ding(i,j,k) + t_dt_lsc(i,j,k) + t_dt_rainout(i,j,k), &
!                       (q(i,j,k,vap) + q_dt_ding(i,j,k,vap))/(1. + q_dt_ding(i,j,k,vap) + q_dt_ding(i,j,k,cond)),&
!                       q_liq(i,j,k)/(1. + q_dt_ding(i,j,k,vap) + q_dt_ding(i,j,k,cond)), &
!                       q_ice(i,j,k)/(1. + q_dt_ding(i,j,k,vap) + q_dt_ding(i,j,k,cond)), &
!                       htmp)
!                  !call calc_enthalpy(pt(i,j,k) + t_dt_conv_ding(i,j,k), q(i,j,k,vap) + q_dt_ding(i,j,k,vap), q_liq(i,j,k),&
!                  !     q_ice(i,j,k), htmp)
!                  h_col_new(i,j) = h_col_new(i,j) + htmp*delp(i,j,k)*(1 + q_dt_ding(i,j,k,vap) + q_dt_ding(i,j,k,cond))
!
!                  ! Note we don't have to bother here because alterations to q and delp cancel out
!                  ! to conserve q*dp (but won't conserve enthalpy calculation because of dry enthalpy!)
!                  m_col_new(i,j) = m_col_new(i,j) + &
!                       (q(i,j,k,vap) + q_dt_ding(i,j,k,vap) + q_liq(i,j,k) + q_ice(i,j,k))*delp(i,j,k)
!               enddo
!            enddo
!         enddo
!
!         write(*,*) 'Rain out mass and enthalpy change'
!         write(*,*) 'Enthalpy: ', maxval(abs((h_col_new(is:ie,js:je) - h_col(is:ie,js:je))/h_col(is:ie,js:je)))
!         write(*,*) 'H2O mass: ', maxval(abs((m_col_new(is:ie,js:je)-m_col(is:ie,js:je))/m_col(is:ie,js:je)))
!         write(*,*) 'maxval qliq,qice', maxval(q_liq(is:ie,js:je,1:npz)), maxval(q_ice(is:ie,js:je,1:npz))

         !! =========================================================================================
         
      q_dt_ding(is:ie,js:je,1:npz,1:nq)= q_dt_ding(is:ie,js:je,1:npz,1:nq)/dt_atmos
      t_dt_conv_ding(is:ie,js:je,1:npz) = t_dt_conv_ding(is:ie,js:je,1:npz)/dt_atmos
      t_dt_lsc(is:ie,js:je,1:npz) = t_dt_lsc(is:ie,js:je,1:npz)/dt_atmos
      t_dt_rainout(is:ie,js:je,1:npz) = t_dt_rainout(is:ie,js:je,1:npz)/dt_atmos
   end if

      

      !-----------------------------------------------------------------------------------------------
      !                                            CLOUDS
      !-----------------------------------------------------------------------------------------------

      ! Not sure if this should be called in exo_Tend.
      ! Choice of scheme (SELECT CASE ?): user input should be in a namelist:
      ! - DIHRT
      ! - lin_cloud_microphys: GFDL global cloud resolving model
      ! - Maybe a more idealised/fast one?

      !-----------------------------------------------------------------------------------------------
      !                                   TOTAL TEMPERATURE TENDENCY
      !-----------------------------------------------------------------------------------------------

      ! t_dt adds up t_dt_rad + t_dt_conv + t_dt_conv_moist to output the total heating. Each component
      ! can be output to a netCDF file.


   
      t_dt(is:ie,js:je,1:npz) = t_dt(is:ie,js:je,1:npz)      &
         + t_dt_rad(is:ie,js:je,1:npz)  &
         + t_dt_conv(is:ie,js:je,1:npz) &
         + t_dt_conv_moist(is:ie,js:je,1:npz) & 
         + t_dt_diff(is:ie,js:je,1:npz)&
         + t_dt_conv_ding(is:ie,js:je,1:npz) &
         + t_dt_lsc(is:ie,js:je,1:npz) &
         + t_dt_rainout(is:ie,js:je,1:npz)

      u_dt(is:ie,js:je,1:npz) = u_dt(is:ie,js:je,:) + u_dt_diff(is:ie,js:je,1:npz)
      v_dt(is:ie,js:je,1:npz) = v_dt(is:ie,js:je,:) + v_dt_diff(is:ie,js:je,1:npz)

      q_dt(is:ie,js:je,1:npz,1:nq) = q_dt(is:ie,js:je,1:npz,1:nq) + q_dt_ding(is:ie,js:je,1:npz,1:nq)

!      t_dt(is:ie,js:je,1:npz)= 0.0
!      q_dt(is:ie,js:je,1:npz,1:nq) = 0.0
!      u_dt(is:ie,js:je,1:npz) = 0.0
!      v_dt(is:ie,js:je,1:npz)= 0.0
      
      !-----------------------------------------------------------------------------------------------
      !                                      SEND DATA TO NETCDF
      !-----------------------------------------------------------------------------------------------

      ! Uses the send_data scheme with the template associated to output X:

      !if ( id_X > 0 ) then
      !  used = send_data ( id_X, X, Time)
      !endif
      call send_data_to_netcdf(is,ie,js,je,npz,npz+1,Time,id_net_F,net_F)
      call send_2D_data_to_netcdf(is,ie,js,je,Time,id_olr,olr)
      call send_2D_data_to_netcdf(is,ie,js,je,Time,id_opr_IR,opr_IR)
      call send_2D_data_to_netcdf(is,ie,js,je,Time,id_opr_W1,opr_W1)
      call send_2D_data_to_netcdf(is,ie,js,je,Time,id_opr_W2,opr_W2)
      call send_2D_data_to_netcdf(is,ie,js,je,Time,id_opr_UV,opr_UV)
      call send_2D_data_to_netcdf(is,ie,js,je,Time,id_opr_VIS1,opr_VIS1)
      call send_2D_data_to_netcdf(is,ie,js,je,Time,id_opr_VIS2,opr_VIS2)
      call send_data_to_netcdf(is,ie,js,je,npz,npz,Time,id_cff,cff)
      call send_2D_data_to_netcdf(is,ie,js,je,Time,id_scff,scff)
      call send_data_to_netcdf(is,ie,js,je,npz,npz+1,Time,id_direct_down,direct_down)
      call send_data_to_netcdf(is,ie,js,je,npz,npz,Time,id_t_dt_rad,t_dt_rad)
      call send_data_to_netcdf(is,ie,js,je,npz,npz,Time,id_t_dt_conv,t_dt_conv)
      call send_data_to_netcdf(is,ie,js,je,npz,npz,Time,id_height,height(:,:,1:npz))
      call send_data_to_netcdf(is,ie,js,je,npz,npz,Time,id_t_dt_conv_moist,t_dt_conv_moist)

      call send_2D_data_to_netcdf(is,ie,js,je,Time,id_surf_sw_down,surf_sw_down)
      call send_2D_data_to_netcdf(is,ie,js,je,Time,id_surf_lw_down,surf_lw_down)
      call send_2D_data_to_netcdf(is,ie,js,je,Time,id_flux_t,flux_t)
      !-----------------------------------------------------------------------------------------------
      !                                         DEALLOCATIONS
      !-----------------------------------------------------------------------------------------------

      ! Deallocate here any array allocated above.

      !stop

   end subroutine Exo_Tend


end module exo_phys_mod

