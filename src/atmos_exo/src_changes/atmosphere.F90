!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module atmosphere_mod

!-----------------------------------------------------------------------
!
!    interface for FV dynamical core with Held-Suarez forcing
!
!-----------------------------------------------------------------------


use constants_mod, only: grav, kappa, cp_air, pi, rdgas, rvgas, SECONDS_PER_DAY
use fms_mod,       only: file_exist, open_namelist_file,   &
                         error_mesg, FATAL,                &
                         check_nml_error, stdlog, stdout,  &
                         write_version_number,             &
                         close_file, set_domain, nullify_domain, mpp_pe, mpp_root_pe
use time_manager_mod, only: time_type, get_time, set_time, operator(+)
use mpp_domains_mod,  only: domain2d
use mpp_mod,             only: FATAL, mpp_error, mpp_pe, stdlog, &
      mpp_npes, mpp_get_current_pelist, get_unit
!------------------
! FV specific codes:
!------------------
use fv_arrays_mod, only: fv_atmos_type
use fv_control_mod,only: fv_init, domain, fv_end, adiabatic, p_ref
use fv_phys_mod,   only: fv_phys, fv_nudge
use fv_diagnostics_mod, only: fv_diag_init, fv_diag, fv_time
use fv_timing_mod,   only: timing_on, timing_off
use fv_restart_mod, only: fv_restart
use fv_dynamics_mod, only: fv_dynamics
use fv_grid_tools_mod, only: grid_type
! MDH commented out
!use lin_cld_microphys_mod, only: lin_cld_microphys_init, lin_cld_microphys_end

! HII added
use exo_phys_init_mod, only : exo_init
use exo_phys_mod, only : Exo_Tend_init
!use          mixed_layer_mod, only: mixed_layer_init

!-----------------------------------------------------------------------

implicit none
private

public   atmosphere_init, atmosphere,  atmosphere_end, atmosphere_domain

namelist /atmosphere_nml/ do_dynamics, p_split, dry_atmosphere

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: atmosphere.F90,v 18.0 2010/03/02 23:26:58 fms Exp $'
character(len=128) :: tag = '$Name: riga $'

!-----------------------------------------------------------------------
!---- private data ----

type        (time_type) :: Time_step_atmos, Time_step_phys
real                    :: dt_atmos
integer :: sec
integer days, seconds

logical :: cold_start      = .false.       ! read in initial condition
logical :: do_dynamics     = .true.
logical :: dry_atmosphere     = .true.
integer :: p_split = 1
integer :: p_count

type(fv_atmos_type), allocatable :: Atm(:)
!-----------------------------------------------------------------------

contains

!#######################################################################

  subroutine atmosphere_init ( Time_init, Time, Time_step )

    type (time_type), intent(in) :: Time_step
    type (time_type), intent(in) :: Time_init
    type (time_type), intent(in) :: Time

    ! local:
    integer :: axes(4)
    integer :: ss, ds
    integer :: ntiles=1
    integer i,j, isc, iec, jsc, jec, isd, ied, jsd, jed
    real pp(2)

      integer :: f_unit, unit, ios ! Namelist reading
      character(80) :: filename

      filename = "input.nml"


      f_unit = get_unit()
      open (f_unit,file=filename)
      rewind (f_unit)
      read (f_unit,atmosphere_nml ,iostat=ios)
      unit = stdlog()
      write(unit, nml=atmosphere_nml)

      close(f_unit)


  !----- write version and namelist to log file -----

    call write_version_number ( version, tag )

  !---- compute physics/atmos time step in seconds ----

    Time_step_atmos = Time_step
    call get_time (Time_step_atmos, sec)
    dt_atmos = real(sec)

  !----- initialize FV dynamical core -----
!   cold_start = (.not.file_exist('INPUT/fv_rst.res.nc').and. .not.file_exist('INPUT/'//fms_tracers_file))
    cold_start = (.not.file_exist('INPUT/fv_core.res.nc'))

    allocate(Atm(ntiles))
    call fv_init(Atm(:),dt_atmos)  ! allocates Atm components

    Atm(1)%moist_phys = .false.

    ! Init model data
         call timing_on('fv_restart')
    call fv_restart(domain, Atm, dt_atmos, seconds, days, cold_start, grid_type)
         call timing_off('fv_restart')

! HI commented out line below, and added isd,ied,jsd,jed
!#ifdef PERTURB_IC
    isc = Atm(1)%isc
    iec = Atm(1)%iec
    jsc = Atm(1)%jsc
    jec = Atm(1)%jec
    isd = Atm(1)%isd
    ied = Atm(1)%ied
    jsd = Atm(1)%jsd
    jed = Atm(1)%jed
    
#ifdef SW_DYNAMICS
! TEST CASE-7
!     if ( .not. cold_start ) then
         do j=jsc,jec
            do i=isc,iec
               pp(1) = Atm(1)%agrid(i,j,1)
               pp(2) = Atm(1)%agrid(i,j,2)
               Atm(1)%delp(i,j,1) = Atm(1)%delp(i,j,1) + 120.*grav*cos(pp(2)) *  &
               exp( -(3.*(pp(1)-pi))**2 ) * exp( -(15.*(pp(2)-pi/4.))**2 )
            enddo
         enddo
!     endif
#endif
!#endif

! Need to implement time in fv_restart files
!   if ( .not. cold_start ) then 
! !
! ! Check consistency in Time
! !
!     fv_time = set_time (seconds, days)
!     call get_time (Time, ss,  ds)
!
!     if(seconds /= ss .or. days /= ds) then
!        unit = stdout()
!       write(unit,*) 'FMS:', ds, ss
!       write(unit,*) 'FV:', days, seconds
!       call error_mesg('FV_init:','Time inconsistent between fv_rst and INPUT/atmos_model.res',FATAL)
!     endif
!   else
      fv_time = time
!   endif



      call fv_diag_init(Atm, axes, Time, Atm(1)%npx, Atm(1)%npy, Atm(1)%npz, p_ref)

      call Exo_Tend_init(isc,iec,jsc,jec,Atm(1)%npz,Atm(1)%ncnst,Time,axes)

      ! Initialise physics
      call exo_init(isc,iec,jsc,jec,Atm(1)%npz, axes, Time, &
           Atm(1)%agrid(isc:iec,jsc:jec,1), Atm(1)%agrid(isc:iec,jsc:jec,2), &
           Time_step, Atm(1)%pt(isc:iec,jsc:jec,1:Atm(1)%npz),(mpp_pe()==mpp_root_pe()), &
           Atm(1)%ncnst, Atm(1)%delp(isc:iec,jsc:jec,1:Atm(1)%npz), &
           Atm(1)%peln(isc:iec,1:Atm(1)%npz+1, jsc:jec), &
           Atm(1)%q(isc:iec,jsc:jec,1:Atm(1)%npz, 1:Atm(1)%ncnst), cold_start,&
           Atm(1)%pk_nd(isc:iec,jsc:jec,1:Atm(1)%npz+1))

! MDH commented out
!    call lin_cld_microphys_init(axes, Time)

!   if( nlev > 1 ) call hs_forcing_init ( axes, Time )

!-----------------------------------------------------------------------

  end subroutine atmosphere_init


!#######################################################################

  subroutine atmosphere (Time)
    type(time_type), intent(in) :: Time

    real:: zvir
    real:: time_total
    real:: tau_winds, tau_press, tau_temp
    
#ifdef NUDGE_IC
    tau_winds =  1. * 3600.
    tau_press = -1.
    tau_temp  = -1.
#else
    tau_winds = -1.
    tau_press = -1.
    tau_temp  = -1.
#endif

    fv_time = Time + Time_step_atmos
    call get_time (fv_time, seconds,  days)

    time_total = days*SECONDS_PER_DAY + seconds

    if ( tau_winds>0. .or. tau_press>0. .or. tau_temp>0. )     &
    call  fv_nudge(Atm(1)%npz, Atm(1)%isc, Atm(1)%iec, Atm(1)%jsc, Atm(1)%jec, Atm(1)%ng, &
                   Atm(1)%u, Atm(1)%v, Atm(1)%delp, Atm(1)%pt, dt_atmos,    &
                   tau_winds, tau_press, tau_temp)

  !---- call fv dynamics -----
    !if ( adiabatic .or. Atm(1)%do_Held_Suarez ) then
    if ( dry_atmosphere) then
         zvir = 0.         ! no virtual effect
    else
         zvir = rvgas/rdgas - 1.
    endif

    call set_domain(Atm(1)%domain)  ! needed for diagnostic output done in fv_dynamics

    if (do_dynamics) then
    call timing_on('fv_dynamics')

    call fv_dynamics(Atm(1)%npx, Atm(1)%npy, Atm(1)%npz, Atm(1)%ncnst-Atm(1)%pnats,         &
                     Atm(1)%ng, dt_atmos, Atm(1)%consv_te,                                  &
                     Atm(1)%fill, Atm(1)%reproduce_sum, kappa, cp_air, zvir,                &
                     Atm(1)%ks, Atm(1)%ncnst, Atm(1)%n_split, Atm(1)%q_split,               &
                     Atm(1)%u, Atm(1)%v, Atm(1)%um, Atm(1)%vm,                              &
                     Atm(1)%w, Atm(1)%delz, Atm(1)%hydrostatic,                             &
                     Atm(1)%pt, Atm(1)%delp, Atm(1)%q, Atm(1)%ps,                           &
                     Atm(1)%pe, Atm(1)%pk, Atm(1)%pk_nd, Atm(1)%peln, Atm(1)%pkz, Atm(1)%phis,            &
                     Atm(1)%omga, Atm(1)%ua, Atm(1)%va, Atm(1)%uc, Atm(1)%vc,               &
                     Atm(1)%ak, Atm(1)%bk, Atm(1)%mfx, Atm(1)%mfy,                          &
                     Atm(1)%cx, Atm(1)%cy, Atm(1)%ze0, Atm(1)%hybrid_z, time_total)
                     
    call timing_off('fv_dynamics')
    end if

    if(Atm(1)%npz /=1 .and. .not. adiabatic)then
       call timing_on('FV_PHYS')

       do p_count = 1,p_split
!       print*, 'fvp'
       call fv_phys(Atm(1)%npx, Atm(1)%npy, Atm(1)%npz, Atm(1)%isc, Atm(1)%iec,          &
                    Atm(1)%jsc, Atm(1)%jec, Atm(1)%ng, Atm(1)%ncnst,                     &
                    Atm(1)%u, Atm(1)%v, Atm(1)%w, Atm(1)%pt, Atm(1)%q, Atm(1)%pe,        &
                    Atm(1)%delp, Atm(1)%peln, Atm(1)%pkz, dt_atmos/p_split,                      &
                    Atm(1)%ua, Atm(1)%va, Atm(1)%phis, Atm(1)%agrid,                     &
                    Atm(1)%ak, Atm(1)%bk, Atm(1)%ks, Atm(1)%ps, Atm(1)%pk, Atm(1)%pk_nd,               &
                    Atm(1)%u_srf, Atm(1)%v_srf, Atm(1)%tsurf, Atm(1)%delz,                  &
                    Atm(1)%hydrostatic, Atm(1)%phys_hydrostatic, Atm(1)%oro, .true.,     &
                    .false., p_ref, Atm(1)%fv_sg_adj, (mpp_pe()==mpp_root_pe()),         &
                    Atm(1)%do_Held_Suarez, fv_time, time_total-dt_atmos+p_count*dt_atmos/p_split, dt_atmos)
       end do
    call timing_off('FV_PHYS')
    endif

    call nullify_domain()

  !---- diagnostics for FV dynamics -----

    call timing_on('FV_DIAG')

    call fv_diag(Atm, zvir, fv_time, Atm(1)%print_freq)

    call timing_off('FV_DIAG')

 end subroutine atmosphere


 subroutine atmosphere_end

    call get_time (fv_time, seconds,  days)

! MDH commented out
!    if ( Atm(1)%nwat==6 )    & 
!    call lin_cld_microphys_end

    call fv_end(Atm)
    deallocate(Atm)

  end subroutine atmosphere_end

 subroutine atmosphere_domain ( fv_domain )
 type(domain2d), intent(out) :: fv_domain

!  returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos
        
   fv_domain = domain
        
 end subroutine atmosphere_domain

end module atmosphere_mod
