module exo_phys_init_mod

  !-----------------------------------------------------------------------
  ! Paths:
  ! constants_mod:      src/shared/constants/constants.F90
  ! mpp_update_domains: src/shared/mpp/mpp_domains.F90
  ! mpp_mod:            src/shared/mpp/mpp.F90
  ! time_manager_mod:   src/shared/time_manager/time_manager.F90
  ! diag_manager_mod:   src/shared/diag_manager/diag_manager.F90
  ! fv_grid_tools_mod:  src/atmos_cubed_sphere/tools/fv_grid_tools.F90
  ! fv_grid_utils_mod:  src/atmos_cubed_sphere/model/fv_grid_utils.F90
  ! fv_mp_mod:          src/atmos_cubed_sphere/tools/fv_mp_mod.F90
  ! fv_diagnostics_mod: src/atmos_cubed_sphere/tools/fv_diagnostics.F90
  ! fv_timing_mod:      src/atmos_cubed_sphere/tools/fv_timing.F90

  !-----------------------------------------------------------------------
  use constants_mod,      only: grav, rdgas, rvgas, RADIAN, kappa, radius, stefan
  use mpp_domains_mod,    only: mpp_update_domains
  use mpp_mod,            only: FATAL, mpp_error, mpp_pe, stdlog, &
                                 mpp_npes, mpp_get_current_pelist, get_unit
  use time_manager_mod,   only: time_type, get_date, get_time
  use diag_manager_mod,   only: send_data
  
  use fv_grid_utils_mod,  only: g_sum
  use fv_diagnostics_mod, only: prt_maxmin
  use fv_timing_mod,      only: timing_on, timing_off
  use rad_coupler_mod,    only: rad_coupler_init
  use moisture_init_mod,  only: init_moisture
  use netcdf_reg_mod,     only: netcdf_reg_RT
  use dry_convection_mod,    only: dry_convection_init
  use exo_phys_mod, only: do_ding_convection
  use fv_mp_mod, only : is_master
  
  implicit none

  private
  public :: Exo_Init


contains

  subroutine Exo_Init(is, ie, js, je,&
                      num_levels, axes, Time, lon, lat, delta_t, T, &
                      nq, delp, peln, q, cold_start, non_dilute)

    integer, intent(in), dimension(4) :: axes
    type(time_type), intent(in)       :: Time
    integer, intent(in)               :: is,ie,js,je,num_levels,nq
    real, intent(in), dimension(:,:)  :: lat, lon
    type(time_type), intent(in)       :: delta_t
    real, intent(inout), dimension(is:ie,js:je,num_levels) :: T,q
    real, intent(in), dimension(is:ie,js:je,num_levels) :: delp
    real, intent(in), dimension(is:ie,num_levels+1, js:je) :: peln
    logical, intent(in) :: cold_start, non_dilute
    
    
    ! Local arguments
    integer :: k
    character(80) :: filename

if (is_master()) then
print*, ' '
print*, ' _____         _________  ___ _____ '
print*, '|  ___|        |  ___|  \/  |/  ___|'
print*, '| |____  _____ | |_  | .  . |\ `--. '
print*, '|  __\ \/ / _ \|  _| | |\/| | `--. \'
print*, '| |___>  < (_) | |   | |  | |/\__/ /'
print*, '\____/_/\_\___/\_|   \_|  |_/\____/ '
print*, ' '
print*, ' '
end if
if (is_master()) write(*,*) 'before initialise condensation'
    ! Initialise condensation scheme

    if (do_ding_convection) then
 if (is_master()) write(*,*) 'before init_moisture'        
       call init_moisture(num_levels,is,ie,js,je,nq, q, peln, delp, T, axes, &
            Time, cold_start, non_dilute)
 if (is_master()) write(*,*) 'after init_moisture'
    endif

    ! Initialise the radiation routine
    !call radiation_init(is, ie, js, je, num_levels, axes, Time)

    ! Initialise the dry adjustment
    !call dry_adj_init()
if (is_master())write(*,*) 'before convection init'
    call dry_convection_init(axes, Time)   
if (is_master()) write(*,*) 'before rad coupler init'
!    call netcdf_reg_RT(axes,Time,scheme)
    call  rad_coupler_init(axes,Time)
    
#ifdef BUILD_SOC
    ! Initialise Socrates
    !call socrates_init(is, ie, js, je, num_levels, axes, Time, lat, lon, lat, delta_t)
#endif
    if (is_master()) write(*,*) 'end of exo_init'
  end subroutine Exo_Init

end module exo_phys_init_mod

