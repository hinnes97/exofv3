module exo_phys_init_mod

  use fv_mp_mod,          only: is_master
  use mpp_domains_mod,    only: domain2d
  use fv_grid_utils_mod,      only: g_sum
  use fv_diagnostics_mod,  only : prt_mass
  use fv_arrays_mod,      only: R_GRID
  use exo_phys_mod,       only: do_ding_convection, do_dry_convection
  use moisture_init_mod,  only: init_moisture
  use rad_coupler_mod,    only: rad_coupler_init
  use dry_convection_mod, only: dry_convection_init
  use time_manager_mod,   only: time_type
  
contains

  subroutine exo_init(is,ie,js,je,npz, nq, pt, delp, q, peln, &
        cold_start, non_dilute,&
       axes, Time)!, domain, area, ps)

    integer, intent(in) :: is,ie,js,je,npz,nq, axes(4)

    real, intent(in), dimension(is:ie,js:je,1:npz) :: pt
    real, intent(in), dimension(is:ie,js:je,1:npz) :: delp
    real, intent(inout), dimension(is:ie,js:je,1:npz,1:nq) :: q
    real, intent(in), dimension(is:ie,1:npz+1,js:je) :: peln
        !real(kind=R_GRID), intent(IN):: area(is-ng:ie+ng,js-ng:je+ng)
    logical, intent(in) :: cold_start, non_dilute

    type(time_type), intent(in) :: Time

    !type(domain2d), intent(in) :: domain
    integer :: k
    real :: avps
    if (is_master()) then
       write(*,*) ' '
       write(*,*) ' ____  _  _  __  ____  _  _  ____  '
       write(*,*) '(  __)( \/ )/  \(  __)/ )( \( __ \ '
       write(*,*) ' ) _)  )  ((  O )) _) \ \/ / (__ ( '
       write(*,*) '(____)(_/\_)\__/(__)   \__/ (____/ '
       write(*,*) ' '
       write(*,*) ' '
    endif
    
    if (do_ding_convection) then
       call init_moisture(npz, is,ie,js,je,nq,q, peln, delp, pt, cold_start, non_dilute)
    else if (do_dry_convection) then
       call dry_convection_init(axes, Time)
    endif

    call rad_coupler_init(axes, Time)
  end subroutine exo_init
  
end module exo_phys_init_mod
