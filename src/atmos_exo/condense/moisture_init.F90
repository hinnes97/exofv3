module moisture_init_mod
  use sat_props_mod, only: q_sat
  use constants_mod,  only: RVGAS, WTMAIR, WTMH2O, CP_AIR, HLV, HLS, CP_VAPOR, RDGAS
  use tracer_manager_mod, only: get_tracer_index
  use field_manager_mod,  only: MODEL_ATMOS
  use fms_mod,        only:  mpp_pe,  check_nml_error
  use time_manager_mod,   only: time_type
  use diag_manager_mod,   only: register_diag_field, send_data
  use fv_mp_mod,        only: is_master
  use mpp_mod,          only: FATAL, mpp_error, mpp_pe, stdlog, &
                                  mpp_npes, mpp_get_current_pelist, get_unit
  use interp_mod, only : get_pk_edge
  
  implicit none
    real, dimension(:,:,:), allocatable :: lheat
  character(len=128) :: version='$Id: moisture_init.f90 $'
  character(len=128) :: tag = 'homemade'

  ! Namelist constants
  real :: tau_c ! Condensation timescale
  real :: tau_e ! Evaporation timescale
  real :: q0    ! Max q
  real :: tau_r ! Vapour replenishment timescale
  real :: p_r   ! Pressure below which replenishment occurs
  
  ! Other parameters
  real :: T0 = 273.1575, p0 = 611.657 ! Triple point of water
  real, parameter :: eps = WTMH2O/WTMAIR
  real :: varpi = 1 - 1/eps

  ! Diagnostics
  integer :: id_lheat, id_qsat

  character(len=8) :: mod_name= 'condense'

  real :: missing_value = -1.e10
  
  private 
  public  :: init_moisture, lheat, cond_end

  namelist/init_moisture_nml/q0, tau_r, p_r
  
contains

  subroutine init_moisture(npz, is, ie, js, je, nq, q, peln, delp, pt, &
       cold_start , non_dilute)
    !----------------------------------------------------------------------------!
    !                                                                            !
    !          Initialise moisture to saturation                                 !
    !                                                                            !
    !----------------------------------------------------------------------------!


    ! Imported variables ---------------------------------------------------------

    integer, intent(in   ) :: npz                    ! No. vertical levels
    integer, intent(in   ) :: is, ie, js, je, nq      ! For array indices
    logical, intent(in   ) :: cold_start, non_dilute
    

    real,    intent(inout) ::     q(is:ie,js:je,1:npz,1:nq) ! tracers
    real,    intent(in   ) ::    pt(is:ie,js:je,1:npz)    ! temperature
    real,    intent(in   ) ::  delp(is:ie,js:je,1:npz)    ! Pressure thickness
    real,    intent(in   ) ::    peln(is:ie,1:npz+1,js:je)! Log interface p levels


    ! Local variables ------------------------------------------------------------
    real :: qsat(is:ie,js:je,1:npz) ! Saturation concentration
    real :: pfull(is:ie,js:je,1:npz) ! Pressure interface levels
    real :: ratio(is:ie,js:je,1:npz)
    real :: qc(is:ie,js:je,1:npz)
    real :: kap_loc(is:ie,js:je,1:npz)
    integer :: i,j,k
    integer :: vap, cond, sphum
    integer :: ierr, ios, f_unit, unit,ind(3)
    logical :: master


    ! Main body of function -----------------------------------------------------
    master = is_master()

    f_unit = get_unit()
    open (f_unit,file='input.nml')
       ! Read initialisation namelist
    rewind (f_unit)
    read (f_unit,init_moisture_nml ,iostat=ios)
    if (ios .gt. 0) then
       if(master) write(6,*) 'initialisation_nml ERROR: reading ',trim('input.nml'),', iostat=',ios
       call mpp_error(FATAL,'FV core terminating')
    endif
    unit = stdlog()
    write(unit, nml=init_moisture_nml)
    close(f_unit)

    open (unit, file='logfile.out',  position='append')
    if ( mpp_pe() ==0 ) then
       write (unit, '(/,80("="),/(a))') trim(version), trim(tag)
       write (unit, nml=init_moisture_nml)
    endif
    close(unit)

!    if (is_master()) write(*,*) '2', minval(pt(is:ie,js:je,1:npz))
    allocate( lheat(is:ie,js:je,1:npz) )

    if(cold_start) then       
       vap = get_tracer_index(MODEL_ATMOS, 'vapour')
       cond   = get_tracer_index(MODEL_ATMOS, 'condensate' )
!
       do k=1,npz
          do j=js,je
             do i=is,ie
                pfull(i,j,k) = delp(i,j,k)/(peln(i,k+1,j) - peln(i,k,j))
             end do
          end do
       end do

!       if (is_master()) write(*,*) '3', minval(pt(is:ie,js:je,1:npz))
!       ! Get specific humidity
       qc(is:ie,js:je,1:npz) = 0.0
!
!       if (is_master()) write(*,*) '4.1', minval(pt(is:ie,js:je,1:npz))
       do j=js,je
          do i=is,ie
             call q_sat(pfull(i,j,:), pt(i,j,:), qc(i,j,:), qsat(i,j,:))
          end do
       end do
!
!       if (is_master()) write(*,*) '4.2', minval(pt(is:ie,js:je,1:npz)), vap, cond, q0, shape(q(is:ie,js:je,:,:))
!       if (is_master()) write(*,*) is,ie,js,je,npz,nq
! Set q equal to a uniform, well mixed dilute limit

!       
       do k=1,npz
          do j=js,je
             do i=is,ie
                q(i,j,k,vap) = min(qsat(i,j,k), q0)
                q(i,j,k,cond) = 0.0
             enddo
          enddo
       enddo

       
       !q(:,:,:,cond) = 0.0
       !q(is:ie,js:je,:,cond) = 0.0
!       if (master) then
!          write(*,*) 'INIT, MAXVAL QVAP', maxval(q(is:ie,js:je,:,vap)),&
!            'INIT, MAXVAL QCOND', maxval(q(is:ie,js:je,:,cond))
!          write(*,*) 'maxval pt, satq, p'
!          !write(*,*) maxval(qsat(is:ie,js:je,1:npz)), maxval(pfull(is:ie,js:je,1:npz)), maxval(pt(is:ie,js:je,1:npz))
!          ind = maxloc(q(is:ie,js:je,1:npz,vap))
!          write(*,*) ind
!          !do k=1,npz
!          !   write(*,*) q(ind(1)+is-1, ind(2)+js-1, k,vap), pt(ind(1)+is-1,ind(2)+js-1,k)
!          !enddo
!          
!       endif
       
!       if (is_master()) write(*,*) '4.3', minval(pt(is:ie,js:je,1:npz))
       ! Cold trap the q value
       do k=npz-1,1,-1
          do j=js,je
             do i = is,ie
                q(i,j,k,vap) = min(q(i,j,k,vap), q(i,j,k+1,vap))
             enddo
          enddo
       enddo
    endif

!    if (is_master()) write(*,*) '5', minval(pt(is:ie,js:je,1:npz))
    if (non_dilute) then
       do k=1,npz
          do j=js,je
             do i=is,ie
                kap_loc(i,j,k) = (q(i,j,k,vap)*rvgas + (1.-q(i,j,k,vap))*rdgas)/&
                     (q(i,j,k,vap)*cp_air + (1.-q(i,j,k,vap))*cp_vapor)
             enddo
          enddo
       enddo
    endif
!
!    if (is_master()) write(*,*) '6', minval(pt(is:ie,js:je,1:npz))
!    if (master) then
!       write(*,*) 'INIT, MAXVAL QVAP 2', maxval(q(is:ie,js:je,:,vap)),&
!            'INIT, MAXVAL QCOND 2', maxval(q(is:ie,js:je,:,cond))
!       if (is_master()) write(*,*) '7', minval(pt(is:ie,js:je,1:npz))
!       endif
!    ! Initialise diagnostics
!    id_lheat = register_diag_field(mod_name, 'lheat', axes(1:3), Time, &
!         'Latent heating rate', 'K/s', missing_value=missing_value)
!    id_qsat = register_diag_field(mod_name, 'qsat', axes(1:3), Time, &
!         'Saturation specific humidity', 'kg/kg', missing_value=missing_value)
!
!if (is_master()) write(*,*) '8', minval(pt(is:ie,js:je,1:npz))
  end subroutine init_moisture

  subroutine cond_end
    deallocate(lheat)
  end subroutine cond_end
end module moisture_init_mod
