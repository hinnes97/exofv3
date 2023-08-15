module condense_mod

  use constants_mod,  only: RVGAS, WTMAIR, WTMH2O, CP_AIR, HLV, HLS
  use tracer_manager_mod, only: get_tracer_index
  use field_manager_mod,  only: MODEL_ATMOS
  use fms_mod,        only:  mpp_pe, check_nml_error
  use fms2_io_mod,    only: open_file, close_file
  use time_manager_mod,   only: time_type
  use diag_manager_mod,   only: register_diag_field, send_data
  use fv_mp_mod,        only: is_master
  use mpp_mod,          only: FATAL, mpp_error, mpp_pe, stdlog, &
                                  mpp_npes, mpp_get_current_pelist, get_unit


  implicit none

  real, dimension(:,:,:), allocatable :: lheat
  character(len=128) :: version='$Id: condense.f90 $'
  character(len=128) :: tag = 'homemade'

  ! Namelist constants
  real :: tau_c ! Condensation timescale
  real :: tau_e ! Evaporation timescale
  real :: q0    ! Max q
  real :: tau_r ! Vapour replenishment timescale
  real :: p_r   ! Pressure below which replenishment occurs
  logical :: nondilute ! Whether to use non-dilute physics

  ! Other parameters
  real :: T0 = 273.1575, p0 = 611.657 ! Triple point of water
  real :: eps = WTMH2O/WTMAIR

  ! Diagnostics
  integer :: id_lheat, id_qsat

  character(len=8) :: mod_name= 'condense'

  real :: missing_value = -1.e10
  
  private 
  public  :: cond_tend, cond_init, lheat, cond_end

  namelist/condense_nml/tau_c, tau_e, q0, tau_r, p_r, nondilute
  
contains

  !-------------------------------------------------------------------------------

  subroutine cond_init(npz, is, ie, js, je, nq, q, peln, delp, pt, axes, &
                       Time, cold_start, non_dilute)

    !----------------------------------------------------------------------------!
    !                                                                            !
    !          Initialise moisture to saturation                                 !
    !                                                                          v  !
    !----------------------------------------------------------------------------!


    ! Imported variables ---------------------------------------------------------

    integer, intent(in   ) :: npz                    ! No. vertical levels
    integer, intent(in   ) :: is, ie, js, je, nq ! For array indices
    integer, intent(in   ) :: axes(4)
    logical, intent(in   ) :: cold_start

    real,    intent(inout) ::     q(is:ie,js:je,npz,nq) ! tracers
    real,    intent(in   ) ::    pt(is:ie,js:je,npz)    ! temperature
    real,    intent(in   ) ::  delp(is:ie,js:je,npz)    ! Pressure thickness
    real,    intent(in   ) ::    peln(is:ie,1:npz+1,js:je)! Log interface p levels

    type(time_type), intent(in) :: Time

    logical, intent(inout) :: non_dilute

    ! Local variables ------------------------------------------------------------
    real :: q_sat(is:ie,js:je,npz) ! Saturation concentration
    real :: pfull(is:ie,js:je,npz) ! Pressure interface levels
    real :: ratio(is:ie,js:je,npz)

    integer :: i,j,k
    integer :: vap, cond, sphum
    integer :: ierr, ios, f_unit, unit
    logical :: master

    ! Main body of function -----------------------------------------------------
    master = is_master()
    f_unit = get_unit()
    open (f_unit,file='input.nml')
       ! Read initialisation namelist
    rewind (f_unit)
    read (f_unit,condense_nml ,iostat=ios)
    if (ios .gt. 0) then
       if(master) write(6,*) 'initialisation_nml ERROR: reading ',trim('input.nml'),', iostat=',ios
       call mpp_error(FATAL,'FV core terminating')
    endif
    unit = stdlog()
    write(unit, nml=condense_nml)
    close(f_unit)

    open (unit, file='logfile.out', position = 'append')
    if ( mpp_pe() ==0 ) then
       write (unit, '(/,80("="),/(a))') trim(version), trim(tag)
       write (unit, nml=condense_nml)
    endif
    close(unit)

    allocate( lheat(is:ie,js:je,1:npz) )

    if(cold_start) then       
       vap = get_tracer_index(MODEL_ATMOS, 'vapour')
       cond   = get_tracer_index(MODEL_ATMOS, 'condensate' )

       do k=1,npz
          do j=js,je
             do i=is,ie
                pfull(i,j,k) = delp(i,j,k)/(peln(i,k+1,j) - peln(i,k,j))
             end do
          end do
       end do

       ! Get specific humidity
       call get_qsat(is,ie,js,je,npz,pt, pfull, ratio, q_sat)

      ! Set q equal to a uniform, well mixed dilute limit
       q(is:ie,js:je,:,vap) = min(q_sat, q0)

       ! Cold trap the q value
       do k=npz-1,1,-1
          do j=js,je
             do i = is,ie
                q(i,j,k,vap) = min(q(i,j,k,vap), q(i,j,k+1,vap))
             enddo
          enddo
       enddo
    endif

    ! For testing, set q to same value as sphum everywhere:
    !cond = get_tracer_index(MODEL_ATMOS, 'sphum')
    !q(:,:,:,vap) = q(:,:,:,cond)
    ! Setup non_dilute to be passed back to atmosphere mod
    non_dilute = nondilute

    ! Initialise diagnostics
    id_lheat = register_diag_field(mod_name, 'lheat', axes(1:3), Time, &
         'Latent heating rate', 'K/s', missing_value=missing_value)
    id_qsat = register_diag_field(mod_name, 'qsat', axes(1:3), Time, &
         'Saturation specific humidity', 'kg/kg', missing_value=missing_value)

    
  end subroutine cond_init

  subroutine cond_tend(npz, is, ie, js, je, nq, q, pfull, q_dt, pt, t_dt, Time,&
       pdt, delp)

    !----------------------------------------------------------------------------!
    !                                                                            !
    !          Compute tracer tendency to relax to saturation                    !
    !                                                                            !
    !----------------------------------------------------------------------------!

    ! Imported variables ---------------------------------------------------------

    integer, intent(in   ) :: npz                    ! No. vertical levels
    integer, intent(in   ) :: is, ie, js, je, nq ! For array indices 

    real,    intent(in   ) ::   pdt                     ! timestep
    real,    intent(in   ) ::     q(is:ie,js:je,npz,nq) ! tracers
    real,    intent(in   ) ::    pt(is:ie,js:je,npz)    ! temperature
    real,    intent(inout) :: pfull(is:ie,js:je,npz)    ! pressure
    real,    intent(inout) ::  t_dt(is:ie,js:je,npz)    ! temp. tend.
    real,    intent(inout) ::  q_dt(is:ie,js:je,npz,nq) ! trace. tend.
    real,    intent(inout) ::  delp(is:ie,js:je,npz)    ! pressure thickness

    type(time_type), intent(in) :: Time
    
    ! Local variables ------------------------------------------------------------

    real :: L   ! Latent heat of vaporisation
    real :: grad    ! Gradient in saturation vapour concentration wrt T
    real :: evap

    real :: q_sat(is:ie,js:je,npz) ! Saturation specific humidity
    real :: ratio(is:ie,js:je,npz) ! Ratio between psat and pfull
    real :: delq(is:ie,js:je,npz)  ! Change in q


    integer :: i,j,k
    integer :: vap, cond ! Vapour and condensate indices

    logical :: used
    
    ! Main body of function ------------------------------------------------------
    ! Get index of vapour and condensate tracers

    vap = get_tracer_index(MODEL_ATMOS, 'vapour')
    cond   = get_tracer_index(MODEL_ATMOS, 'condensate'  )
        
    call get_qsat(is,ie,js,je,npz,pt,pfull,ratio,q_sat)

    ! Relax towards saturation if exceeding saturation concentration
    ! Use Frierson 2006 relaxation so that latent heat release doesn't
    ! cause undersaturation

    lheat = 0. ! Just in case no if statements are triggered
    
    do k=1,npz
       do j=js,je
          do i=is,ie
             if (pt(i,j,k) > T0) then
                L = HLV
             else
                L = HLS
             endif

             if (q(i,j,k,vap) > q_sat(i,j,k) ) then
                ! Condensation if supersaturated
                grad = (L/(RVGAS*pt(i,j,k)**2))*eps*ratio(i,j,k) &
                     /(1+(eps-1)*ratio(i,j,k))**2
                
                delq(i,j,k) = (q_sat(i,j,k) - q(i,j,k,vap))/(1 + L/CP_AIR*grad) * pdt/tau_c
                q_dt(i,j,k,vap) = q_dt(i,j,k,vap) + delq(i,j,k)/pdt
                lheat(i,j,k) = -delq(i,j,k)*L/pdt/CP_AIR
                t_dt(i,j,k) = t_dt(i,j,k) + lheat(i,j,k)

             else
                delq(i,j,k) = (q_sat(i,j,k) - q(i,j,k,vap))/(1 + L/CP_AIR*grad) * pdt/tau_e
             endif

             if (pfull(i,j,k) > p_r) then
                q_dt(i,j,k,vap) = q_dt(i,j,k,vap) + (q0 - q(i,j,k,vap))/tau_r
             endif
             
          end do
       end do
    end do

    call rain(is,ie,js,je,npz,nq,q_dt, delq, delp, lheat, pt, pdt)

    ! Do diagnostics
    if ( id_lheat > 0 ) then
       used = send_data ( id_lheat, lheat, Time)
    endif
    
    if ( id_qsat > 0 ) then
       used = send_data (id_qsat, q_sat, Time)
    endif

    
  end subroutine cond_tend

  subroutine get_qsat(is,ie,js,je,npz,pt, pfull, x_sat,q_sat)
    integer, intent(in) :: is,ie,js,je,npz
    real,    intent(in) :: pt(is:ie,js:je,npz)
    real,    intent(in) :: pfull(is:ie,js:je,npz)

    real,    intent(out) :: q_sat(is:ie,js:je,npz) ! Saturation mass concentration
    real,    intent(out) :: x_sat(is:ie,js:je,npz) ! Saturation molar concentration

    ! Local variables --------------------------------------------------

    real :: p_sat(is:ie,js:je,npz)

    where( pt(is:ie,js:je,:) > T0 )
       p_sat =  p0*exp(-HLV/RVGAS*(1/pt(is:ie,js:je,:) - 1/T0))
    elsewhere
       p_sat =  p0*exp(-HLS/RVGAS*(1/pt(is:ie,js:je,:) - 1/T0))
    endwhere

    x_sat = p_sat/pfull

    q_sat = eps*x_sat/(1 + (eps - 1)*x_sat)
    
  end subroutine get_qsat

  subroutine rain(is,ie,js,je,npz,nq,q_dt, delq, delp, lheat, pt, pdt)
    integer, intent(in) :: is,ie,js,je,npz,nq
    real, intent(in)    :: pdt
    real, intent(in)    :: delp(is:ie,js:je,npz)
    real, intent(in)    :: delq(is:ie,js:je,npz)
    real, intent(inout) :: q_dt(is:ie,js:je,npz,nq)
    real, intent(inout) :: lheat(is:ie,js:je,npz)
    real, intent(in   ) :: pt(is:ie,js:je,npz)
    
    integer:: k, vap, cond
    real :: avail_mass(is:ie,js:je)
    real :: q_change(is:ie,js:je,npz)
    real :: L(is:ie,js:je,npz)
    
    vap = get_tracer_index(MODEL_ATMOS, 'vapour')

    where (pt > T0)
       L = HLV
    elsewhere
       L = HLS
    endwhere

    avail_mass = 0.
    q_change = 0.
    
    do k=1,npz
      
       ! Add lost vapour to available evaporation supply
       where (delq(:,:,k) < 0.0)
          avail_mass(:,:) = avail_mass(:,:) &
               - delq(:,:,k)*delp(:,:,k)
       endwhere

       where ((delq(:,:,k) > 0.0) .and. (avail_mass(:,:) > 0.0))
          q_change(:,:,k) = min( avail_mass(:,:)/delp(:,:,k), delq(:,:,k) ) 
          q_dt(:,:,k,vap) = q_dt(:,:,k,vap) + q_change(:,:,k)/pdt
          avail_mass(:,:) = avail_mass(:,:) - q_change(:,:,k)*delp(:,:,k)
          lheat(:,:,k) = lheat(:,:,k) - q_change(:,:,k) * L(:,:,k)/CP_AIR/pdt
       endwhere
       
    end do
  end subroutine rain
  
  subroutine cond_end
    deallocate(lheat)
  end subroutine cond_end
  
end module condense_mod
