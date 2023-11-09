module ding_convection
  use phys, only : T_TP => H2O_TriplePointT, P_TP => H2O_TriplePointP, L_sub => H2O_L_sublimation, &
       L_vap => H2O_L_vaporization_TriplePoint, CP_v => H2O_cp, &
       mu_d => H2He_solar_MolecularWeight, mu_v => H2O_MolecularWeight, Rstar
  use sat_props_mod, only: q_sat, q_sat_single, r_sat_single, r_sat, sat_vp, cp_liquid, hv_tp, hl_tp,&
       get_f_ice, mv_mvmd
  
!  use tables, only : find_var_lin, lheat, hv, cpv_new, cpl_new, Ttp, find_var_simplified  use constants, only : grav, cp_vapor
  use tracer_manager_mod, only: get_tracer_index
  use constants_mod, only: grav, rvgas, rdgas, cp_d => cp_air
  use field_manager_mod,  only: MODEL_ATMOS
  use time_manager_mod, only: time_type
  use fv_mp_mod, only: is_master
  use netcdf_reg_mod, only: netcdf_reg_ding
  use omp_lib
  
  implicit none
  

  real :: precision = 1.e-12
  real :: max_dt = 5.0
  integer :: n_iter_conv = 500
  integer :: n_iter_newt = 100
  real :: tolerance = 1.e-6
  integer :: iv, ic  ! Tracer indices for vapour and condensates

  logical:: master
  
contains

  subroutine ding_adjust_init(Time, axes)
    type(time_type), intent(in)       :: Time
    integer, dimension(4), intent(in) :: axes
    
    iv = get_tracer_index(MODEL_ATMOS, 'vapour')
!ic = get_tracer_index(MODEL_ATMOS, 'condensate')

    call netcdf_reg_ding(axes, Time)
    
  end subroutine ding_adjust_init
  
  subroutine ding_adjust(p, delp, T, q, t_dt_ding, q_dt_ding, q_liq, q_ice)
    !==========================================================================
    ! Description
    !==========================================================================
    ! Performs moist adiabatic adjustment pairwise, conserving non-dilute moist
    ! enthalpy using Newton iteration, as in Ding & Pierrehumbert 2016

    !==========================================================================
    ! Input variables
    !==========================================================================    
    real, intent(in)   , dimension(:) :: p! Pressure
    real, intent(in), dimension(:) :: delp ! Pressure thickness will change only in
                                           ! condensation steps
    real, intent(in), dimension(:) :: T ! Temperature and specific humid.
    real, intent(in), dimension(:,:) :: q ! Tracers
        
    !==========================================================================
    ! Output variables 
    !==========================================================================

    real, intent(inout) :: t_dt_ding(:)   ! Temperature tendency
    real, intent(inout) :: q_dt_ding(:) ! Condensate and 
    real, intent(inout) :: q_liq(:), q_ice(:)

    !==========================================================================
    ! Local variables
    !==========================================================================
    real :: T_tmp(size(p)), qv_tmp(size(p)), qc_tmp(size(p)), qi_tmp(size(p))
    
    real :: dlnTdlnp, temp, T_lift,rho_lift,q_lift, pfact, grads(size(p)), q_lifts(size(p))
    real :: qsats(size(p)), rhos(size(p)),  qt(size(p)),  q_cond(size(p)),&
            T_virt(size(p)), T_virt_lift, R_lift,T_virt_lifts(size(p)), T_lifts(size(p)), R_lifts(size(p))
    real :: R_local(size(p)), cp_local(size(p)), qsatlift, R_loaded(size(p))
    integer :: k, npz, m,n
    logical ::conv, check

    !write(*,*) 'BEGINNING OF DING ADJUST'
    !write(*,*) maxval(q(:,iv)), maxval(q(:,ic)), iv, ic
    
    !==========================================================================
    ! Main body
    !==========================================================================


    npz = size(p)
    
    qt = q(:,iv) !+q(:,ic)
    
    ! ENTHALPY TEST
!    h_before=0.0
!    m_bef=0.0
!    
!    do k=1,npz
!       call calc_enthalpy(T(k), q(k,iv), 0.0, 0.0, h_temp)
!       hvert0(k) = h_temp
!       h_before = h_before + h_temp*delp(k)
!       m_bef = m_bef + q(k,iv)*delp(k)
!    enddo

    ! Only do adjustment on copy of profiles
    T_tmp = T
    qv_tmp = q(:,iv)
    qc_tmp = q_liq!qc_tmp = q(:,ic)
    qi_tmp = q_ice
    q_cond = q_liq + q_ice

    conv = .false.
    n= 1

    check = .true.
    do while (n < N_iter_conv+1 .and. check)
       check  = .false.
       ! Find saturation vapour specific humidity
       call q_sat(p,T_tmp,q_cond,qsats)

       cp_local = cp_d*(1-qv_tmp-qc_tmp-qi_tmp) + cp_v*qv_tmp + cp_liquid*qc_tmp + cp_v*qi_tmp
       R_local = Rstar*(qv_tmp/mu_v/(1-qc_tmp-qi_tmp) + (1-qv_tmp - qc_tmp - qi_tmp)/(1-qc_tmp-qi_tmp)/mu_d)
       R_loaded = R_local*(1 - qc_tmp - qi_tmp)
       rhos = p/R_local/T
       T_virt = R_local*T_tmp/rdgas

       ! Upwards sweep
       do k=npz-1,1,-1

          ! Do dry lift to check for condensation
          dlnTdlnp = R_local(k+1)/cp_local(k+1)
          pfact = exp(dlnTdlnp*log(p(k)/p(k+1)))
          T_lift = T_tmp(k+1)*pfact

          call q_sat_single(p(k), T_lift, qi_tmp(k+1) + qc_tmp(k+1), qsatlift)
! If saturated, gradient should be moist gradient
!          if (is_master() .and. omp_get_thread_num().eq. 0)  then
!             write(*,*) k, p(k) ,qsats(k), qsatlift
!             write(*,*) 'qv_tmp(k)', qv_tmp(k)
!             write(*,*) 'qsats(k)', qsats(k)
!             write(*,*) 'qsatlift', qsatlift
!             write(*,*) k, qv_tmp(k)/qsats(k)
!             write(*,*) T(k)
!             call sat_vp(T(k), temp)
!             write(*,*) temp, p(k)
!             conv = .true.
!          endif
          if ((qv_tmp(k+1) .gt. qsats(k+1)  .or. qv_tmp(k) .gt. qsats(k)) .and.&
               (qsatlift < 1.e19 .and. qsats(k+1) < 1.e19 .and. qsats(k) < 1.e19)) then
                       if (qsats(k+1) .lt. 0) then
                write(*,*) 'QSATS<0', qsats(k+1)
             endif
             
             ! get moist gradient
             call gradient(p(k+1), T_tmp(k+1), dlnTdlnp,temp)
             grads(k+1) = dlnTdlnp
             ! Calculate lifted temperature from adiabatic expansion
             pfact = exp(dlnTdlnp*log(p(k)/p(k+1)))
             T_lift = T_tmp(k+1)*pfact
             T_lifts(k) = T_lift
             ! Calculate q here as proportion of gas phase rv/(1+rv)
             call r_sat_single(p(k), T_lift, q_lift)
             q_lift = q_lift/(1 + q_lift)
             q_lifts(k) = q_lift
             !call q_sat_single(p(k), T_lift, 0.0, q_lift)
             R_lift = Rstar*(q_lift/mu_v + (1. - q_lift)/mu_d)
             R_lifts(k) = R_lift
             T_virt_lift = T_lift * R_lift/rdgas
             T_virt_lifts(k) = T_virt_lift
!             if (is_master() .and. omp_get_thread_num() .eq. 0) then
!                write(*,*) 'IN CONVECTION'
!                   write(*,*) k, T_virt(k), T_virt_lift
!                   write(*,*) k, T_tmp(k), T_lift
!                   conv = .true.
!             endif
             
             if (T_virt_lift .gt. T_virt(k)*(1 + tolerance)) then
                check = .true.
! Routine does pairwise moist adjustment
                !if (is_master() .and. omp_get_thread_num().eq. 0)  write(*,*) 'moist conv, k', k
                call conserve_moist_enthalpy(p(k), p(k+1), T_tmp(k), T_tmp(k+1), qv_tmp(k), qv_tmp(k+1), &
                     delp(k), delp(k+1), qc_tmp(k), qc_tmp(k+1), qi_tmp(k), qi_tmp(k+1))
                
! Recalculate the new rho, qsat etc
!                if (is_master() .and. omp_get_thread_num().eq. 0)  then
!                   write(*,*) 'CONVECTION', n
!                   call gradient(p(k+1), T_tmp(k+1), dlnTdlnp)
!                   write(*,*) k, dlnTdlnp, (log(T_tmp(k+1)) - log(T_Tmp(k)))/(log(p(k+1)) - log(p(k)))
!                   conv = .true.
!                endif
                do m=k,k+1
                   R_local(m) = Rstar*(qv_tmp(m)/mu_v/(1 - qc_tmp(m) - qi_tmp(m)) + &
                        (1 - qv_tmp(m) - qi_tmp(m) - qc_tmp(m))/(1 - qc_tmp(m) - qi_tmp(m))/mu_d)
                   R_loaded(m) = R_local(m)*(1 - qc_tmp(m) - qi_tmp(m))
                   rhos(m) = p(m)/R_local(m)/T(m)
                   T_virt(m) = T_tmp(m) * R_local(m) / rdgas
                   cp_local(m) = cp_d*(1-qv_tmp(m)-qc_tmp(m)-qi_tmp(m)) + cp_v*qv_tmp(m) + &
                        cp_liquid*qc_tmp(m) + cp_v*qi_tmp(m)
                   q_cond(m) = qi_tmp(m) + qc_tmp(m)
                   call q_sat_single(p(m), T_Tmp(m), q_cond(m), qsats(m))

                enddo
             endif
             
          else

             ! Do dry expansion with value of R_local
             dlnTdlnp = R_loaded(k+1)/cp_local(k+1)
            pfact = exp(dlnTdlnp*log(p(k)/p(k+1)))
            T_lift = T_tmp(k+1)*pfact

             T_virt_lift = T_lift*R_local(k+1)/rdgas
             rho_lift = p(k)/T_lift/R_local(k+1)
!             if (is_master() .and. omp_get_thread_num() .eq. 0) then
!                if (k==35) then
!                   write(*,*) 'BEFORE'
!                   write(*,*) T(k), T_lift, T_virt(k), T_virt_lift
!                endif
!             endif

             if (T_virt_lift .gt. T_virt(k)*(1 + tolerance)) then
                check = .true.
                !if (is_master() .and. omp_get_thread_num() .eq. 0) write(*,*) 'k, dry', k
                call conserve_dry_enthalpy(p(k), p(k+1), T_tmp(k), T_tmp(k+1), qv_tmp(k), qv_tmp(k+1), &
                     delp(k), delp(k+1), qc_tmp(k), qc_tmp(k+1), qi_tmp(k), qi_tmp(k+1))
!

! Recalculate rhos
                do m=k,k+1
                   R_loaded(m) = Rstar*(qv_tmp(m)/mu_v + &
                        (1. - qv_tmp(m) - qi_tmp(m) - qc_tmp(m))/mu_d)
                   R_local(m) = R_loaded(m)/(1. - qc_tmp(m) - qi_tmp(m))
                   rhos(m) = p(m)/R_local(m)/T(m)
                   T_virt(m) = T_tmp(m)*R_loaded(m)/rdgas
                   cp_local(m) = cp_d*(1-qv_tmp(m)-qc_tmp(m)-qi_tmp(m)) + cp_v*qv_tmp(m) + &
                        cp_liquid*qc_tmp(m) + cp_v*qi_tmp(m)
                   q_cond(m) = qi_tmp(m) + qc_tmp(m)
                   call q_sat_single(p(m), T_Tmp(m), q_cond(m), qsats(m))
                enddo
!                if (is_master() .and. omp_get_thread_num() .eq. 0) then
!                   if (k==35) then
!                   dlnTdlnp = R_loaded(k+1)/cp_local(k+1)
!                   pfact = exp(dlnTdlnp*log(p(k)/p(k+1)))
!                   T_lift = T_tmp(k+1)*pfact
!                   T_virt_lift = T_lift*R_local(k+1)/rdgas
!                   rho_lift = p(k)/T_lift/R_local(k+1)
!                   write(*,*) 'AFTER', T_lift, T_tmp(k)
!                endif
!                endif

             endif
          endif ! if (q(k+1) .gt. qsats(k+1))

       enddo ! do k=npz-1,1,-1

       
       do k=1,npz-1

          ! Do dry lift to check for condensation
          dlnTdlnp = R_local(k+1)/cp_local(k+1)
          pfact = exp(dlnTdlnp*log(p(k)/p(k+1)))
          T_lift = T_tmp(k+1)*pfact

          call q_sat_single(p(k), T_lift, qi_tmp(k+1) + qc_tmp(k+1), qsatlift)
! If saturated, gradient should be moist gradient
!          if (is_master() .and. omp_get_thread_num().eq. 0)  then
!             write(*,*) k, p(k) ,qsats(k), qsatlift
!             write(*,*) 'qv_tmp(k)', qv_tmp(k)
!             write(*,*) 'qsats(k)', qsats(k)
!             write(*,*) 'qsatlift', qsatlift
!             write(*,*) k, qv_tmp(k)/qsats(k)
!             write(*,*) T(k)
!             call sat_vp(T(k), temp)
!             write(*,*) temp, p(k)
!             conv = .true.
!          endif

          if ((qv_tmp(k+1) .gt. qsats(k+1) .or. qv_tmp(k) .gt. qsats(k)) .and.&
               (qsatlift < 1.e19 .and. qsats(k+1) < 1.e19 .and. qsats(k) .lt. 1.e19)) then
             if (qsats(k+1) .lt. 0) then
                write(*,*) 'QSATS<0', qsats(k+1)
             endif
             
             ! get moist gradient
             call gradient(p(k+1), T_tmp(k+1), dlnTdlnp,temp)
             grads(k+1) = dlnTdlnp
             ! Calculate lifted temperature from adiabatic expansion
             pfact = exp(dlnTdlnp*log(p(k)/p(k+1)))
             T_lift = T_tmp(k+1)*pfact
             T_lifts(k) = T_lift
             ! Calculate q here as proportion of gas phase rv/(1+rv)
             call r_sat_single(p(k), T_lift, q_lift)
             q_lift = q_lift/(1 + q_lift)
             q_lifts(k) = q_lift
             !call q_sat_single(p(k), T_lift, 0.0, q_lift)
             R_lift = Rstar*(q_lift/mu_v + (1. - q_lift)/mu_d)
             R_lifts(k) = R_lift
             T_virt_lift = T_lift * R_lift/rdgas
             T_virt_lifts(k) = T_virt_lift
!             if (is_master() .and. omp_get_thread_num() .eq. 0) then
!                write(*,*) 'IN CONVECTION'
!                   write(*,*) k, T_virt(k), T_virt_lift
!                   write(*,*) k, T_tmp(k), T_lift
!                   conv = .true.
!             endif
             
             if (T_virt_lift .gt. T_virt(k)*(1 + tolerance)) then
                check = .true.
! Routine does pairwise moist adjustment
               ! if (is_master() .and. omp_get_thread_num().eq. 0)  write(*,*) 'moist conv, k', k, n
                call conserve_moist_enthalpy(p(k), p(k+1), T_tmp(k), T_tmp(k+1), qv_tmp(k), qv_tmp(k+1), &
                     delp(k), delp(k+1), qc_tmp(k), qc_tmp(k+1), qi_tmp(k), qi_tmp(k+1))
                
! Recalculate the new rho, qsat etc
!                if (is_master() .and. omp_get_thread_num().eq. 0)  then
!                   write(*,*) 'CONVECTION', n
!                   call gradient(p(k+1), T_tmp(k+1), dlnTdlnp)
!                   write(*,*) k, dlnTdlnp, (log(T_tmp(k+1)) - log(T_Tmp(k)))/(log(p(k+1)) - log(p(k)))
!                   conv = .true.
!                endif
                do m=k,k+1
                   R_local(m) = Rstar*(qv_tmp(m)/mu_v/(1 - qc_tmp(m) - qi_tmp(m)) + &
                        (1 - qv_tmp(m) - qi_tmp(m) - qc_tmp(m))/(1 - qc_tmp(m) - qi_tmp(m))/mu_d)
                   R_loaded(m) = R_local(m)*(1 - qc_tmp(m) - qi_tmp(m))
                   rhos(m) = p(m)/R_local(m)/T(m)
                   T_virt(m) = T_tmp(m) * R_local(m) / rdgas
                   cp_local(m) = cp_d*(1-qv_tmp(m)-qc_tmp(m)-qi_tmp(m)) + cp_v*qv_tmp(m) + &
                        cp_liquid*qc_tmp(m) + cp_v*qi_tmp(m)
                   q_cond(m) = qi_tmp(m) + qc_tmp(m)
                   call q_sat_single(p(m), T_Tmp(m), q_cond(m), qsats(m))

                enddo
             endif
             
          else

             ! Do dry expansion with value of R_local
             dlnTdlnp = R_loaded(k+1)/cp_local(k+1)
            pfact = exp(dlnTdlnp*log(p(k)/p(k+1)))
            T_lift = T_tmp(k+1)*pfact

             T_virt_lift = T_lift*R_local(k+1)/rdgas
             rho_lift = p(k)/T_lift/R_local(k+1)
!             if (is_master() .and. omp_get_thread_num() .eq. 0) then
!                if (k==35) then
!                   write(*,*) 'BEFORE'
!                   write(*,*) T(k), T_lift, T_virt(k), T_virt_lift
!                endif
!             endif

             if (T_virt_lift .gt. T_virt(k)*(1 + tolerance)) then
                check = .true.
                
                !if (is_master() .and. omp_get_thread_num() .eq. 0) write(*,*) 'k, dry', k, n
                call conserve_dry_enthalpy(p(k), p(k+1), T_tmp(k), T_tmp(k+1), qv_tmp(k), qv_tmp(k+1), &
                     delp(k), delp(k+1), qc_tmp(k), qc_tmp(k+1), qi_tmp(k), qi_tmp(k+1))
!

! Recalculate rhos
                do m=k,k+1
                   R_loaded(m) = Rstar*(qv_tmp(m)/mu_v + &
                        (1. - qv_tmp(m) - qi_tmp(m) - qc_tmp(m))/mu_d)
                   R_local(m) = R_loaded(m)/(1. - qc_tmp(m) - qi_tmp(m))
                   rhos(m) = p(m)/R_local(m)/T(m)
                   T_virt(m) = T_tmp(m)*R_loaded(m)/rdgas
                   cp_local(m) = cp_d*(1-qv_tmp(m)-qc_tmp(m)-qi_tmp(m)) + cp_v*qv_tmp(m) + &
                        cp_liquid*qc_tmp(m) + cp_v*qi_tmp(m)
                   q_cond(m) = qi_tmp(m) + qc_tmp(m)
                   call q_sat_single(p(m), T_Tmp(m), q_cond(m), qsats(m))
                enddo
!                if (is_master() .and. omp_get_thread_num() .eq. 0) then
!                   if (k==35) then
!                   dlnTdlnp = R_loaded(k+1)/cp_local(k+1)
!                   pfact = exp(dlnTdlnp*log(p(k)/p(k+1)))
!                   T_lift = T_tmp(k+1)*pfact
!                   T_virt_lift = T_lift*R_local(k+1)/rdgas
!                   rho_lift = p(k)/T_lift/R_local(k+1)
!                   write(*,*) 'AFTER', T_lift, T_tmp(k)
!                endif
!                endif

             endif
          endif ! if (q(k+1) .gt. qsats(k+1))

       enddo ! do k=npz-1,1,-1

       n = n+1
    enddo ! do n=1,N_iter

    !if (is_master() .and. omp_get_threaD_num() .eq. 0) write (*,*) 'n', n
!    if (is_master().and. omp_get_thread_num() .eq. 0 .and. n .gt. 3) then
!       do k=1,npz-1
!          write(*,*) '-----------------------------------------------------'
!          write(*,*) k, p(k), T_tmp(k)
!          write(*,*) 'q', qv_tmp(k)
!          write(*,*) 'theory grad', grads(k+1)
!          write(*,*) 'real_grad', (log(T_tmp(k+1)) - log(T_Tmp(k)))/(log(p(k+1)) - log(p(k))), (log(T_virt(k+1)) - log(T_virt(k)))/(log(p(k+1)) - log(p(k)))
!          write(*,*) 'Tlift', T_tmp(k+1)*(p(k)/p(k+1))**(grads(k+1)), T_lifts(k)
!          write(*,*) "T virt", T_virt(k)
!          write(*,*)  'T_virt_lift', T_tmp(k+1)*(p(k)/p(k+1))**(grads(k+1))*(q_lifts(k)*rvgas/rdgas + (1-q_lifts(k))), T_lifts(k)*(q_lifts(k)*rvgas/rdgas + (1-q_lifts(k)))
!          write(*,*) 'R_lift', R_lifts(k), (q_lifts(k)*rvgas + (1-q_lifts(k))*rdgas)
!          write(*,*) 'T_virt_lifts(k)', T_virt_lifts(k)
!          write(*,*) 'Relative hum', qv_tmp(k)/qsats(k)
!          call gradient(p(k), T_tmp(k), grads(k), temp)
!          write(*,*) 'qcrit', 1/temp/(1. - mu_d/mu_v)
!          write(*,*) '-----------------------------------------------------'
!       enddo
       
!    endif
    
!    if (is_master() .and. omp_get_thread_num().eq. 0)  then
!       do k=1,npz
!          !dlnTdlnp = R_loaded(k+1)/cp_local(k+1)
!          call gradient(p(k+1), T_tmp(k+1), dlnTdlnp)
!          write(*,*) k, dlnTdlnp, (log(T_tmp(k+1)) - log(T_Tmp(k)))/(log(p(k+1)) - log(p(k))), R_loaded(k+1)/cp_local(k+1)
!          write(*,*) qv_tmp(k), qsats(k), q_cond(k)
!       
!    enddo
! endif

       
    t_dt_ding = T_tmp - T
    q_dt_ding(:) = qv_tmp - q(:,iv)
    q_ice = qi_tmp
    q_liq = qc_tmp

    call q_sat(p, T_tmp, qc_tmp + qi_tmp, qsats)
    do k=1,npz
       if (qv_tmp(k) .lt. 0) then
          write(*,*) 'VAPOUR LESS THAN 0!!', k, q(k,iv), qv_tmp(k)
          stop
       endif
    enddo
    
!    ! ENTHALPY TEST
!    h_aft=0.0
!    m_aft = 0.0
!    do k=1,npz
!       call calc_enthalpy(T_tmp(k), qv_tmp(k), q_liq(k), q_ice(k), h_temp)
!       h_aft = h_aft + h_temp*delp(k)
!       m_aft = m_aft + (qv_tmp(k)+q_ice(k) + q_liq(k))*delp(k)
!    enddo
!    if (conv_triggered) write(*,*) 'conservation',  abs((h_aft - h_before)/h_before), abs((m_bef - m_aft)/m_bef)!!'ICE', maxval(qi_tmp), maxval(q_ice)


  end subroutine ding_adjust

  subroutine conserve_dry_enthalpy(p1, p2, T1, T2, q1, q2, dp1, dp2, qc1, qc2, qi1, qi2)
    !==========================================================================
    ! Description
    !==========================================================================
    ! Calculates adjusted state with a Newton iteration by conserving dry enthalpy
    ! in the absence of condensation on adiabatic lifting

    !==========================================================================
    ! Input variables
    !==========================================================================    
    real, intent(in) :: p1, p2 ! Pressure of upper/lower layer
    real, intent(inout) :: T1, T2 ! Temperature of upper/lower layer
    real, intent(inout) :: q1, q2,qc1,qc2,qi1,qi2 ! Specific humidity of lower/upper layer
    real, intent(in)    :: dp1, dp2 ! Pressure thicknesses of the layers
    
    !==========================================================================
    ! Local variables
    !==========================================================================    
    real :: k_old, T2_guess, T1_new,k_new, k1, k2

    real :: q1_new, q2_new,k_diff, qc1_new, qc2_new, qi1_new, qi2_new
    real :: temp1,temp2,temp4, temp
    integer ::  n, m
    
    !==========================================================================
    ! Main body
    !==========================================================================    

    ! Enthalpy of the old state
    
    call calc_enthalpy(T1, q1, qc1,qi1,k1)
    call calc_enthalpy(T2, q2, qc2,qi2,k2)

    k_old = k1*dp1 + k2*dp2
    
    ! Calculate constant q of the after state by conserving mass
    q1_new = (q1*dp1 + q2*dp2)/(dp1 + dp2)
    q2_new = q1_new

    qc1_new = (qc1*dp1 + qc2*dp2)/(dp1+dp2)
    qc2_new = qc1_new

    qi1_new = (qi1*dp1 + qi2*dp2)/(dp1 + dp2)
    qi2_new = qi1_new
    

    T2_guess = T2
    n = 0
    temp= 1.e13
    do while ((abs(temp) .gt. precision) .and. (n .lt. n_iter_newt))
       ! Get R/cp of the lower layer
       
       call k_adj_dry(T2_guess, p1, p2, dp1, dp2, q1_new,  qc1_new,qi1_new, T1_new, k_new)
       call dk_adj_dry(T2_guess, p1,p2,dp1,dp2,q1_new, qc1_new,qi1_new, k_diff)

       temp = (k_new - k_old)/k_diff
       
       if (n .gt. n_iter_newt - 10) then
          
          write(*,*) 'NEAR END OF NEWT'
          call k_adj_dry(T2_guess, p1, p2, dp1, dp2, q1_new,  qc1_new, qi1_new, T1_new, k_new, .true.)
          write(*,*) T2_guess, temp, k_new, k_old, k_diff, (k_new-k_old), (k_new-k_old)/k_diff
       endif
       
       T2_guess = T2_guess  - temp
       n = n+1
       
    enddo

    if (abs(temp) .gt. precision) then
       write(*,*) 'max iteration reached, DRY'
       write(*,*) temp4, T2_guess

       do m=1,20
          call k_adj_dry(T2_guess+temp, p1,p2,dp1,dp2,q1_new,qc1_new,qi1_new,temp1, temp2, .true.)
          write(*,*) 'm=', m, temp1 - k_old, k_new-k_old, T2_guess-1.+m*0.2
       enddo

       return
    endif
    
    call k_adj_dry(T2_guess, p1, p2, dp1, dp2, q1_new, qc1_new, qi1_new,T1_new, k_new)

!if ((abs(T2_guess-T2) .gt. 20.0) .or. (abs(T1-T1_new).gt.20.0) ) return
    q1 = q1_new
    q2 = q2_new

    qc1 = qc1_new
    qc2 = qc2_new
    qi1 = qi1_new
    qi2 = qi2_new
    
    T2 = T2_guess
    T1 = T1_new
  end subroutine conserve_dry_enthalpy

  subroutine k_adj_dry(T2_guess, p1, p2, dp1, dp2, q1, qc1, qi1,T1_new, k_new, printme)
    real, intent(in) :: T2_guess, p1, p2, dp1, dp2, q1, qc1, qi1
    real, intent(out) :: T1_new, k_new
    logical, intent(in), optional :: printme

    real :: R_new, cpv_2,cpc_2,cpi_2, k1, k2, pfact,cp_new, f_ice
    
    R_new = Rstar*(q1/mu_v + (1-q1-qc1-qi1)/mu_d)

    if (present(printme)) write(*,*) 'k_adj_dry, R_new', R_new

    call get_f_ice(T2_guess, f_ice)
    
    cpv_2 = cp_v
    cpc_2 = cp_liquid
    cpi_2 = cp_v

    cp_new = q1*cpv_2 + qc1*cpc_2*(1-f_ice) + qi1*cpi_2*f_ice + (1 - q1 - qc1 - qi1)*cp_d
    if (present(printme)) write(*,*) 'k_adj_dry, cp_new', cp_new
    pfact = exp(R_new/cp_new * log(p1/p2))
    T1_new = T2_guess*pfact

    call calc_enthalpy(T1_new, q1,qc1, qi1, k1)
    call calc_enthalpy(T2_guess, q1,qc1, qi1, k2)

    if (present(printme)) write(*,*) 'k_adj_dry, k1, k2', k1, k2
    
    k_new = k1*dp1 + k2*dp2
    
    if (present(printme)) write(*,*) 'k_adj_dry', k_new
  end subroutine k_adj_dry

  subroutine dk_adj_dry(T2_guess, p1, p2, dp1, dp2,q1, qc1,qi1, dk_adj)
    real, intent(in) :: T2_guess, p1, p2, dp1, dp2,q1,qc1, qi1
    real, intent(out) :: dk_adj

    real :: t, kp, km
    real :: eps = 1.e-8

    call k_adj_dry(T2_guess+eps/2, p1, p2, dp1, dp2, q1, qc1, qi1, t, kp)
    call k_adj_dry(T2_guess-eps/2, p1, p2, dp1, dp2, q1, qc1, qi1, t, km)

    dk_adj = (kp - km)/eps
  end subroutine dk_adj_dry
  
  subroutine conserve_moist_enthalpy(p1, p2, T1, T2, q1, q2, dp1,dp2, qc_1, qc_2, qi_1, qi_2)
                                     
    !==========================================================================
    ! Description
    !==========================================================================
    ! Calculates adjusted state with a Newton iteration by conserving moist enthalpy

    !==========================================================================
    ! Input variables
    !==========================================================================    
    real, intent(in) :: p1, p2 ! Pressure of upper/lower layer
    real, intent(inout) :: T1, T2 ! Temperature of upper/lower layer
    real, intent(inout) :: q1, q2 ! Specific humidity of lower/upper layer
    real, intent(in)    :: dp1, dp2 ! Pressure thicknesses of the layers
    real, intent(inout)   :: qc_1, qc_2 ! Vapour and condensate in adjusted state
    real, intent(inout) :: qi_1, qi_2 ! Ice in adjusted state
    
    !==========================================================================
    ! Local variables
    !==========================================================================    
    real :: k_old, T2_guess,  q1_new, q2_new, k_new

    real :: k_diff, dT, T1_new, grad
    real :: k1_old, k2_old,mq_old, mq_new, qt_10, qt_20
    real :: qc_1_new, qc_2_new, qi_1_new, qi_2_new, out, tvirt2, tvirt1,tvirtlift, tlift, qlift, qsat1,qsat2
    integer :: n

    ! Calculate moist enthalpy of old state
    !write(*,*) 'IN MOIST ADJUSTMENT'
    call calc_enthalpy(T1, q1, qc_1, qi_1, k1_old)
    call calc_enthalpy(T2, q2, qc_2, qi_2, k2_old)

    k_old = k1_old*dp1 + k2_old*dp2
    mq_old = dp1*(q1+qc_1 + qi_1) + dp2*(q2+qc_2 + qi_2)
!!$    write(*,*) 'IN CONSMOISTENTH'
!!$    write(*,*) q1, qc_1, qi_1, dp1
!!$    write(*,*) q2, qc_2, qi_2, dp2
    qt_10 = q1 + qc_1 + qi_1
    qt_20 = q2 + qc_2 + qi_2
    
    
    n = 0
    dT = 1.e13 ! Arbitrary high number
    T2_guess = (T2+T1)*0.5
    !write(*,*) 'before newt', T1, T2
    do while ((abs(dT) .gt. precision) .and. (n .lt. n_iter_newt))
       !write(*,*) 'n = ', n
       call k_adj(T2_guess, p1, p2, dp1, dp2,qt_10,qt_20, q1_new, q2_new, &
            qc_1_new, qc_2_new, qi_1_new, qi_2_new, T1_new, k_new)
       call dk_adj(T2_guess, p1,p2,dp1,dp2, qt_10, qt_20, k_diff)

       !write(*,*) 'IN CONSMOISTENTH', T1,T2,T1_new, T2_guess, qi_1_new, qi_2_new
       dT = (k_new - k_old)/k_diff
       dT = min(dT, 5.)
       dT = max(dT, -5.)
       T2_guess = T2_guess  - 0.5*dT

       n = n+1

       if (T1 .ne. T1 .or. T2 .ne. T2 .or. T1_new .ne. T1_new .or. T2_guess .ne. T2_guess) then
          write(*,*) 'NAN', T1, T2, T1_new, T2_guess
          write(*,*) dT, k_new, k_old, k_diff
       endif
       
       if (n .gt. n_iter_newt-5) then
          write(*,*) 'NEAR MAX NEWT, dT, T2guess, T1, T2', 0.5*dT, T2_guess, T1, T2
          write(*,*) q1, q2, qi_1, qi_2
          call q_sat_single(p1, T1, qi_1 + qc_1, out)
          write(*,*) 'qsat1', out
          call sat_vp(T1, out)
          write(*,*) 'psat1, p1', out, p1
          call q_sat_single(p2, T2, qi_2 + qc_2, out)
          write(*,*) 'qsat2', out
          call sat_vp(T2, out)
          write(*,*) 'psat2, p2', out, p2
       endif
       
       if (n .eq. n_iter_newt-1) then
          write(*,*) 'Max newton iteration reached, MOIST'
          return
       endif
       
    enddo
    
    call k_adj(T2_guess, p1, p2, dp1, dp2, qt_10,qt_20,q1_new, q2_new, &
         qc_1_new, qc_2_new, qi_1_new, qi_2_new, T1_new, k_new, .true.)
    if (abs(dT) .gt. precision) then
       write(*,*) 'MAX newton, MOIST', T1, T2
       return
    endif

    if (q1_new .lt. 0. .or. q2_new .lt. 0) then
       write(*,*) 'LESS THAN 0', q1_new, q2_new, q1, q2
       write(*,*) T1, T2,T1_new, T2_guess
    endif
    
    q1 = q1_new
    q2 = q2_new
    qc_1 = qc_1_new
    qc_2=  qc_2_new
    qi_1 = qi_1_new
    qi_2 = qi_2_new
    
    T2 = T2_guess
    T1 = T1_new
    
    mq_new = (q1 +qc_1 + qi_1)*dp1 + (q2 +qc_2+qi_2)*dp2

!    if (is_master() .and. omp_get_thread_num() .eq. 0) then
!       call gradient(p2, T2_guess, grad)
!       write(*,*) 'grad, dlntdlnp', grad, (log(T2) - log(T1))/(log(p2) - log(p1))
!       tvirt2 = T2*(rvgas*q2 + (1 - q2)*rdgas)/rdgas
!       tvirt1 = T1*(rvgas*q1 + (1 - q1)*rdgas)/rdgas
!       tlift = T2*exp(grad*log(p1/p2))
!       call r_sat_single(p1, tlift, qlift)
!       qlift = qlift/(1 + qlift)
!       tvirtlift = tlift*(1. + (rvgas/rdgas - 1)*qlift)
!       write(*,*) 'T2, T1', t2, t1
!       write(*,*) 'tvirt2, tvirt1', tvirt2,tvirt1
!       write(*,*) 'tlift, tvirtlift', tlift, tvirtlift
!       call q_sat_single(p1, T1, qi_1 + qc_1, qsat1)
!       call q_sat_single(p1, T2, qi_2 + qc_2, qsat2)
!       write(*,*) 'RH1, RH2', q1/qsat1, q2/qsat2
!       write(*,*) 'q2, qlift', q2, qlift, q2/qlift
!    endif
  end subroutine conserve_moist_enthalpy

  subroutine k_adj(T2_guess, p1, p2, dp1, dp2, q1_0,q2_0,qv_1,qv_2,&
       qc_1, qc_2,qi_1, qi_2, T1_new,k_new, printme)
    real, intent(in):: T2_guess, p1, p2, dp1, dp2, q1_0,q2_0
    real, intent(out) ::  qc_1, qc_2, k_new, T1_new,qv_1, qv_2, qi_1, qi_2
    logical, intent(in), optional :: printme

    real :: dlnTdlnp, rc_1, rc_2, rv_1, rv_2
    real :: k1, k2, eta, f1, pfact, f_ice, m_i
    
    call gradient(p2, T2_guess, dlnTdlnp)
    
    pfact = exp(dlnTdlnp*log(p1/p2))
    T1_new = T2_guess*pfact

    call r_sat_single(p2, T2_guess,rv_2)
    call r_sat_single(p1, T1_new,rv_1)
    
    ! Ratio of upper condensate to lower condensate mixing ratio
    eta = (1 + rv_1)/(1+ rv_2)

    ! Initial mass of vapour+condensate
    m_i = q1_0*dp1 + q2_0*dp2

    ! Calculate rc_2 according to Ding and Pierrehumbert (note our lower + upper layer)
    ! indices are swapped wrt this paper
    rc_2 = (m_i*(1 + rv_2) - dp2*rv_2 - dp1*rv_1/eta)/(dp1 + dp2 - m_i)

    if (rc_2 .lt. 0) then
!       if (present(printme)) then
!          if (printme .and. is_master() .and. omp_get_thread_num() .eq. 0) write(*,*) 'NOT ENOUGH VAP FOR CONDENSATE'
!       endif
       ! Condensate is negative (not enough initial mass of water)
       qv_1 = rv_1/(1 + rv_1)
       qv_2 = rv_2/(1 + rv_2)

       f1 = m_i / (qv_2*dp2 + qv_1*dp1)

       qv_1 = qv_1 * f1
       qv_2 = qv_2 * f1

       qc_1 = 0.
       qc_2 = 0.
       qi_1 = 0.
       qi_2 = 0.
       
    else
       ! Condensate is positive
       call get_f_ice(T1_new, f_ice)
       

       rc_1 = rc_2 * eta

       qv_1 = rv_1/(1+rv_1 + rc_1)
       qc_1 = rc_1/(1+rv_1 + rc_1) * (1-f_ice)
       qi_1 = rc_1/(1+rv_1 + rc_1) * f_ice

       call get_f_ice(T2_guess, f_ice)

       qv_2 = rv_2/(1+rv_2 + rc_2)
       qc_2 = rc_2/(1+rv_2 + rc_2)*(1 - f_ice)
       qi_2 = rc_2/(1+rv_2 + rc_2)*f_ice

    endif

    ! Find enthalpy of the new structure
    ! Dry component = cp*T, moist component = tabulated

    call calc_enthalpy(T1_new, qv_1, qc_1, qi_1, k1)
    call calc_enthalpy(T2_guess, qv_2, qc_2,qi_2, k2)

    k_new = k1*dp1 + k2*dp2

    if (k_new .ne. k_new) then
       write(*,*) 'NAN', T2_guess, T1_new, qv_1, qv_2, qc_1, qc_2
       write(*,*) k1, k2, rv_1, rv_2, dlnTdlnp
    endif
    
  end subroutine k_adj

    
  subroutine large_scale_cond(p, T, qv, qc, qi, q_dt_lsc, t_dt_lsc)
    real, intent(in) :: p(:)
    real, intent(in) :: T(:), qv(:)
    real, intent(inout) :: qc(:), qi(:)
    real, intent(inout) :: q_dt_lsc(:)
    real, intent(out) :: t_dt_lsc(:)
    ! Local variables
    integer:: k, n
    integer :: npz

    real :: qsat,dT , qt, q_guess, T_guess, f_ice_guess
    real :: kold, knew, Tgp, qgp, Tgm,qgm, kp, km, diff, qc_guess, qi_guess
    real :: q_cond
    real :: T_tmp(size(p)), qv_tmp(size(p)), qc_tmp(size(p)), qi_tmp(size(p))

    real :: mbef, maft
    ! Take into account tendencies from Ding Convection
    T_tmp = T
    qv_tmp = qv 
    qc_tmp = qc !q(:,ic) + q_dt(:,ic)
    qi_tmp = qi
    
    npz = size(p)

    do k=1,npz

       qt = qc_tmp(k) + qv_tmp(k) + qi_tmp(k)
       mbef = qt
       q_cond = qc_tmp(k) + qi_tmp(k)
       !call sat_vp(T_tmp(k), psat)
       
       call q_sat_single(p(k), T_tmp(k), q_cond,qsat)

       if (qv_tmp(k) .gt. qsat) then
       !if (x*p(k) .gt. psat) then
          
          ! Do in terms of the enthalpy of beginning vs final state
          call calc_enthalpy(T_tmp(k), qv_tmp(k), qc_tmp(k),qi_tmp(k), kold)
          dT = 0.0
          T_guess = T_tmp(k)
          q_guess = qsat
          n = 0
          dT = 1.e13
          do while ( (abs(dT) .gt. precision) .and. (n .lt. n_iter_newt))
             call mv_mvmd(p(k), T_guess, q_guess)
             q_guess = q_guess/(1-q_guess) * (1 - qt)
             !call q_sat_single(p(k), T_guess,q_cond, q_guess)
             call get_f_ice(T_guess, f_ice_guess)

             qc_guess = (qt - q_guess)*(1-f_ice_guess)
             qi_guess = (qt - q_guess)*f_ice_guess
             call calc_enthalpy(T_guess, q_guess, qc_guess, qi_guess, knew)

             
             Tgp = T_guess+1.e-8
             call mv_mvmd(p(k), Tgp, qgp)
             qgp = qgp/(1-qgp) * (1 - qt)

!             call q_sat_single(p(k), Tgp, q_cond,qgp)
             call get_f_ice(Tgp, f_ice_guess)
             
             qc_guess = (qt - qgp)*(1-f_ice_guess)
             qi_guess = (qt - qgp)*f_ice_guess
             call calc_enthalpy(Tgp, qgp,qc_guess,qi_guess, kp)

             Tgm = T_guess - 1.e-8
             call mv_mvmd(p(k), Tgm, qgm)
             qgm = qgm/(1-qgm) * (1 - qt)

             !call q_sat_single(p(k), Tgm, q_cond, qgm)
             call get_f_ice(Tgm, f_ice_guess)
             
             qc_guess = (qt - qgm)*(1-f_ice_guess)
             qi_guess = (qt - qgm)*f_ice_guess
             
             call calc_enthalpy(Tgm, qgm,qc_guess, qi_guess,km)

             diff = (kp-km)/(2.e-8)
             dT = (kold - knew)/diff
             
             T_guess = T_guess + dT
             n = n+1
          enddo

          ! Update T(k) with the value from large scale condensation
          if (abs(dT) .lt. precision) then
             T_tmp(k) = T_guess
             !call q_sat_single(p(k), T_guess,q_cond, q_guess)
             call mv_mvmd(p(k), T_guess, q_guess)
             q_guess = q_guess/(1-q_guess)* (1 - qt)
             !qgp = qgp/(1-qgp) * (1 - qt)

             call get_f_ice(T_guess, f_ice_guess)

             
             qc_tmp(k) = (qt - q_guess)*(1-f_ice_guess)
             
             qi_tmp(k) = (qt - q_guess)*f_ice_guess

             if (q_guess .lt. 0 ) then
                write(*,*) 'q<0 in cond: q_guess, qi_tmp(k), k'
                write(*,*)  q_guess, qi_tmp(k), k
                write(*,*) 'T_guess, qsat, qv_tmp(k)'
                write(*,*) T_guess, qsat, qv_tmp(k)
                write(*,*) qv(k), qsat
                stop
             endif
             
             qv_tmp(k) = q_guess   
          else
             write(*,*) 'MAX COND ITERATION REACHED'
          endif
          
          
       endif

       maft = qi_tmp(k) + qc_tmp(k) + qv_tmp(k)
       if (abs(mbef - maft)/mbef .gt. 1.e-10) then
          write(*,*) 'k, mbef, maft'
          write(*,*) k, mbef, maft
       endif
       
    enddo

    ! Update tendencies with large scale condensation
    t_dt_lsc = T_tmp - T
    q_dt_lsc = qv_tmp - qv
    !q_dt(:,ic) = qi_tmp + qc_tmp - (qc + qi)
    qi = qi_tmp
    qc = qc_tmp
    !q_dt(:,ic) = qc_tmp - q(:,ic)
    
  end subroutine large_scale_cond
  
  subroutine calc_enthalpy(T, q_v, q_c, q_i, kout)
    real, intent(in)  :: T, q_v, q_c, q_i ! Temp, vapour and condensate ratio
    real, intent(out) :: kout ! output enthalpy

    real :: h_d, h_v, h_c, h_i

    h_d = cp_d*T
    h_v = hv_tp + cp_v*(T-T_tp)
    h_i = h_v - L_sub
    h_c = hl_tp + cp_liquid*(T-T_tp)

    ! Per unit mass, multiply by dp out of routine!
    kout = h_d*(1 - q_v - q_c - q_i) + h_c*q_c + h_v*q_v + h_i*q_i
  end subroutine calc_enthalpy

  subroutine calc_enthalpy_cond(T, q_l, q_i, kout)
    real, intent(in) :: T, q_l, q_i ! Temp, liquid and ice
    real, intent(out) :: kout ! Enthalpy of these components

    real  :: h_l, h_i

    h_l = hl_tp + cp_liquid*(T-T_tp)
    h_i = hv_tp + cp_v*(T-T_tp) - L_sub

    kout = h_l*q_l + h_i*q_i
  end subroutine calc_enthalpy_cond

  subroutine calc_enthalpy_mass(T, m_d, m_v, m_l, m_i, kout, printit)
    ! Mass weighted enthalpy, use when mass of layer will be changing!
    real, intent(in) :: T, m_d, m_v, m_l, m_i
    real, intent(out) :: kout
    logical, optional, intent(in) :: printit
    real :: h_d, h_v, h_i, h_l, h_v_specific, h_i_specific

    h_d = cp_d*T*m_d
    !h_v = (hv_tp + cp_v*(T-T_tp))*m_v
    h_v_specific = hv_tp + cp_v*(T - T_tp)
    h_v = h_v_specific*m_v
    
    h_i_specific= h_v_specific - L_sub
    h_i = h_i_specific*m_i
    h_l = (hl_tp + cp_liquid*(T-T_tp))*m_l

    if (present(printit)) then
       if (printit) then
          
          write(*,*) 'DRY, VAPOUR, LIQUID, ICE'
          write(*,*) h_d, h_v, h_l, h_i
       endif
    endif
    
   kout = h_d + h_v + h_i + h_l
  end subroutine calc_enthalpy_mass
  
  subroutine rain_out_revap(T, p, delp, qv, ql, qi, q_dt_rainout, dT, ua, va, du,dv)
    real, intent(in) :: T(:), p(:), delp(:), qv(:)
    real, intent(inout) :: q_dt_rainout(:), ql(:), qi(:)
    real, intent(out) :: dT(:)
    real, intent(in) :: ua(:), va(:)
    real, intent(inout) :: du(:), dv(:)
    real, dimension(size(p)) :: ql_tmp, qi_tmp, T_tmp, qv_tmp
    real :: small, T_cond, satq, delT, T_guess, Tgp, Tgm
    real :: h0, hcond, hgas, hnew, hp, hm, diff, md, mv, mv_new, m_cond
    real:: fi
    
    real, dimension(size(p)) :: m_ds, m_vs, m_ls, m_is, delqi, delql, delqv, umom, vmom
    real :: umom_tot, umom_aft
    integer :: npz, k, l,n, m
    
    logical :: l_evap
    ! Add on tendencies from previous parts of the condensation scheme

    npz = size(T)
    dT = 0.0
    
    do k=1,npz
       T_tmp(k) = T(k)
       ql_tmp(k) = ql(k)
       qi_tmp(k) = qi(k)
       qv_tmp(k) = qv(k)
       m_ds(k) = (1. - qv_tmp(k) - ql_tmp(k) - qi_tmp(k))*delp(k)
       m_vs(k) = qv_tmp(k)*delp(k)
       m_ls(k) = ql_tmp(k)*delp(k)
       m_is(k) = qi_tmp(k)*delp(k)
    enddo

    small = 1.e-20

    l_evap = .false.
    
!    if (is_master() .and. omp_get_thread_num() .eq. 0) then
!       write(*,*) 'REVAP'
!       do k=1,npz
!          write(*,*) k, ql_tmp(k), qi_tmp(k)
!       enddo
!    endif
!!$ =====================================================================================
!!$ Uncomment for enthalpy diagnostics
!!$ =====================================================================================    
!    hbef = 0.0
!    
!    mbef_i = 0.0
!    mbef_v = 0.0
!    mbef_l = 0.0
!    do k=1,npz
!       mbef_i = mbef_i + m_is(k)
!       mbef_l = mbef_l + m_ls(k)
!       mbef_v = mbef_v + m_vs(k)
!    
!    enddo
    

!if (is_master() ) write(*,*) shape(m_is), shape(m_vs), shape(m_ls)
!    do k=1,npz
!       !call calc_enthalpy_mass(T_tmp(k), m_ds(k), m_vs(k), m_ls(k), m_is(k), htmp)
!       call calc_enthalpy(T_tmp(k), qv_tmp(k), ql(k), qi(k), htmp)
!       hbef = hbef + htmp*delp(k)
!    enddo
!!$ ====================================================================================

    delqv = 0.0
    delql = 0.0
    delqi = 0.0
    umom = 0.0
    vmom = 0.0
    umom_tot = 0.0
    do k=1,npz
       umom_tot = umom_tot + delp(k)*ua(k)
    enddo
    
    do k=1,npz
       ! From top down, check if there is condensate in the layer. Then search for the next
       ! undersaturated layer at < 90% RH (Xianyu Tan Pure steam conv. paper)
       ! and attempt re-evaporation here
       if ( ql_tmp(k) .gt. small .or. qi_tmp(k) .gt. small) then
          ! Condensate! Attempt re-evaporation

          T_cond = T_tmp(k)
          m_cond = (ql_tmp(k) + qi_tmp(k))*delp(k)
          m_cond = m_ls(k) + m_is(k)

          fi = m_is(k)/(m_is(k) + m_ls(k))
          ! fraction of ice in this condensate
          

          do l=k,npz
             
             call q_sat_single(p(l), T_Tmp(l), ql_tmp(l) + qi_tmp(l), satq)

             if (qv_tmp(l) .lt. satq*0.9 .and. qi_tmp(l) .lt. small .and. ql_tmp(l) .lt. small) then
                ! Re-evaporate!
                ! First attempt to re-evaporate all the condensate in this layer

                ! Condensate falling out of layer
                l_evap = .true.
                delqi(k) = -qi(k)
                delql(k) = -ql(k)
                ql(k) = 0.0
                qi(k) = 0.0

                md = m_ds(l)
                mv = m_vs(l)
                mv_new  = mv + m_cond

                ! Initial enthalpy of gas component
                call calc_enthalpy_mass(T_tmp(l), md, mv, m_ls(l), m_is(l), hgas)
                ! Initial enthalpy of the condensate
                call calc_enthalpy_mass(T_tmp(k), 0.0, 0.0, m_cond*(1-fi), m_cond*fi, hcond)

                ! Total starting enthalpy
                h0 = hgas + hcond

                delT = 1.e13
                T_guess = T_tmp(l)
                ! Newton iteration!
                n = 1
                do while((abs(delT) .gt. precision) .and. n .lt. n_iter_newt)
                   ! Calculate enthalpy at guess temperature
                   ! Mass of dry substance stays constant!
                   call calc_enthalpy_mass(T_guess, md, mv_new, 0.0, 0.0, hnew)

                   ! Calculate enthalpy at perturbed upwards temp
                   Tgp = T_guess + 1.e-8
                   call calc_enthalpy_mass(Tgp, md, mv_new, 0.0, 0.0, hp)
                   
                   ! Calculate enthalpy at perturbed downwards temp
                   Tgm = T_guess - 1.e-8
                   call calc_enthalpy_mass(Tgm, md, mv_new, 0.0, 0.0, hm)
                   
                   ! Find perturbed temperature by doing newton iteration
                   diff = (hp - hm)/2.e-8 ! Derivative
                   delT = (h0 - hnew)/diff

                   T_guess = T_guess + delT
                   n = n+1
                enddo

                if (abs(delT) .lt. precision) then
                   ! Check to see if qv now over q_sat, if so, redo iteration to get to
                   ! Saturation and move down to next layer

                   call q_sat_single(p(l), T_guess, 0.0, satq)

                   if (mv_new/(mv_new + md) .lt. 0.9*satq) then
                      ! All condensate succesfully evaporated

                      delqv(l) = delqv(l) + m_cond/delp(l)
                      
                      T_tmp(l) = T_guess
                      dT(l) = T_tmp(l) - T(l)
                      qv_tmp(l) = qv_tmp(l) + m_cond/delp(l)
                      m_vs(l) = m_vs(l) + m_cond
                      q_dt_rainout(l) = qv_tmp(l) - qv(l)

                      umom(l) = umom(l) + m_cond*ua(k)
                      vmom(l) = vmom(l) + m_cond*va(k)

                      umom(k) = umom(k) - m_cond*ua(k)
                      
                      m_is(k) = 0.0
                      m_ls(k) = 0.0

                     ! Move onto next layer with condensate
                      exit
                      
                   else if (l .eq. npz) then
                      write(*,*) 'Cannot evaporate in lower atmosphere!!'

                      do m=1,npz
                         call q_sat(p, T_tmp, ql_tmp + qi_tmp, dT)
                         write(*,*) m, T_tmp(m), qv_tmp(m), dT(m), qv_tmp(m)/dT(m)
                      enddo
                      stop
                   endif
                   ! Below is scheme to re-evaporate up to saturation. More complicated because it requires a newton
                   ! iteration to determine T change that will bring up to saturation. Seems more complicated than necesssary
                   ! in an already complicated scheme. Leaving here in case useful in future
                   
!                    else
!                       l_evap=.true.
!                       write(*,*) 're-evap up to saturation'
!                       write(*,*) '=================================================='
!                       !write(*,*) 'MAX newton'
              
!                       ! Evaporate up to saturation, subtract condensate, then
!                       ! move to next undersaturated layer
!                       delT = 1.e13
!                       T_guess = T_tmp(l)
!                       n = 1
!                       do while((abs(delT) .gt. precision) .and. n .lt. n_iter_newt)
!                          !write(*,*) l, 'Evaporating up to saturation'
!                          ! Calculate enthalpy at guess temperature
!                          ! Mass of dry substance stays constant!
!                          call q_sat_single(p(l), T_guess, 0.0, satq)
!                          !m_evap = (satq*delp(l) - mv)/(1 - satq)
!                          m_evap = md*satq/(1-satq) - mv
!                          mv_new = mv + m_evap
!                          call calc_enthalpy_cond(T_cond, m_evap*(1-fi)/delp(k), m_evap*fi/delp(k), hcond)
!                          call calc_enthalpy_mass(T_cond, 0.0, 0.0, m_evap*(1-fi), m_evap*fi, hcond)
!                          h0 = hgas + hcond
!                          !h0 = hgas*delp(l) + hcond*delp(k)

!                          call calc_enthalpy(T_guess, mv_new/(mv_new+md), 0.0,0.0, hnew)
!                          call calc_enthalpy_mass(T_guess, md, mv_new, 0.0, 0.0, hnew)
!                          !hnew = hnew*delp(l)
                         
!                          f = hnew - h0
                        
!                          ! Calculate enthalpy at perturbed upwards temp
!                          Tgp = T_guess + 1.e-8
!                          call q_sat_single(p(l), Tgp, 0.0, satq)
!                          !m_evap = (satq*delp(l) - mv)/(1 - satq)
!                          m_evap = md*satq/(1-satq) - mv
                         
!                          mv_new = mv + m_evap
!                          call calc_enthalpy_cond(T_cond, m_evap*(1-fi)/delp(k), m_evap*fi/delp(k), hcond)
!                          call calc_enthalpy_mass(T_cond, 0.0, 0.0, m_evap*(1-fi), m_evap*fi, hcond)
!                          !h0 = hgas*delp(l) + hcond*delp(k)
!                          h0 = hgas + hcond
                         
!                          call calc_enthalpy(Tgp, mv_new/(mv_new+md), 0.0,0.0, hnew)
!                          call calc_enthalpy_mass(Tgp, md, mv_new, 0.0, 0.0, hnew)
!                          !hnew = hnew*delp(l)
!                          fp = hnew - h0

!                         ! Calculate enthalpy at perturbed downwards temp
!                          Tgm = T_guess - 1.e-8
!                          call q_sat_single(p(l), Tgm, 0.0, satq)
!                          !m_evap = (satq*delp(l) - mv)/(1 - satq)
!                          m_evap = md*satq/(1-satq) - mv
!                          mv_new = mv + m_evap
!                          call calc_enthalpy_cond(T_cond, m_evap*(1-fi)/delp(k), m_evap*fi/delp(k), hcond)
!                          call calc_enthalpy_mass(T_cond, 0.0, 0.0, m_evap*(1-fi), m_evap*fi, hcond)
!                          !h0 = hgas*delp(l) + hcond*delp(k)
!                          h0 = hgas + hcond

!                          call calc_enthalpy(Tgm, mv_new/(mv_new+md), 0.0,0.0, hnew)
!                          call calc_enthalpy_mass(Tgm, md, mv_new, 0.0, 0.0, hnew)
!                          !hnew = hnew*delp(l)
!                         fm = hnew - h0
                        
!                           ! Find perturbed temperature by doing newton iteration
!                          diff = (fp - fm)/2.e-8 ! Derivative
!                          delT = -f/diff

! !!$                         write(*,*) 'fp, fm, f, diff'
! !!$                         write(*,*) fp, fm, f, diff
! !!$                         write(*,*) 'T_guess, delT'
! !!$                         write(*,*) T_guess, delT, T_guess - T_tmp(l)
! !!$                         write(*,*) 'm_evap, satq, fi'
! !!$                         write(*,*) m_evap, satq, fi
! !!$                         
!                          T_guess = T_guess + 0.5*delT
!                          if (n .gt. n_iter_newt - 10) then
!                             write(*,*) 'Approaching max newton iters'
!                             write(*,*) n_iter_newt - n, T_guess, delT, m_evap, mv, md
!                             n = n+1
!                          endif
                         
!                       enddo

!                       if (abs(delT) .lt. precision) then
                         
!                          write(*,*) 'newton converged'
!                          call calc_enthalpy_mass(T_tmp(k), m_ds(k), m_vs(k), m_evap*(1-fi), m_evap*fi, h_k0)
!                          call calc_enthalpy_mass(T_tmp(l), m_ds(l), m_vs(l), m_ls(l), m_is(l), h_l0)

!                          call q_sat_single(p(l), T_guess, 0.0, satq)
!                          m_evap = md*satq/(1-satq) - mv!(satq*delp(l) - mv)/(1 - satq)
!                          mv_new = mv + m_evap
!                          m_vs(l) = mv_new
! !                         write(*,*) 'masses, m_cond, m_evap, mv'
! !                         write(*,*) m_cond, m_evap, mv, mv_new/(mv_new+md)
!                          hsmall = hsmall + hnew
!                          qv_tmp(l) = qv_tmp(l) + m_evap/delp(l)
!                          q_dt(l,iv) = qv_tmp(l) - qv(l)

!                          !write(*,*) 'T_guess, delT'
!                          !write(*,*) T_guess, T_guess - T_tmp(l), T_tmp(l), T_cond
!                          !write(*,*) k,l
! !                         write(*,*) 'hnew, hgas bef, hcond bef, relhchang'

!                          !                         write(*,*) hnew, hgas, hcond, hgas + hcond, abs((hnew - h0)/h0)
!                          call calc_enthalpy(T_tmp(k), qv_tmp(k), 0.0, 0.0, h_k)
!                          call calc_enthalpy_mass(T_tmp(k), m_ds(k), m_vs(k), 0.0,0.0,h_k)
!                          call calc_enthalpy(T_guess, qv_tmp(l), 0.0,0.0, h_l)
!                          call calc_enthalpy_mass(T_guess, md, mv_new, 0.0, 0.0, h_l)

!                          !write(*,*) 'BEFORE, AFTER, CHANGE'
!                          !write(*,*) h_k0, h_l0, h_k, h_l, h_k0+h_l0, h_k + h_l
!                          !write(*,*) abs((h_k + h_l-h_k0-h_l0)/(h_l0+h_k0))
!                          T_tmp(l) = T_guess
!                          dT(l) = T_tmp(l) - T(l)

! !!$                         haft = 0.0
! !!$                         do m=1,npz
! !!$                            call calc_enthalpy_mass(T_tmp(m), m_ds(m), m_vs(m), m_ls(m), m_is(m), htmp)
! !!$                            haft = haft+htmp
! !!$                         enddo
! !!$                         write(*,*) 'enthalpy change ', hbef, haft, abs((haft - hbef)/hbef)

!                          m_cond = m_cond - m_evap
!                          !write(*,*) 'f = ', f
!                       else
!                          write(*,*) 'max newton iter', abs(delT), precision, n, n_iter_newt
!                       endif
                                            
!                    endif
                   
                else
                   write(*,*) 'max newton'
                endif
                
             endif
             
          enddo
       endif
       
    enddo

! Update momentum! Will this cause error if layer is both rained out of and re-evaporated into? Don't think it's possible
    do k=1,npz
       du(k) = (ua(k) + umom(k)/delp(k))/(1. + delqv(k) + delql(k) + delqi(k)) - ua(k)
       dv(k) = (va(k) + vmom(k)/delp(k))/(1. + delqv(k) + delql(k) + delqi(k)) - va(k)
    enddo
    umom_aft = 0.0
    do k=1,npz
       umom_aft = umom_aft + (ua(k) + du(k))*delp(k)*(1 + delqv(k)  + delqi(k) + delql(k))
    enddo
    
    ! ========================================================
    ! Uncomment for enthalpy diagnostics
!    if (is_master() ) then
!    haft = 0.0
!    h_also = 0.0
!    maft_v = 0.0
!    do m=1,npz
!       call calc_enthalpy_mass(T_tmp(m), m_ds(m), m_vs(m), m_ls(m), m_is(m), htmp)
!       h_also  = h_also + htmp
!       call calc_enthalpy(T_tmp(m), qv_tmp(m)/(1. + delqv(m) + delql(m) + delqi(m)),&
!            ql(l)/(1. + delqv(m) + delql(m) + delqi(m)),&
!            qi(l)/(1. + delqv(m) + delql(m) + delqi(m)),&
!            htmp)
!       haft = haft+htmp*delp(m)*(1 + delqv(m) + delql(m) + delqi(m))
!       !write(*,*) 'm', m, qv_tmp(m)/(1+ delqv(m) + delql(m) + delqi(m)), m_vs(m)/(m_ds(m) + m_ls(m) + m_is(m) + m_vs(m)), qv_tmp(m)
!       maft_v = maft_v + delp(m)*qv_tmp(m)
!    enddo
!    endif
    ! ========================================================


  end subroutine rain_out_revap
  
  subroutine dk_adj(T2_guess, p1, p2, dp1, dp2, q1_0, q2_0, dk_a)
    real, intent(in) :: T2_guess, p1, p2, dp1, dp2, q1_0, q2_0
    real, intent(out) :: dk_a

    real :: qc_1, qc_2, qv_1, qv_2, qi_1, qi_2, k_out_up, k_out_dn, T_temp

    real :: eps = 1.e-8
    ! Estimate derivative of dk_adj wrt T2_guess

    call k_adj(T2_guess + eps/2, p1,p2,dp1,dp2,q1_0, q2_0, &
               qv_1,qv_2,qc_1,qc_2,qi_1, qi_2, T_temp,k_out_up)
    call k_adj(T2_guess - eps/2, p1,p2,dp1,dp2,q1_0, q2_0, &
               qv_1,qv_2,qc_1,qc_2,qi_1, qi_2, T_temp,k_out_dn)

    dk_a = (k_out_up - k_out_dn)/eps
  end subroutine dk_adj


     subroutine cold_trap(p, T, q, delp, ktrop)
       real, intent(inout) :: T(:)
       real, intent(inout) :: q(:)
       real, intent(inout) :: p(:), delp(:)
       integer, intent(in) :: ktrop

       integer ::k,m
       real,dimension(size(T)) :: rs,qs, r0s, qs2
       call r_sat(p,T,rs)
       r0s = q/(1-q)
       qs= rs/(1 + rs)
       qs2 = rs*(1-q)

       !write(*,*) 'delp', delp
       do k=ktrop,2,-1
          if( r0s(k) .lt. rs(k)) then
             !delp(k-1) = delp(k-1) + delp(k-1)*(qs2(k-1) - q(k-1))
             delp(k-1) = delp(k-1) + delp(k-1)*(rs(k-1) - r0s(k-1))/(1+r0s(k-1))
             q(k-1) = qs(k-1)
          else
             do m=1,k-1
                !delp(m) = delp(m) + delp(m)*(qs2(k) - q(m))
                delp(m) = delp(m) + delp(m)*(rs(k) - r0s(m))/(1 + r0s(m))
             enddo
             q(1:k-1) = qs(k)
             exit
          endif
       enddo

     end subroutine cold_trap
     
  subroutine gradient(p,T, dlnTdlnp, dlnpsat_dlnt)
    !==========================================================================
    ! Description
    !==========================================================================
    ! Calculates the moist adiabatic gradient d(log(T))/d(log(p)) as in
    ! Ding and Pierrehumbert 2016
    
    !==========================================================================
    ! Input variables
    !==========================================================================
    real, intent(in) :: p, T !Pressure and temperature

    !==========================================================================
    ! Output variables
    !==========================================================================
    real, intent(out) :: dlnTdlnp ! Moist adiabatic gradient

    real,intent(inout),optional :: dlnpsat_dlnt
    !==========================================================================
    ! Local variables
    !==========================================================================

    real :: L, qsat, rsat, num, denom, R_av, cp_av, varpi, f_ice, t2

    !==========================================================================
    ! Main body
    !==========================================================================

    call get_f_ice(T, f_ice)

    L = L_sub*f_ice + (hv_tp - hl_tp + (cp_v-cp_liquid)*(T-T_tp))*(1-f_ice)
    !L = hv_tp - hl_tp

    call r_sat_single(p,T,rsat)
    qsat = rsat/(1 + rsat)

    R_av  = (1 - qsat)*Rstar/mu_d + qsat*Rstar/mu_v
    cp_av = (1 - qsat)*cp_d + qsat*cp_v
    varpi = 1 - mu_d/mu_v

    ! Changed for more accurate version
    t2 = L/rvgas/T

    if (present(dlnpsat_dlnt)) then
       dlnpsat_dlnt = t2
    endif
    
    
    num = 1 - qsat + qsat*L/R_av/T * (1 - varpi*qsat)
    denom = 1- qsat + qsat*(L/cp_av/T)*(L*mu_v/Rstar/T)*(1 - varpi*qsat)

    dlnTdlnp = (num/denom) * R_av/cp_av
   !==========================================================================
    ! Sanity check of above expression
    !==========================================================================
    ! When qsat=0, tends to R_d/cp_d -- correct
    ! When qsat<<1, tends to (R_d/cp_d)* (1 + q*(L/R_d/T))/(1 + qsat*(L^2/cp_d/R_v/T^2))
    ! this is equal to expression in Pierrehumbert 2010 textbook + Oxford C5 course notes
    ! When qsat=1, only RH terms in num and denom remain to give R_v*T/L = pure steam
    ! adiabat -- correct
    !============================================================================
  end subroutine gradient
end module ding_convection
