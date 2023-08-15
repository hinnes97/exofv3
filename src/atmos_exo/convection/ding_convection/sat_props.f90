module sat_props_mod
  use phys, only : T_TP => H2O_TriplePointT, P_TP => H2O_TriplePointP, L_sub => H2O_L_sublimation, &
       cp_v => H2O_cp, mu_d => H2He_solar_MolecularWeight, &
       mu_v => H2O_MolecularWeight

  use constants_mod, only: rvgas
  implicit none

  ! Some physical constants for water
  real, parameter :: cp_liquid = 4219.9 ! Specific heat capacity of liquid water at TP [J/kg]
  real, parameter :: hv_tp = 2500920.   ! Specific enthalpy of water vapour at triple point
  real, parameter :: hl_tp = 1.0000000  ! Specific enthalpy of liquid water at the triple point

  real, parameter :: eps = mu_v/mu_d
  real, parameter :: T_low = T_tp - 20.0
contains

  subroutine get_xh2o(q_h2o, q_cond, x_h2o)
    real, intent(in) :: q_h2o
    real, intent(in) :: q_cond
    real, intent(out) :: x_h2o

    x_h2o = q_h2o/(eps*(1 - q_cond) - (eps-1)*q_h2o)
    
  end subroutine get_xh2o
  
    subroutine sat_vp(T, sat_p)
    real, intent(in) :: T
    real, intent(out) :: sat_p

    real :: L0, log_p_p0, f_ice, sat_p_l, sat_p_i

    call get_f_ice(T, f_ice)
    
    ! Over vapour
    L0 = hv_tp - hl_tp
    log_p_p0 = (1./T- 1./T_TP)*((cp_v - cp_liquid) - L0)/rvgas + (cp_v-cp_liquid)/rvgas*log(T/T_TP)
    sat_p_l = exp(log_p_p0)*P_TP

    ! Over ice
    L0 = L_sub
    sat_p_i = p_tp*exp(-L0/rvgas * (1./T - 1./T_tp) )

    
    sat_p = sat_p_i*f_ice + sat_p_l*(1-f_ice)

  end subroutine sat_vp

    subroutine r_sat_single(p, T, r)
    real, intent(in) :: p, T
    real, intent(out) :: r

    integer :: k
    real :: psat

    call sat_vp(T, psat)

    ! Make sure p - psat is not close to 0, limit minimum value
    r = eps*psat/(p - psat) 
    !q= eps*psat/p/(1 + (eps - 1)*psat/p) * (1-qc)
       
  end subroutine r_sat_single


    subroutine r_sat(p,T,r)
    real, intent(in) :: p(:), T(:)
    real, intent(out) :: r(:)

    integer :: k,npz
    real :: psat
    npz = size(p)

    do k=1,npz
       call r_sat_single(p(k),T(k),r(k))
    enddo
  end subroutine r_sat


  subroutine q_sat_single(p, T, qc, q)
    real, intent(in) :: p, T, qc
    real, intent(out) :: q

    real :: psat
    
    call sat_vp(T, psat)

    q = (1 - qc) * eps* (psat/p)/(1 + (eps-1)*(psat/p))
  end subroutine q_sat_single

  subroutine mv_mvmd(p, T, val)
    real, intent(in) :: p, T
    real, intent(out) :: val

    real :: psat

    call sat_vp(T, psat)

    val = eps*(psat/p)/(1 + (eps-1)*(psat/p))
    
  end subroutine mv_mvmd
  
!!$  subroutine q_sat_single(p, T, qt, q)
!!$    real, intent(in) :: p, T, qt
!!$    real, intent(out) :: q
!!$
!!$    integer :: k
!!$    real ::rsat
!!$
!!$    call r_sat_single(p,T,rsat)
!!$    q = rsat*(1-qt)
!!$    !q= eps*psat/p/(1 + (eps - 1)*psat/p) * (1-qc)
!!$       
!!$  end subroutine q_sat_single


  subroutine q_sat(p, T, qt, q)
    real, intent(in) :: p(:), T(:), qt(:)
    real, intent(out) :: q(:)

    integer :: k
    real :: psat

    do k=1,size(q)

       call q_sat_single(p(k), T(k), qt(k), q(k))
       !psat = p_sat(T(k))
    enddo
  end subroutine q_sat

  subroutine get_f_ice(T, f_ice)
    real, intent(in) :: T
    real, intent(out) :: f_ice

    if (T .lt. T_low) then
       f_ice = 1.0
    else if (T .gt. T_low .and. T .lt. T_tp) then
       f_ice = (T_tp - T)/(T_tp - T_low)
    else
       f_ice = 0.0
    endif
  end subroutine get_f_ice

end module sat_props_mod
