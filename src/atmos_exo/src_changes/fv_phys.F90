module fv_phys_mod

use time_manager_mod,   only: time_type
use fv_update_phys_mod, only: fv_update_phys
use fv_timing_mod,      only: timing_on, timing_off
use  exo_phys_mod,     only:   Exo_Tend

#ifdef MARS_GCM
use  Mars_phys_mod,     only:   Mars_phys
! MDH commented out 20/10
!#else
!use hswf_mod,           only: Held_Suarez_Strat, Held_Suarez_Tend, Sim_phys
#endif

implicit none
  logical:: nudge_initialized = .false.
  real, allocatable:: u0(:,:,:), v0(:,:,:), t0(:,:,:), dp(:,:,:)

public :: fv_phys, fv_nudge

contains

!-----------------------------------------------------------------------

 subroutine fv_phys(npx, npy, npz, is, ie, js, je, ng, nq,       &
                    u, v, w, pt, q, pe, delp, peln, pkz, pdt,    &
                    ua, va, phis, grid, ak, bk, ks, ps, pk, pk_nd,    &
                    u_srf, v_srf, ts, delz, hydrostatic, phys_hydrostatic, &
                    oro, strat, rayf, p_ref, fv_sg_adj, master,  &
                    do_Held_Suarez, Time, time_total, dt_atmos)


    integer, INTENT(IN   ) :: npx, npy, npz
    integer, INTENT(IN   ) :: is, ie, js, je, ng, nq
    integer, INTENT(IN   ) :: fv_sg_adj
    real, INTENT(IN) :: p_ref, dt_atmos
    real, INTENT(IN) :: oro(is:ie,js:je)

    real   , INTENT(INOUT) ::    u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
    real   , INTENT(INOUT) ::    v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) ::    w(is-ng:ie+  ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) ::   pt(is-ng:ie+  ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) ::   ts(is-ng:ie+  ng,js-ng:je+  ng)
    real   , INTENT(INOUT) :: delp(is-ng:ie+  ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) ::    q(is-ng:ie+  ng,js-ng:je+  ng,npz, nq)
    real   , INTENT(INOUT) ::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
    real   , INTENT(INOUT) :: peln(is   :ie     ,1:npz+1,js   :je     )
    real   , INTENT(INOUT) ::  pkz(is   :ie     ,js   :je     ,1:npz)
    real   , INTENT(INOUT) ::  pk (is   :ie     ,js   :je     ,npz+1)
    real,    INTENT(INOUT) :: pk_nd(is  :ie,     js:   je      ,npz+1)
    real   , INTENT(INOUT) ::   ua(is-ng:ie+ng,js-ng:je+ng,npz)
    real   , INTENT(INOUT) ::   va(is-ng:ie+ng,js-ng:je+ng,npz)
    real   , INTENT(INOUT) ::   ps(is-ng:ie+ng,js-ng:je+ng)

    ! Added by MDH 26/10/21
    !real :: t_dt_conv(is:ie,js:je,1:npz)

    real   , INTENT(IN   ) :: phis(is-ng:ie+ng,js-ng:je+ng)
    real   , INTENT(IN   ) :: grid(is-ng:ie+ng,js-ng:je+ng, 1:2)
    real   , INTENT(IN   ) :: ak(npz+1), bk(npz+1)
    integer, INTENT(IN   ) :: ks

    real   , INTENT(IN   )    :: pdt
    logical, INTENT(IN   )    :: strat, rayf, master, do_Held_Suarez
    real, INTENT(inout):: u_srf(is:ie,js:je)
    real, INTENT(inout):: v_srf(is:ie,js:je)


    type (time_type), intent(in) :: Time
    real, INTENT(IN), optional:: time_total
    logical, intent(in) ::  hydrostatic, phys_hydrostatic
    real, intent(inout) ::  delz(is:ie,js:je,npz)
    real, allocatable:: u_dt(:,:,:), v_dt(:,:,:), t_dt(:,:,:), q_dt(:,:,:,:), ts_dt(:,:)
    logical moist_phys

    integer  isd, ied, jsd, jed

#ifdef MARS_GCM
    real, allocatable:: qratio(:,:,:)
    integer :: i, j, k, m
#endif MARS_GCM


    isd = is-ng;   ied = ie + ng
    jsd = js-ng;   jed = je + ng

       allocate ( u_dt(isd:ied,jsd:jed,npz) )
       allocate ( v_dt(isd:ied,jsd:jed,npz) )
       allocate ( t_dt(is:ie,js:je,npz) )
       allocate ( ts_dt(is:ie,js:je) )

       allocate ( q_dt(is:ie,js:je,npz,nq) )
       u_dt = 0.
       v_dt = 0.
       t_dt = 0.
       q_dt = 0.

#ifdef MARS_GCM
       allocate ( qratio(is:ie,js:je,npz) )

       call timing_on('Mars_PHYS')
       call Mars_phys(npx, npy, npz, is, ie, js, je, ng, nq,   &
                     u_dt, v_dt, t_dt, q_dt, ua, va, pt, q,   &
                     phis, pe, delp, peln, pdt, grid, ak, bk,       &
                     qratio, rayf, master, Time, time_total )

       call timing_off('Mars_PHYS')

!!!      Be sure to set moist_phys to true
       call timing_on('UPDATE_PHYS')

       call fv_update_phys (pdt, is, ie, js, je, isd, ied, jsd, jed, ng, nq,   &
                            u, v, delp, pt, q, ua, va, ps, pe, peln, pk, pk_nd,pkz,  &
                            ak, bk, phis, u_srf, v_srf, ts, delz, hydrostatic, &
                            u_dt, v_dt, t_dt, q_dt, .true., Time )


!-----------------------------------------
! Adjust mass mixing ratio of all tracers:
!!!! May NEED TO FIX THIS?!?!  Perhaps do this for prognostic tracers only, and
!!!!      let the physics code figure out how to handle diagnostic tracers
!-----------------------------------------
     do m=1,nq
       do k=1,npz
         do j=js,je
             do i= is,ie
                q(i,j,k,m) = q(i,j,k,m) / qratio(i,j,k)
             enddo
          enddo
       enddo
     enddo

       deallocate ( u_dt )
       deallocate ( v_dt )
       deallocate ( t_dt )
       deallocate ( q_dt )
       deallocate ( qratio )
       call timing_off('UPDATE_PHYS')

! ------------------------------------------------
#else

! MDH commented out 20/10
!    if( do_Held_Suarez ) then
!       moist_phys = .false.
!      call Held_Suarez_Strat(npx, npy, npz, is, ie, js, je, ng, nq,  &
!                             u, v, pt, q, pe, delp, peln, pkz, pdt,  &
!                             ua, va, grid, ak, bk, ks, strat,  &
!!                             rayf, master, Time, time_total)
!       call Held_Suarez_Tend(npx, npy, npz, is, ie, js, je, ng, nq,   &
!                              u, v, pt, q, pe, delp, peln, pkz, pdt,  &
!                              ua, va, u_dt, v_dt, t_dt, q_dt, grid,   &
!                              delz, hydrostatic, ak, bk, ks,          &
!                              strat, rayf, master, Time, time_total)
!    else
!       moist_phys = .true.
!                                             call timing_on('SIM_PHYS')
!       call Sim_phys(npx, npy, npz, is, ie, js, je, ng, nq,                &
!                     u_dt, v_dt, t_dt, q_dt, u, v, w, ua, va, pt, delz, q, &
!                     pe, delp, peln, oro, hydrostatic, phys_hydrostatic,   &
!                     pdt, grid, ak, bk, p_ref, fv_sg_adj, master, Time, time_total)
!                                            call timing_off('SIM_PHYS')
!    endif
       write(*,*) 'BEFORE EXOTEND'
       call Exo_Tend(npx, npy, npz, is, ie, js, je, ng, nq,   &
                              u, v, w, pt, q, ts, pe, delp, peln, pkz, pdt,  &
                              ua, va, u_dt, v_dt, t_dt, ts_dt, q_dt, grid,   &
                              delz, hydrostatic, ak, bk, ks,          &
                              strat, rayf, master, Time, dt_atmos)

       call timing_on('UPDATE_PHYS')
    call fv_update_phys (pdt, is, ie, js, je, isd, ied, jsd, jed, ng, nq,   &
                         u, v, delp, pt, q, ua, va, ps, pe, peln, pk, pk_nd,pkz,  &
                         ak, bk, phis, u_srf, v_srf, ts, delz, hydrostatic, &
                         u_dt, v_dt, t_dt, ts_dt, q_dt, moist_phys, Time )

    deallocate ( u_dt )
    deallocate ( v_dt )
    deallocate ( t_dt )
    deallocate ( q_dt )

       call timing_off('UPDATE_PHYS')

#endif

 end subroutine fv_phys


 subroutine fv_nudge( npz, is, ie, js, je, ng, u, v, delp, pt, dt,    &
                      tau_winds, tau_press, tau_temp )
!
! Nudge the prognostic varaibles toward the IC
! This is only useful for generating balanced steady state IC

    real   , INTENT(IN   ) :: tau_winds, tau_press, tau_temp, dt
    integer, INTENT(IN   ) :: npz, is, ie, js, je, ng
    real   , INTENT(INOUT) ::    u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
    real   , INTENT(INOUT) ::    v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) :: delp(is-ng:ie+  ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) ::   pt(is-ng:ie+  ng,js-ng:je+  ng,npz)
    real c_w, c_p, c_t
    integer  i, j, k

    if ( .not. nudge_initialized ) then

         if( tau_winds > 0. ) then
             allocate ( u0(is:ie,  js:je+1,npz) )
             allocate ( v0(is:ie+1,js:je  ,npz) )
             do k=1,npz
                do j=js,je+1
                   do i=is,ie
                      u0(i,j,k) = u(i,j,k) 
                   enddo
                enddo
                do j=js,je
                   do i=is,ie+1
                      v0(i,j,k) = v(i,j,k)
                   enddo
                enddo
             enddo
         endif

         if( tau_press > 0. ) then
             allocate ( dp(is:ie,js:je,npz) )
             do k=1,npz
                do j=js,je
                   do i=is,ie
                      dp(i,j,k) = delp(i,j,k)
                   enddo
                enddo
             enddo
         endif

         if( tau_temp > 0. ) then
             allocate ( t0(is:ie,js:je,npz) )
             do k=1,npz
                do j=js,je
                   do i=is,ie
                      t0(i,j,k) = pt(i,j,k)
                   enddo
                enddo
             enddo
         endif

         nudge_initialized = .true.
         return
    endif

! Nudge winds to initial condition:

      do k=1,npz
         if( tau_winds > 0. ) then
            c_w = dt/tau_winds
            do j=js,je+1
               do i=is,ie
                  u(i,j,k) = (u(i,j,k)+c_w*u0(i,j,k)) / (1.+c_w)
               enddo
            enddo
            do j=js,je
               do i=is,ie+1
                  v(i,j,k) = (v(i,j,k)+c_w*v0(i,j,k)) / (1.+c_w)
               enddo
            enddo
         endif

         if ( tau_press > 0. ) then
             c_p = dt/tau_press
             do j=js,je
                do i=is,ie
                   delp(i,j,k) = (delp(i,j,k)+c_p*dp(i,j,k)) / (1.+c_p)
                enddo
             enddo
         endif

         if ( tau_temp > 0. ) then
             c_t = dt/tau_temp
             do j=js,je
                do i=is,ie
                   pt(i,j,k) = (pt(i,j,k)+c_t*t0(i,j,k)) / (1.+c_t)
                enddo
             enddo
         endif
      enddo

 end subroutine fv_nudge


end module fv_phys_mod
