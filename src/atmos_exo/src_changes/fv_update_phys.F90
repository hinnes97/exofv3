module fv_update_phys_mod

  use constants_mod,      only: kappa, rdgas, rvgas, grav, cp_air, pi, cp_air, cp_vapor
  use field_manager_mod,  only: MODEL_ATMOS
  use mpp_domains_mod,    only: mpp_update_domains
  use mpp_parameter_mod,  only: AGRID_PARAM=>AGRID
  use mpp_mod,            only: FATAL, mpp_error
  use time_manager_mod,   only: time_type
  use tracer_manager_mod, only: get_tracer_index

  use fv_arrays_mod,      only: fv_atmos_type
  use fv_control_mod,     only: npx, npy, npz, ncnst, k_top, nwat, fv_debug, &
                                tau_h2o, phys_hydrostatic, dwind_2d, filter_phys,&
                                non_dilute, master
  use fv_mp_mod,          only: domain, gid
  use fv_eta_mod,         only: get_eta_level
  use fv_grid_utils_mod,  only: edge_vect_s,edge_vect_n,edge_vect_w,edge_vect_e, &
                                es, ew, vlon, vlat, z11, z12, z21, z22, &
                                sina_u, sina_v, da_min
  use fv_grid_tools_mod,  only: dx, dy, rdxc, rdyc, rarea, dxa, dya, grid_type
  use fv_timing_mod,      only: timing_on, timing_off
  use fv_diagnostics_mod, only: prt_maxmin
  use interp_mod,         only: get_pk_edge
#ifdef CLIMATE_NUDGE
  use atmos_nudge_mod,    only: get_atmos_nudge, do_ps
#else
  use fv_nwp_nudge_mod,   only: fv_nwp_nudge
#endif

  implicit none

  public :: fv_update_phys, del2_phys

  contains

  subroutine fv_update_phys ( dt, is, ie, js, je, isd, ied, jsd, jed, ng, nq,     &
                              u, v, delp, pt, q, ua, va, ps, pe,  peln, pk, pk_nd,pkz,  &
                              ak, bk, phis, u_srf, v_srf, ts, delz, hydrostatic,  &
                              u_dt, v_dt, t_dt, ts_dt, q_dt, moist_phys, Time, nudge )
    real, intent(in)   :: dt
    integer, intent(in):: is,  ie,  js,  je, ng
    integer, intent(in):: isd, ied, jsd, jed
    integer, intent(in):: nq            ! tracers modified by physics 
                                        ! ncnst is the total nmber of tracers
    logical, intent(in):: moist_phys
    logical, intent(in):: hydrostatic
    logical, optional, intent(in):: nudge

    type (time_type), intent(in) :: Time

    real, intent(in), dimension(npz+1):: ak, bk
    real, intent(in) :: phis(isd:ied,jsd:jed)
    real, intent(inout):: delz(is:ie,js:je,npz)

! Winds on lat-lon grid:
    real, intent(inout), dimension(isd:ied,jsd:jed,npz):: ua, va

! Tendencies from Physics:
    real, intent(inout), dimension(isd:ied,jsd:jed,npz):: u_dt, v_dt
    real, intent(inout):: t_dt(is:ie,js:je,npz)
    real, intent(inout):: ts_dt(is:ie,js:je)
    real, intent(inout):: q_dt(is:ie,js:je,npz,nq)

! Saved Bottom winds for GFDL Physics Interface
    real, intent(out), dimension(is:ie,js:je):: u_srf, v_srf

    real, intent(inout), dimension(isd:ied,jsd:jed):: ts

    real, intent(inout):: u(isd:ied  ,jsd:jed+1,npz)  ! D grid zonal wind (m/s)
    real, intent(inout):: v(isd:ied+1,jsd:jed  ,npz)  ! D grid meridional wind (m/s)
    real, intent(inout), dimension(isd:ied,jsd:jed,npz):: pt, delp
    real, intent(inout):: q(isd:ied,jsd:jed,npz, ncnst)   ! specific humidity and constituents

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, intent(inout):: ps  (isd:ied  ,jsd:jed)           ! Surface pressure (pascal)
    real, intent(inout):: pe  (is-1:ie+1, npz+1,js-1:je+1)  ! edge pressure (pascal)
    real, intent(inout):: pk  (is:ie,js:je  , npz+1)        ! pe**cappa
    real, intent(inout):: peln(is:ie,npz+1,js:je)           ! ln(pe)
    real, intent(inout):: pkz (is:ie,js:je,npz)             ! finite-volume mean pk

    real, intent(inout) :: pk_nd(is:ie,js:je,npz+1)

!**********
! Halo Data
!**********
    real, parameter::    q1_h2o = 2.2E-6
    real, parameter::    q7_h2o = 3.8E-6
    real, parameter::  q100_h2o = 3.8E-6
    real, parameter:: q1000_h2o = 3.1E-6
    real, parameter:: q2000_h2o = 2.8E-6
    real, parameter:: q3000_h2o = 3.0E-6

! Local arrays:
    real  ps_dt(is:ie,js:je)
    real  phalf(npz+1), pfull(npz)
    real kap_loc(is:ie,js:je,1:npz)
    real kap_edge(is:ie,js:je,1:npz+1)
    real pkz_test(is:ie,js:je,1:npz)
    integer  i, j, k, m
    integer  sphum, liq_wat, ice_wat, cld_amt   ! GFDL AM physics
    integer vap
    integer  rainwat, snowwat, graupel          ! Lin Micro-physics
    real     qstar, dbk, rdg, gama_dt, zvir
    integer ind(3)
    if ( filter_phys ) then
         call del2_phys(t_dt, delp, 0.2, npx, npy, npz, is, ie, js, je, &
                        isd, ied, jsd, jed, 0)
         do m=1,nq
            call del2_phys(q_dt(:,:,:,m), delp, 0.2, npx, npy, npz, is, ie, js, je, &
                           isd, ied, jsd, jed, 0)
         enddo
    endif

    rdg = -rdgas / grav
    gama_dt = dt * cp_air / (cp_air-rdgas)

#if defined(MARS_GCM) || defined(VENUS_GCM)
!$omp parallel do default(shared) private(i, j, k, m, qstar, ps_dt)
    do k=1, npz
       do j=js,je
          do i=is,ie
             ua(i,j,k) = ua(i,j,k) + dt*u_dt(i,j,k)
             va(i,j,k) = va(i,j,k) + dt*v_dt(i,j,k)
          enddo
       enddo

       if ( hydrostatic ) then
          do j=js,je
             do i=is,ie
                pt(i,j,k) = pt(i,j,k) + dt*t_dt(i,j,k)
             enddo
          enddo
       else
         if ( phys_hydrostatic ) then
! Heating/cooling from physics is assumed to be isobaric hydrostatic proc
! "nagative" definiteness of delz is maintained.
             do j=js,je
                do i=is,ie
                   delz(i,j,k) = delz(i,j,k) / pt(i,j,k)
               pt(i,j,k) = pt(i,j,k) + dt*t_dt(i,j,k)
                   delz(i,j,k) = delz(i,j,k) * pt(i,j,k)
                enddo
             enddo
         else
! Convert tendency from constant-p to constant-volume
             do j=js,je
                do i=is,ie
                   pt(i,j,k) = pt(i,j,k) + t_dt(i,j,k)*gama_dt
                enddo
             enddo
         endif
       endif

!----------------
! Update tracers:
!----------------
       do m=1,nq
          do j=js,je
             do i=is,ie
                q(i,j,k,m) = q(i,j,k,m) + dt*q_dt(i,j,k,m)
             enddo
          enddo
       enddo
    enddo  ! openmp k-loop

#else
    sphum   = 1

    if ( moist_phys ) then
        cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')
           zvir = rvgas/rdgas - 1.
    else
        cld_amt = 7
           zvir = 0.
    endif

    if ( nwat>=3 ) then
        sphum   = get_tracer_index (MODEL_ATMOS, 'sphum')
        liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
        ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
    endif

    if ( nwat==6 ) then
! Micro-physics:
        rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
        snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
        graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
        if ( cld_amt<7 ) call mpp_error(FATAL,'Cloud Fraction allocation error') 
    endif

    if ( fv_debug ) then
       if ( gid==0 ) write(*,*) nq, nwat, sphum, liq_wat, ice_wat, rainwat, snowwat, graupel
       call prt_maxmin('delp_b_update', delp, is, ie, js,  je, ng, npz, 0.01, gid==0)
       do m=1,nq
          call prt_maxmin('q_dt', q_dt(is,js,1,m), is, ie, js, je, 0, npz, 1., gid==0)
       enddo
    endif

    call get_eta_level(npz, 1.0E5, pfull, phalf, ak, bk)

!$omp parallel do default(shared) private(i, j, k, m, qstar, ps_dt)
    do k=1, npz

! Do idealized Ch4 chemistry
       if ( tau_h2o>0.0 .and. pfull(k) < 3000. ) then

           if ( pfull(k) < 1. ) then
               qstar = q1_h2o
           elseif ( pfull(k) <   7. .and. pfull(k) >=    1. ) then
               qstar = q1_h2o + (q7_h2o-q1_h2o)*log(pfull(k)/1.)/log(7.)
           elseif ( pfull(k) <  100. .and. pfull(k) >=    7. ) then
               qstar = q7_h2o + (q100_h2o-q7_h2o)*log(pfull(k)/7.)/log(100./7.)
           elseif ( pfull(k) < 1000. .and. pfull(k) >=  100. ) then
               qstar = q100_h2o + (q1000_h2o-q100_h2o)*log(pfull(k)/1.E2)/log(10.)
           elseif ( pfull(k) < 2000. .and. pfull(k) >= 1000. ) then
               qstar = q1000_h2o + (q2000_h2o-q1000_h2o)*log(pfull(k)/1.E3)/log(2.)
           else
               qstar = q3000_h2o
           endif

           do j=js,je
              do i=is,ie
                 q_dt(i,j,k,sphum) = q_dt(i,j,k,sphum) + (qstar-q(i,j,k,sphum))/(tau_h2o*86400.)
              enddo
           enddo
       endif

       do j=js,je
          do i=is,ie
             ua(i,j,k) = ua(i,j,k) + dt*u_dt(i,j,k)
             va(i,j,k) = va(i,j,k) + dt*v_dt(i,j,k)
          enddo
       enddo

       if ( hydrostatic ) then
          do j=js,je
             do i=is,ie
!               pt(i,j,k) = pt(i,j,k) + dt*t_dt(i,j,k) /     &
!                           (1.-kappa*ak(1)/delp(i,j,k)*(peln(i,k+1,j)-peln(i,k,j)))
                pt(i,j,k) = pt(i,j,k) + dt*t_dt(i,j,k)
                if (k==1) then
                   ts(i,j) = ts(i,j) + dt*ts_dt(i,j)
                end if
             enddo
          enddo
       else
         if ( phys_hydrostatic ) then
! Heating/cooling from physics is assumed to be isobaric hydrostatic proc
! "nagative" definiteness of delz is maintained.
             do j=js,je
                do i=is,ie
                   delz(i,j,k) = delz(i,j,k) / pt(i,j,k)
!                     pt(i,j,k) = pt(i,j,k) + dt*t_dt(i,j,k) /     &
!                              (1.-kappa*ak(1)/delp(i,j,k)*(peln(i,k+1,j)-peln(i,k,j)))
               pt(i,j,k) = pt(i,j,k) + dt*t_dt(i,j,k)
                   delz(i,j,k) = delz(i,j,k) * pt(i,j,k)
                enddo
             enddo
         else
! Convert tendency from constant-p to constant-volume
             do j=js,je
                do i=is,ie
                   pt(i,j,k) = pt(i,j,k) + t_dt(i,j,k)*gama_dt
                enddo
             enddo
         endif
       endif

!----------------
! Update tracers:
!----------------
       do m=1,nq
          do j=js,je
             do i=is,ie
                q(i,j,k,m) = q(i,j,k,m) + dt*q_dt(i,j,k,m)
             enddo
          enddo
       enddo

!--------------------------------------------------------
! Adjust total air mass due to changes in water substance
!--------------------------------------------------------

      if ( nwat==6 ) then
! micro-physics with 6 water substances
        do j=js,je
           do i=is,ie
              ps_dt(i,j)  = 1. + dt * ( q_dt(i,j,k,sphum  ) +    &
                                        q_dt(i,j,k,liq_wat) +    &
                                        q_dt(i,j,k,rainwat) +    &
                                        q_dt(i,j,k,ice_wat) +    &
                                        q_dt(i,j,k,snowwat) +    &
                                        q_dt(i,j,k,graupel) )
              delp(i,j,k) = delp(i,j,k) * ps_dt(i,j)
           enddo
        enddo
      elseif( nwat==3 ) then
! GFDL AM2/3 phys (cloud water + cloud ice)
        do j=js,je
           do i=is,ie
               ps_dt(i,j) = 1. + dt*(q_dt(i,j,k,sphum  ) +    &
                                     q_dt(i,j,k,liq_wat) +    &
                                     q_dt(i,j,k,ice_wat) )
              delp(i,j,k) = delp(i,j,k) * ps_dt(i,j)
           enddo
        enddo
      elseif ( nwat>0 ) then
        do j=js,je
           do i=is,ie
              ps_dt(i,j)  = 1. + dt*sum(q_dt(i,j,k,1:nwat))
              delp(i,j,k) = delp(i,j,k) * ps_dt(i,j)
           enddo
        enddo
      endif

!-----------------------------------------
! Adjust mass mixing ratio of all tracers 
!-----------------------------------------
      if ( nwat /=0 ) then
        do m=1,ncnst   
          if( m /= cld_amt ) then  ! cloud fraction in GFDL physics
            do j=js,je
               do i=is,ie
                  q(i,j,k,m) = q(i,j,k,m) / ps_dt(i,j)
               enddo
            enddo
          endif
        enddo
      endif
   enddo ! openmp k-loop

#endif ! Mars and Venus GCM

! [delp, (ua, va), pt, q] updated. Perform nudging if requested

!------- nudging of atmospheric variables toward specified data --------

    ps_dt(:,:) = 0.

    if (present(nudge)) then
#ifdef CLIMATE_NUDGE
!--------------------------------------------
! All fields will be updated; tendencies added
!--------------------------------------------
      if (nudge) then
        call get_atmos_nudge ( Time, dt, beglon, endlon, beglat, endlat,    &
             npz, ng, ps(beglon:endlon,:), ua(beglon:endlon,:,:), &
             va(beglon:endlon,:,:), pt(beglon:endlon,:,:), &
             q(beglon:endlon,:,:,:), ps_dt(beglon:endlon,:), u_dt(beglon:endlon,:,:),  & 
             v_dt(beglon:endlon,:,:), t_dt(beglon:endlon,:,:), &
             q_dt(beglon:endlon,:,:,:) )

        if (do_ps) then
!--------------
! Update delp
!--------------
            do k=1,npz
               dbk = dt * (bk(k+1) - bk(k))
               do j=js,je
                  do i=is,ie
                     delp(i,j,k) = delp(i,j,k) + dbk*ps_dt(i,j)
                  enddo
               enddo
            enddo
        endif
     endif
#else
! All fields will be updated except winds; wind tendencies added
      if (nudge) then
!$omp parallel do default(shared) private(i, j, k)
        do j=js,je
         do k=2,npz+1                                                                             
          do i=is,ie
            pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
          enddo
         enddo
         do i=is,ie
           ps(i,j) = pe(i,npz+1,j)
         enddo
        enddo
        call fv_nwp_nudge ( Time, dt, npz,  ps_dt, u_dt, v_dt, t_dt, q_dt,   &
                            zvir, ak, bk, ts, ps, delp, ua, va, pt, nwat, q,  phis )
      endif
#endif
   endif

!----------------------------------------
! Update pe, peln, pkz, and surface winds
!----------------------------------------
  if ( fv_debug ) then
       call prt_maxmin('PS_b_update',     ps, is, ie, js,  je, ng,   1, 0.01, gid==0)
!       call prt_maxmin('pd_dt',  ps_dt, is, ie, js,  je, 0,   1, 1., gid==0)
       call prt_maxmin('delp_a_update', delp, is, ie, js,  je, ng, npz, 0.01, gid==0)
  endif

!$omp parallel do default(shared) private(i, j, k)
   do j=js,je
      do k=2,npz+1                                                                             
         do i=is,ie
              pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
            peln(i,k,j) = log( pe(i,k,j) )
            pk(i,j,k) = exp( kappa*peln(i,k,j) )
         enddo
      enddo

      do i=is,ie
            ps(i,j) = pe(i,npz+1,j)
         u_srf(i,j) = ua(i,j,npz)
         v_srf(i,j) = va(i,j,npz)
      enddo

      if ( hydrostatic) then! .and. .not. non_dilute) then! .and. .not. non_dilute) then
         do k=1,npz
            do i=is,ie
               pkz_test(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(kappa*(peln(i,k+1,j)-peln(i,k,j)))
            enddo
         enddo
      endif
   enddo      ! j-loop

! Update pk_nd and pkz
   if (non_dilute) then
      vap = get_tracer_index(MODEL_ATMOS, 'vapour')
      do k=1,npz
         do j=js,je
            do i=is,ie
               kap_loc(i,j,k) = (q(i,j,k,vap)*rvgas + (1.-q(i,j,k,vap))*rdgas)/&
                    (q(i,j,k,vap)*cp_vapor + (1.-q(i,j,k,vap))*cp_air)
            enddo
         enddo
      enddo

      do k=1,npz
         do j=js,je
            do i=is,ie
               pkz(i,j,k) = exp(kap_loc(i,j,k)*log(delp(i,j,k)/(peln(i,k+1,j) - peln(i,k,j))))
               !pkz(i,j,k) = (pk_nd(i,j,k+1) - pk_nd(i,j,k))/(kap_loc(i,j,k)*(peln(i,k+1,j)-peln(i,k,j)))
            enddo
         enddo
      enddo
!      if (gid==2) then
!         write(*,*) '=-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-'
!         write(*,*) 'PKZ COMP'
!         write(*,*) '=-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-'
!         ind =  maxloc(abs(pk_nd(is:ie,js:je,1:npz+1) - pk(is:ie,js:je,1:npz+1))/pk(is:ie,js:je,1:npz+1))
!         write(*,*) pk(is+ind(1)-1,js+ind(2)-1,ind(3)), pk_nd(is+ind(1)-1,js+ind(2)-1,ind(3)), kap_edge(is+ind(1)-1,js+ind(2)-1,ind(3))
!         write(*,*) ind(1)-1, ind(2)-1, ind(3)
!         write(*,*) is+ind(1)-1, js+ind(2)-1, ind(3)
!         write(*,*) ie, je, npz
!         write(*,*) peln(is+ind(1)-1,ind(3),js+ind(2)-1), exp(peln(is+ind(1)-1, ind(3), js+ind(2)-1)), pe(is+ind(1)-1,ind(3),js+ind(2)-1)
!         write(*,*) ps(is+ind(1)-1,js+ind(2)-1)
!         ind =  maxloc(abs(pkz_test(is:ie,js:je,1:npz) - pkz(is:ie,js:je,1:npz))/pkz(is:ie,js:je,1:npz))
!         write(*,*) pkz_test(is+ind(1)-1,js+ind(2)-1,ind(3)), pkz(is+ind(1)-1,js+ind(2)-1,ind(3))
!         write(*,*) ind(3)
!         !write(*,*) pk(ind(1),ind(2),ind(3)+1),pk_nd(ind(1),ind(2),ind(3)+1)
!         !write(*,*) kap_loc(ind(1),ind(2),ind(3)), kappa
!         !write(*,*) kap_edge(ind(1), ind(2), ind(3)),kappa
!         write(*,*) '=-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-'
!      endif
      
   endif

!   if (gid==2) then
!      write(*,*) 'hello world', is, ie, js,je 
!      write(*,*) peln(8,51,15), pe(8,51,15)
!   endif
   
!-------------------------------------------------------------------------
! Re-compute the full (nonhydrostatic) pressure due to temperature changes
!-------------------------------------------------------------------------
    if ( .not.hydrostatic ) then
      if ( k_top>1 ) then
         do k=1,k_top-1
            do j=js,je
               do i=is,ie
                  pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k)) /     &
                         (kappa*(peln(i,k+1,j)-peln(i,k,j)))
               enddo
            enddo
         enddo
      endif

      do k=k_top,npz
         do j=js,je
            do i=is,ie
! perfect gas law: p = density * rdgas * virtual_temperature
!              pkz(i,j,k) = ( rdg*delp(i,j,k)*pt(i,j,k)/delz(i,j,k) )**kappa
               pkz(i,j,k) = exp( kappa*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
                                           (1.+zvir*q(i,j,k,sphum))/delz(i,j,k)) )
            enddo
         enddo
      enddo
   endif
                                                    call timing_on(' Update_dwinds')
  if ( dwind_2d ) then
    call update2d_dwinds_phys(is, ie, js, je, isd, ied, jsd, jed, dt, u_dt, v_dt, u, v)
  else
    call update_dwinds_phys(is, ie, js, je, isd, ied, jsd, jed, dt, u_dt, v_dt, u, v)
  endif
                                                    call timing_off(' Update_dwinds')

  if ( fv_debug ) then
       call prt_maxmin('PS_a_update', ps, is, ie, js, je, ng,   1, 0.01, gid==0)
  endif

  end subroutine fv_update_phys


  subroutine del2_phys(qdt, delp, cd, npx, npy, km, is, ie, js, je, &
                       isd, ied, jsd, jed, ngc)
! This routine is for filtering the physics tendency
   integer, intent(in):: npx, npy, km
   integer, intent(in):: is, ie, js, je, isd, ied, jsd, jed, ngc
   real,    intent(in):: cd            ! cd = K * da_min;   0 < K < 0.25
   real, intent(in   ):: delp(isd:ied,jsd:jed,km)
   real, intent(inout):: qdt(is-ngc:ie+ngc,js-ngc:je+ngc,km)
!
   real :: q(isd:ied,jsd:jed,km)
   real :: fx(is:ie+1,js:je), fy(is:ie,js:je+1)
   real :: mask(is:ie+1,js:je+1)
   real :: f1(is:ie+1), f2(js:je+1)
   real :: damp
   integer i,j,k

! Applying mask to cd, the damping coefficient?
   damp = 0.25 * cd * da_min

! Mask defined at corners
   do i=is,ie+1
      f1(i) = (1. - sin(real(i-1)/real(npx-1)*pi))**2
   enddo

   do j=js,je+1
      f2(j) = (1. - sin(real(j-1)/real(npy-1)*pi))**2
      do i=is,ie+1
         mask(i,j) = damp * (f1(i) + f2(j))
      enddo
   enddo

! mass weighted tendency from physics is filtered
   do k=1,km
      do j=js,je
         do i=is,ie
            q(i,j,k) = qdt(i,j,k)*delp(i,j,k)
         enddo
      enddo
   enddo
                     call timing_on('COMM_TOTAL')
   call mpp_update_domains(q, domain, complete=.true.)
                     call timing_off('COMM_TOTAL')

   do k=1,km
      do j=js,je
         do i=is,ie+1
            fx(i,j) = (mask(i,j)+mask(i,j+1))*dy(i,j)*sina_u(i,j)*(q(i-1,j,k)-q(i,j,k))*rdxc(i,j)
         enddo
      enddo
      do j=js,je+1
         do i=is,ie
            fy(i,j) = (mask(i,j)+mask(i+1,j))*dx(i,j)*sina_v(i,j)*(q(i,j-1,k)-q(i,j,k))*rdyc(i,j)
         enddo
      enddo
      do j=js,je
         do i=is,ie
            qdt(i,j,k) = qdt(i,j,k) + rarea(i,j)*(fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))/delp(i,j,k)
         enddo
      enddo
   enddo

  end subroutine del2_phys


  subroutine update_dwinds_phys(is, ie, js, je, isd, ied, jsd, jed, dt, u_dt, v_dt, u, v)

! Purpose; Transform wind tendencies on A grid to D grid for the final update
 
  integer, intent(in):: is,  ie,  js,  je
  integer, intent(in):: isd, ied, jsd, jed
  real,    intent(in):: dt
  real, intent(inout):: u(isd:ied,  jsd:jed+1,npz)
  real, intent(inout):: v(isd:ied+1,jsd:jed  ,npz)
  real, intent(inout), dimension(isd:ied,jsd:jed,npz):: u_dt, v_dt

! local:
  real v3(is-1:ie+1,js-1:je+1,3)
  real ue(is-1:ie+1,js:je+1,3)    ! 3D winds at edges
  real ve(is:ie+1,js-1:je+1,  3)    ! 3D winds at edges
  real, dimension(is:ie):: ut1, ut2, ut3
  real, dimension(js:je):: vt1, vt2, vt3
  real dt5, gratio
  integer i, j, k, m, im2, jm2


       call timing_on('COMM_TOTAL')
  call mpp_update_domains(u_dt, domain, complete=.false.)
  call mpp_update_domains(v_dt, domain, complete=.true.)
       call timing_off('COMM_TOTAL')

    dt5 = 0.5 * dt
    im2 = (npx-1)/2
    jm2 = (npy-1)/2

!$omp parallel do default(shared) private(i, j, k, ut1, ut2, ut3, vt1, vt2, vt3, ue, ve, v3)
    do k=1, npz

     if ( grid_type > 3 ) then    ! Local & one tile configurations

       do j=js,je+1
          do i=is,ie
             u(i,j,k) = u(i,j,k) + dt5*(u_dt(i,j-1,k) + u_dt(i,j,k))
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             v(i,j,k) = v(i,j,k) + dt5*(v_dt(i-1,j,k) + v_dt(i,j,k))
          enddo
       enddo

     else
! Compute 3D wind tendency on A grid
       do j=js-1,je+1
          do i=is-1,ie+1
             v3(i,j,1) = u_dt(i,j,k)*vlon(i,j,1) + v_dt(i,j,k)*vlat(i,j,1)
             v3(i,j,2) = u_dt(i,j,k)*vlon(i,j,2) + v_dt(i,j,k)*vlat(i,j,2)
             v3(i,j,3) = u_dt(i,j,k)*vlon(i,j,3) + v_dt(i,j,k)*vlat(i,j,3)
          enddo
       enddo

! Interpolate to cell edges
       do j=js,je+1
          do i=is-1,ie+1
             ue(i,j,1) = v3(i,j-1,1) + v3(i,j,1)
             ue(i,j,2) = v3(i,j-1,2) + v3(i,j,2)
             ue(i,j,3) = v3(i,j-1,3) + v3(i,j,3)
          enddo
       enddo

       do j=js-1,je+1
          do i=is,ie+1
             ve(i,j,1) = v3(i-1,j,1) + v3(i,j,1)
             ve(i,j,2) = v3(i-1,j,2) + v3(i,j,2)
             ve(i,j,3) = v3(i-1,j,3) + v3(i,j,3)
          enddo
       enddo

! --- E_W edges (for v-wind):
     if ( is==1 ) then
       i = 1
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = edge_vect_w(j)*ve(i,j-1,1)+(1.-edge_vect_w(j))*ve(i,j,1)
             vt2(j) = edge_vect_w(j)*ve(i,j-1,2)+(1.-edge_vect_w(j))*ve(i,j,2)
             vt3(j) = edge_vect_w(j)*ve(i,j-1,3)+(1.-edge_vect_w(j))*ve(i,j,3)
        else
             vt1(j) = edge_vect_w(j)*ve(i,j+1,1)+(1.-edge_vect_w(j))*ve(i,j,1)
             vt2(j) = edge_vect_w(j)*ve(i,j+1,2)+(1.-edge_vect_w(j))*ve(i,j,2)
             vt3(j) = edge_vect_w(j)*ve(i,j+1,3)+(1.-edge_vect_w(j))*ve(i,j,3)
        endif
       enddo
       do j=js,je
          ve(i,j,1) = vt1(j)
          ve(i,j,2) = vt2(j)
          ve(i,j,3) = vt3(j)
       enddo
     endif
     if ( (ie+1)==npx ) then
       i = npx
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = edge_vect_e(j)*ve(i,j-1,1)+(1.-edge_vect_e(j))*ve(i,j,1)
             vt2(j) = edge_vect_e(j)*ve(i,j-1,2)+(1.-edge_vect_e(j))*ve(i,j,2)
             vt3(j) = edge_vect_e(j)*ve(i,j-1,3)+(1.-edge_vect_e(j))*ve(i,j,3)
        else
             vt1(j) = edge_vect_e(j)*ve(i,j+1,1)+(1.-edge_vect_e(j))*ve(i,j,1)
             vt2(j) = edge_vect_e(j)*ve(i,j+1,2)+(1.-edge_vect_e(j))*ve(i,j,2)
             vt3(j) = edge_vect_e(j)*ve(i,j+1,3)+(1.-edge_vect_e(j))*ve(i,j,3)
        endif
       enddo
       do j=js,je
          ve(i,j,1) = vt1(j)
          ve(i,j,2) = vt2(j)
          ve(i,j,3) = vt3(j)
       enddo
     endif
! N-S edges (for u-wind):
     if ( js==1 ) then
       j = 1
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = edge_vect_s(i)*ue(i-1,j,1)+(1.-edge_vect_s(i))*ue(i,j,1)
             ut2(i) = edge_vect_s(i)*ue(i-1,j,2)+(1.-edge_vect_s(i))*ue(i,j,2)
             ut3(i) = edge_vect_s(i)*ue(i-1,j,3)+(1.-edge_vect_s(i))*ue(i,j,3)
        else
             ut1(i) = edge_vect_s(i)*ue(i+1,j,1)+(1.-edge_vect_s(i))*ue(i,j,1)
             ut2(i) = edge_vect_s(i)*ue(i+1,j,2)+(1.-edge_vect_s(i))*ue(i,j,2)
             ut3(i) = edge_vect_s(i)*ue(i+1,j,3)+(1.-edge_vect_s(i))*ue(i,j,3)
        endif
       enddo
       do i=is,ie
          ue(i,j,1) = ut1(i)
          ue(i,j,2) = ut2(i)
          ue(i,j,3) = ut3(i)
       enddo
     endif
     if ( (je+1)==npy ) then
       j = npy
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = edge_vect_n(i)*ue(i-1,j,1)+(1.-edge_vect_n(i))*ue(i,j,1)
             ut2(i) = edge_vect_n(i)*ue(i-1,j,2)+(1.-edge_vect_n(i))*ue(i,j,2)
             ut3(i) = edge_vect_n(i)*ue(i-1,j,3)+(1.-edge_vect_n(i))*ue(i,j,3)
        else
             ut1(i) = edge_vect_n(i)*ue(i+1,j,1)+(1.-edge_vect_n(i))*ue(i,j,1)
             ut2(i) = edge_vect_n(i)*ue(i+1,j,2)+(1.-edge_vect_n(i))*ue(i,j,2)
             ut3(i) = edge_vect_n(i)*ue(i+1,j,3)+(1.-edge_vect_n(i))*ue(i,j,3)
        endif
       enddo
       do i=is,ie
          ue(i,j,1) = ut1(i)
          ue(i,j,2) = ut2(i)
          ue(i,j,3) = ut3(i)
       enddo
     endif
       do j=js,je+1
          do i=is,ie
             u(i,j,k) = u(i,j,k) + dt5*( ue(i,j,1)*es(1,i,j,1) +  &
                                         ue(i,j,2)*es(2,i,j,1) +  &
                                         ue(i,j,3)*es(3,i,j,1) )
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             v(i,j,k) = v(i,j,k) + dt5*( ve(i,j,1)*ew(1,i,j,2) +  &
                                         ve(i,j,2)*ew(2,i,j,2) +  &
                                         ve(i,j,3)*ew(3,i,j,2) )
          enddo
       enddo
! Update:
      endif   ! end grid_type
 
    enddo         ! k-loop

  end subroutine update_dwinds_phys 


  subroutine update2d_dwinds_phys(is, ie, js, je, isd, ied, jsd, jed, dt, u_dt, v_dt, u, v)

! Purpose; Transform wind tendencies on A grid to D grid for the final update

  integer, intent(in):: is,  ie,  js,  je
  integer, intent(in):: isd, ied, jsd, jed
  real,    intent(in):: dt
  real, intent(inout):: u(isd:ied,  jsd:jed+1,npz)
  real, intent(inout):: v(isd:ied+1,jsd:jed  ,npz)
  real, intent(inout), dimension(isd:ied,jsd:jed,npz):: u_dt, v_dt

! local:
  real ut(isd:ied,jsd:jed)
  real:: dt5, gratio
  integer i, j, k

! Transform wind tendency on A grid to local "co-variant" components:
!$omp parallel do private (i,j,k, ut)
    do k=1,npz
       do j=js,je
          do i=is,ie
                 ut(i,j) = z11(i,j)*u_dt(i,j,k) + z12(i,j)*v_dt(i,j,k)
             v_dt(i,j,k) = z21(i,j)*u_dt(i,j,k) + z22(i,j)*v_dt(i,j,k)
             u_dt(i,j,k) = ut(i,j)
          enddo
       enddo
    enddo
! (u_dt,v_dt) are now on local coordinate system
       call timing_on('COMM_TOTAL')
  call mpp_update_domains(u_dt, v_dt, domain, gridtype=AGRID_PARAM)
       call timing_off('COMM_TOTAL')

    dt5 = 0.5 * dt

!$omp parallel do private (i,j,k, gratio)
    do k=1, npz

     if ( grid_type > 3 ) then    ! Local & one tile configurations

       do j=js,je+1
          do i=is,ie
             u(i,j,k) = u(i,j,k) + dt5*(u_dt(i,j-1,k) + u_dt(i,j,k))
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             v(i,j,k) = v(i,j,k) + dt5*(v_dt(i-1,j,k) + v_dt(i,j,k))
          enddo
       enddo

     else

!--------
! u-wind
!--------
! Edges:
    if ( js==1 ) then
       do i=is,ie
          gratio = dya(i,2) / dya(i,1)
          u(i,1,k) = u(i,1,k) + dt5*((2.+gratio)*(u_dt(i,0,k)+u_dt(i,1,k))  &
                   -(u_dt(i,-1,k)+u_dt(i,2,k)))/(1.+gratio)
       enddo
    endif

! Interior
    do j=max(2,js),min(npy-1,je+1)
       do i=is,ie
          u(i,j,k) = u(i,j,k) + dt5*(u_dt(i,j-1,k)+u_dt(i,j,k))
       enddo
    enddo

    if ( (je+1)==npy ) then
       do i=is,ie
          gratio = dya(i,npy-2) / dya(i,npy-1)
          u(i,npy,k) = u(i,npy,k) + dt5*((2.+gratio)*(u_dt(i,npy-1,k)+u_dt(i,npy,k)) &
                     -(u_dt(i,npy-2,k)+u_dt(i,npy+1,k)))/(1.+gratio)
       enddo
    endif

!--------
! v-wind
!--------
! West Edges:
    if ( is==1 ) then
       do j=js,je
          gratio = dxa(2,j) / dxa(1,j)
          v(1,j,k) = v(1,j,k) + dt5*((2.+gratio)*(v_dt(0,j,k)+v_dt(1,j,k)) &
                   -(v_dt(-1,j,k)+v_dt(2,j,k)))/(1.+gratio)
       enddo
    endif

! Interior
    do j=js,je
       do i=max(2,is),min(npx-1,ie+1)
          v(i,j,k) = v(i,j,k) + dt5*(v_dt(i-1,j,k)+v_dt(i,j,k))
       enddo
    enddo

! East Edges:
    if ( (ie+1)==npx ) then
       do j=js,je
          gratio = dxa(npx-2,j) / dxa(npx-1,j)
          v(npx,j,k) = v(npx,j,k) + dt5*((2.+gratio)*(v_dt(npx-1,j,k)+v_dt(npx,j,k)) &
                     -(v_dt(npx-2,j,k)+v_dt(npx+1,j,k)))/(1.+gratio)
       enddo
    endif

    endif   ! end grid_type

    enddo         ! k-loop

  end subroutine update2d_dwinds_phys


end module fv_update_phys_mod
