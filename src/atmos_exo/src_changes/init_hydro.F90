! $Id: init_hydro.F90,v 17.0 2009/07/21 02:52:36 fms Exp $

module init_hydro_mod


      use constants_mod, only: grav, rdgas, stefan, kappa
      use fv_grid_utils_mod,    only: g_sum
      use fv_grid_tools_mod,    only: area
      use fv_mp_mod,        only: gid, masterproc
      use mpp_mod,          only: FATAL, mpp_error, mpp_pe, stdlog, &
                                  mpp_npes, mpp_get_current_pelist, get_unit
!     use fv_diagnostics_mod, only: prt_maxmin

      implicit none
      private

      public :: p_var, hydro_eq


    ! HI added:
    !---------------------------------------------------------------------------------------
    ! Namelist values to read in 
    !---------------------------------------------------------------------------------------
    
    ! Default values of namelist arguments
    integer :: init_type = 1   ! Select initialisation type:
                               ! 1: Isothermal atmosphere
                               ! 2: Guillot 2010 radiative equilibrium profile
                               ! 3: Dry adiabat initialisation
    real    :: ts = 300.       ! Surface temperature
    real    :: T_strat = 200.  ! Minimum stratospheric temperature. If init_type = 1, this is
                               ! the atmospheric temperature

    ! Required for pressure level setup
    real :: surface_pressure
  
    ! Required for Guillot 2010 Initialisation
    real  :: S = 1360.       ! Stellar constant
    real  :: f = 1./4.       ! Flux redistribution constant (1/2 = dayside, 1/4 = global redist.)
    real  :: tau_sw, tau_lw  ! Total SW and LW fluxes for the atmosphere  
    real  :: Tint            ! Internal temperature
    

    namelist /initialisation_nml/ init_type, ts, T_strat, S, f, tau_sw, tau_lw, Tint, &
                                surface_pressure


contains

!-------------------------------------------------------------------------------
 subroutine p_var(km, ifirst, ilast, jfirst, jlast, ptop, ptop_min,    &
                  delp, delz, pt, ps,  pe, peln, pk, pkz, cappa, q, ng, nq,    &
                  dry_mass, adjust_dry_mass, mountain, moist_phys,      &
                  hydrostatic, ktop, nwat, make_nh)
               
! Given (ptop, delp) computes (ps, pk, pe, peln, pkz)
! Input:
   integer,  intent(in):: km
   integer,  intent(in):: ifirst, ilast            ! Longitude strip
   integer,  intent(in):: jfirst, jlast            ! Latitude strip
   integer,  intent(in):: nq, nwat
   integer,  intent(in):: ng
   integer,  intent(in):: ktop
   logical, intent(in):: adjust_dry_mass, mountain, moist_phys, hydrostatic
   real, intent(in):: dry_mass, cappa, ptop, ptop_min
   real, intent(in   )::   pt(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng, km)
   real, intent(inout):: delz(ifirst   :ilast   ,jfirst   :jlast   , km)
   real, intent(inout):: delp(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng, km)
   real, intent(inout)::    q(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng, km, nq)
   logical, optional:: make_nh
! Output:
   real, intent(out) ::   ps(ifirst-ng:ilast+ng, jfirst-ng:jlast+ng)
   real, intent(out) ::   pk(ifirst:ilast, jfirst:jlast, km+1)
   real, intent(out) ::   pe(ifirst-1:ilast+1,km+1,jfirst-1:jlast+1) ! Ghosted Edge pressure
   real, intent(out) :: peln(ifirst:ilast, km+1, jfirst:jlast)    ! Edge pressure
   real, intent(out) ::  pkz(ifirst:ilast, jfirst:jlast, km)

! Local
   real ratio(ifirst:ilast)
   real pek, lnp, ak1, rdg, dpd
   integer i, j, k


! Check dry air mass & compute the adjustment amount:
   call drymadj(km, ifirst, ilast,  jfirst,  jlast, ng, cappa, ptop, ps, &
                delp, q, nq, nwat, dry_mass, adjust_dry_mass, moist_phys, dpd)

   pek = ptop ** cappa

   do j=jfirst,jlast
      do i=ifirst,ilast
         pe(i,1,j) = ptop
         pk(i,j,1) = pek
      enddo

      if ( adjust_dry_mass ) then
         do i=ifirst,ilast
            ratio(i) = 1. + dpd/(ps(i,j)-ptop)
         enddo 
         do k=1,km
            do i=ifirst,ilast
               delp(i,j,k) = delp(i,j,k) * ratio(i)
            enddo
         enddo
      endif

      do k=2,km+1
         do i=ifirst,ilast
            pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
            peln(i,k,j) = log(pe(i,k,j))
            pk(i,j,k) = exp( cappa*peln(i,k,j) )
!            pk(i,j,k) = pe(i,k,j)**cappa
         enddo
      enddo

      do i=ifirst,ilast
         ps(i,j) = pe(i,km+1,j)
      enddo

      if( ptop < ptop_min ) then
!---- small ptop modification -------------
          ak1 = (cappa + 1.) / cappa
          do i=ifirst,ilast
             peln(i,1,j) = peln(i,2,j) - ak1
          enddo
      else
             lnp = log( ptop )
          do i=ifirst,ilast
             peln(i,1,j) = lnp
          enddo
      endif

      if ( hydrostatic ) then
         do k=1,km
            do i=ifirst,ilast
               pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(cappa*(peln(i,k+1,j)-peln(i,k,j)))
            enddo
         enddo
      endif
   enddo

!--------
! Caution:
!------------------------------------------------------------------
! The following form is the same as in "fv_update_phys.F90"
! Therefore, restart reproducibility is only enforced in diabatic cases
!------------------------------------------------------------------
! For adiabatic runs, this form is not exactly the same as in mapz_module;
! Therefore, rounding differences will occur with restart!

   if ( .not.hydrostatic ) then

      if ( ktop>1 ) then
! Compute pkz using hydrostatic formular:
         do k=1,ktop-1
            do j=jfirst,jlast
            do i=ifirst,ilast
               pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(cappa*(peln(i,k+1,j)-peln(i,k,j)))
            enddo
            enddo
         enddo
      endif

      rdg = -rdgas / grav
      if ( present(make_nh) ) then
          if ( make_nh ) then
             do k=1,km
                do j=jfirst,jlast
                   do i=ifirst,ilast
                      delz(i,j,k) = rdg*pt(i,j,k)*(peln(i,k+1,j)-peln(i,k,j))
                   enddo
                enddo
             enddo
             if(gid==0) write(*,*) 'delz computed from hydrostatic state'
!            call prt_maxmin('delz', delz, ifirst, ilast, jfirst, jlast, 0, km, 1., gid==masterproc)
          endif
      endif

      do k=ktop,km
         do j=jfirst,jlast
            do i=ifirst,ilast
               pkz(i,j,k) = (rdg*delp(i,j,k)*pt(i,j,k)/delz(i,j,k))**cappa
            enddo
         enddo
      enddo
   endif

 end subroutine p_var



 subroutine drymadj(km,  ifirst, ilast, jfirst,  jlast,  ng, &  
                    cappa,   ptop, ps, delp, q,  nq,  nwat,  &
                    dry_mass, adjust_dry_mass, moist_phys, dpd)

! !INPUT PARAMETERS:
      integer km
      integer ifirst, ilast  ! Long strip
      integer jfirst, jlast  ! Latitude strip    
      integer nq, ng, nwat
      real, intent(in):: dry_mass
      real, intent(in):: ptop
      real, intent(in):: cappa
      logical, intent(in):: adjust_dry_mass
      logical, intent(in):: moist_phys

! !INPUT/OUTPUT PARAMETERS:     
      real, intent(in)::   q(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng,km,nq)
      real, intent(in)::delp(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng,km)     !
      real, intent(inout):: ps(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng)        ! surface pressure
      real, intent(out):: dpd
! Local
      real  psd(ifirst:ilast,jfirst:jlast)     ! surface pressure  due to dry air mass
      real  psmo, psdry
      integer i, j, k

!$omp parallel do default(shared) private(i, j, k)
      do j=jfirst,jlast

         do i=ifirst,ilast
             ps(i,j) = ptop
            psd(i,j) = ptop
         enddo

         do k=1,km
            do i=ifirst,ilast
               ps(i,j) = ps(i,j) + delp(i,j,k)
            enddo
         enddo

       if ( nwat>=1 ) then
          do k=1,km
             do i=ifirst,ilast
                psd(i,j) = psd(i,j) + delp(i,j,k)*(1. - sum(q(i,j,k,1:nwat)))
             enddo
          enddo
        else
          do i=ifirst,ilast
             psd(i,j) = ps(i,j)
          enddo
        endif
      enddo

! Check global maximum/minimum
#ifndef QUICK_SUM
      psdry = g_sum(psd, ifirst, ilast, jfirst, jlast, ng, area, 1, .true.) 
       psmo = g_sum(ps(ifirst:ilast,jfirst:jlast), ifirst, ilast, jfirst, jlast,  &
                     ng, area, 1, .true.) 
#else
      psdry = g_sum(psd, ifirst, ilast, jfirst, jlast, ng, area, 1) 
       psmo = g_sum(ps(ifirst:ilast,jfirst:jlast), ifirst, ilast, jfirst, jlast,  &
                     ng, area, 1) 
#endif

      if(gid==masterproc) then
         write(6,*) 'Total surface pressure (mb) = ', 0.01*psmo
         if ( moist_phys ) then
              write(6,*) 'mean dry surface pressure = ', 0.01*psdry
              write(6,*) 'Total Water (kg/m**2) =', real(psmo-psdry,4)/GRAV
         endif
      endif

      if( adjust_dry_mass ) Then
          dpd = real(dry_mass - psdry,4)
          if(gid==masterproc) write(6,*) 'dry mass to be added (pascals) =', dpd
      endif

 end subroutine drymadj



 subroutine hydro_eq(km, is, ie, js, je, ps, hs, drym, delp, ak, bk,  &
                     pt, tsurf, delz, ng, mountain, hybrid_z, agrid)
! Input: 
  integer, intent(in):: is, ie, js, je, km, ng
  real, intent(in):: ak(km+1), bk(km+1)
  real, intent(in):: hs(is-ng:ie+ng,js-ng:je+ng)
  real, intent(in):: drym
  logical, intent(in):: mountain
  logical, intent(in):: hybrid_z
  real, intent(in) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
! Output
  real, intent(out):: ps(is-ng:ie+ng,js-ng:je+ng)
  real, intent(out)::   pt(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(out)::   tsurf(is-ng:ie+ng,js-ng:je+ng)
  real, intent(out):: delp(is-ng:ie+ng,js-ng:je+ng,km)
  real, intent(inout):: delz(is:ie,js:je,km)
! Local
  real   gz(is:ie,km+1)
  real   ph(is:ie,km+1)

  real   pf(km)

  real mslp, z1, t1, p1, t0, a0, psm
  real ztop, c0
#ifdef INIT_4BYTE
  real(kind=4) ::  dps 
#else
  real dps    ! note that different PEs will get differt dps during initialization
              ! this has no effect after cold start
#endif
  real p0, gztop, ptop

  integer  i,j,k, f_unit, unit, ios, exists
  logical master
  character(80) :: filename

   ! HI Added following:

   ! Read in initialisation_nml, using block altered from control.f90
    filename = "input.nml"
    inquire(file=filename,exist=exists)

    master = gid == 0
    if (.not. exists) then  ! This will be replaced with fv_error wrapper
       if(master) write(6,*) "file ",trim(filename)," doesn't exist"
       call mpp_error(FATAL,'FV core terminating')
    else

       f_unit = get_unit()
       open (f_unit,file=filename)
       ! Read initialisation namelist
       rewind (f_unit)
       read (f_unit,initialisation_nml ,iostat=ios)
       if (ios .gt. 0) then
          if(master) write(6,*) 'initialisation_nml ERROR: reading ',trim(filename),', iostat=',ios
          call mpp_error(FATAL,'FV core terminating')
       endif
       unit = stdlog()
       write(unit, nml=initialisation_nml)

       close(f_unit)

    end if

    ! end

  if ( gid==masterproc ) write(*,*) 'Initializing ATM hydrostatically'


  if ( gid==masterproc ) write(*,*) 'Initializing Earth'
! Given p1 and z1 (250mb, 10km)
        p1 = 25000.
        z1 = 10.E3 * grav
        t1 = 215.
        t0 = 280.            ! sea-level temp.
        a0 = (t1-t0)/z1
        c0 = t0/a0

     if ( hybrid_z ) then
          ptop = 100.   ! *** hardwired model top *** 
     else
          ptop = ak(1)
     endif

     ztop = z1 + (rdgas*t1)*log(p1/ptop)
     if(gid==masterproc) write(6,*) 'ZTOP is computed as', ztop/grav*1.E-3

     !mslp = surface_pressure
 !    mslp = 200000.E2

     do j=js,je
        do i=is,ie
           ps(i,j) = surface_pressure
        enddo
     enddo
     dps = 0.

  do j=js,je
     do i=is,ie
        ps(i,j) = ps(i,j) + dps
        gz(i,   1) = ztop
        gz(i,km+1) = hs(i,j)
        ph(i,   1) = ptop                                                     
        ph(i,km+1) = ps(i,j)                                               
     enddo

     if ( hybrid_z ) then
!---------------
! Hybrid Z
!---------------
        do k=km,2,-1
           do i=is,ie
              gz(i,k) = gz(i,k+1) - delz(i,j,k)*grav 
           enddo
        enddo
! Correct delz at the top:
        do i=is,ie
            delz(i,j,1) = (gz(i,2) - ztop) / grav
        enddo
 
        do k=2,km
           do i=is,ie
              if ( gz(i,k) >= z1 ) then
! Isothermal
                 ph(i,k) = ptop*exp( (gz(i,1)-gz(i,k))/(rdgas*t1) )
              else
! Constant lapse rate region (troposphere)
                 ph(i,k) = ps(i,j)*((hs(i,j)+c0)/(gz(i,k)+c0))**(1./(a0*rdgas))
              endif
           enddo
        enddo
     else
!---------------
! Hybrid sigma-p
!---------------
       do k=2,km+1
          do i=is,ie
             ph(i,k) = ak(k) + bk(k)*ps(i,j)
          enddo
       enddo

       do k=2,km
          do i=is,ie
             if ( ph(i,k) <= p1 ) then
! Isothermal
                 gz(i,k) = ztop + (rdgas*t1)*log(ptop/ph(i,k))
             else
! Constant lapse rate region (troposphere)
                 gz(i,k) = (hs(i,j)+c0)/(ph(i,k)/ps(i,j))**(a0*rdgas) - c0
             endif
          enddo
       enddo
     endif  ! end hybrid_z

! Convert geopotential to Temperature
      do k=1,km
         do i=is,ie
              pt(i,j,k) = (gz(i,k)-gz(i,k+1))/(rdgas*(log(ph(i,k+1)/ph(i,k))))
              pt(i,j,k) = max(t1, pt(i,j,k))
              delp(i,j,k) = ph(i,k+1) - ph(i,k)
              pf(k) = (ph(i,k+1) - ph(i,k))/(log(ph(i,k+1)) - log(ph(i,k)))
         enddo
      enddo
      if (gid==masterproc .and. j==js) then
         i = is
         do k=1,km
            write(*,*) k, pt(i,j,k), gz(i,k+1), (gz(i,k)-gz(i,k+1)), ph(i,k)
         enddo
      endif

   enddo    ! j-loop

   if ( hybrid_z ) then 
!      call prt_maxmin('INIT_hydro: delz', delz, is, ie, js, je,  0, km, 1., gid==masterproc)
!      call prt_maxmin('INIT_hydro: DELP', delp, is, ie, js, je, ng, km, 1., gid==masterproc)
   endif
!  call prt_maxmin('INIT_hydro: PT  ', pt,   is, ie, js, je, ng, km, 1., gid==masterproc)

   ! HI Added following:


    call tp_init(is-ng, ie+ng, js-ng, je+ng, ng, km, pf, pt, agrid)
    

   do j=js,je
      do i=is,ie
         tsurf(i,j) = pt(i,j,km)
      enddo
   enddo

 end subroutine hydro_eq


 subroutine tp_init(is, ie, js, je, ng, num_levels, pf, Tf, agrid)
    !---------------------------------------------------------------------------------------
    ! Purpose: Initialising temperature grid at beginning of run                          
    ! --------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------
    ! Input arguments
    !---------------------------------------------------------------------------------------
    integer, intent(in) :: is,ie,js,je,ng,num_levels
    real,    intent(in) :: pf(num_levels)
    real, intent(in) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
    !---------------------------------------------------------------------------------------
    ! Output arguments
    !---------------------------------------------------------------------------------------
    real, intent(out) :: Tf(is:ie, js:je, num_levels)


    !---------------------------------------------------------------------------------------
    ! Local arguments
    !---------------------------------------------------------------------------------------
    
    integer :: i,j,k
    logical :: master
    character(80) :: filename

    ! Required for Guillot 2010 setup
    real, dimension(num_levels) :: tau ! Optical depth = tau_lw*(p/ps)
    real  :: gamma, mu, lat_tl 
    real  :: Tirr ! Irradiation temperature
    
    !---------------------------------------------------------------------------------------
    ! Main body of function
    !---------------------------------------------------------------------------------------

    ! Select type of initialisation profile
    select case(init_type)
       
    ! Isothermal profile
    case(1)
       Tf(is:ie,js:je,1:num_levels) = T_strat

    ! Guillot 2010, equation 29
    case(2)
       Tirr = (S/stefan)**0.25
       tau = pf/surface_pressure * tau_lw
       gamma = tau_sw/tau_lw

       do k=1,num_levels
          Tf(is:ie,js:je,k) = ( 3./4.*Tint**4.*(2./3. + tau(k)) + &
               3./4.*Tirr**4*f*(2./3. + 1./gamma/sqrt(3.) + &
               (gamma/sqrt(3.) - 1/gamma/sqrt(3.))*exp(-gamma*tau(k)*sqrt(3.)) ) ) ** (0.25)
      enddo

    ! Dry adiabat
    case(3)
       do k=1, num_levels
          do j=js,je
             do i=is,ie
                Tf(i,j,k) = ts* (pf(k)/pf(num_levels)) ** kappa
                Tf(i,j,k) = max(Tf(i,j,k), T_strat)
             end do
          end do
       end do

    ! Guillot 2010, equation 27, + Koll 2015 equation B3
    case(4)
       Tirr = (S/stefan)**0.25
       tau = pf/surface_pressure * tau_lw
       gamma = tau_sw/tau_lw
       do i=is,ie
       do j=js,je
       lat_tl = asin(cos(agrid(i,j,2))*cos(agrid(i,j,1)))
       mu = cos(lat_tl)
       do k=1,num_levels
          Tf(i,j,k) = ( 3./4.*Tint**4.*(2./3. + tau(k)) + &
               3./4.*Tirr**4*mu*(2./3. + mu/gamma + &
                  ((gamma/(3*mu)) - (mu/gamma)) * exp( -gamma*tau(k)/mu) ) )** (0.25)
       enddo
       enddo
       enddo

       
    end select

  end subroutine tp_init

end module init_hydro_mod
