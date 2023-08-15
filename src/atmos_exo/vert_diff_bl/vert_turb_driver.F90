
module vert_turb_driver_mod

!-----------------------------------------------------------------------
!
!       driver for compuing vertical diffusion coefficients
!
!         choose either:
!              1) mellor-yamada 2.5 (with tke)
!              2) non-local K scheme
!
!-----------------------------------------------------------------------
!---------------- modules ---------------------


!use      my25_turb_mod, only: my25_turb_init, my25_turb_end,  &
!                              my25_turb, tke_surf, tke

use    diffusivity_mod, only: diffusivity, molecular_diff

use   shallow_conv_mod, only: shallow_conv_init, shallow_conv

use   diag_manager_mod, only: register_diag_field, send_data

use   time_manager_mod, only: time_type, get_time, operator(-)

use      constants_mod, only: rdgas, rvgas, kappa

use            fms_mod, only: error_mesg,  check_nml_error, &
                             mpp_root_pe,      &
                              mpp_pe, input_nml_file, FATAL, stdlog, &
                              write_version_number
use                fms2_io_mod, only : file_exists, open_file, close_file
implicit none
private

!---------------- interfaces ---------------------

public   vert_turb_driver_init, vert_turb_driver_end, vert_turb_driver


!-----------------------------------------------------------------------
!--------------------- version number ----------------------------------

character(len=128) :: version = '$Id: vert_turb_driver.f90,v 1.6 2002/07/16 22:37:45 fms Exp $'
character(len=128) :: tag = '$Name: havana $'

!-----------------------------------------------------------------------
 real, parameter :: p00    = 1000.0E2
 real, parameter :: p00inv = 1./p00
 real, parameter :: d622   = rdgas/rvgas
 real, parameter :: d378   = 1.-d622
 real, parameter :: d608   = d378/d622

!---------------- private data -------------------

 logical :: do_init = .true.

 real :: gust_zi = 1000.   ! constant for computed gustiness (meters)

!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------

 logical :: do_shallow_conv  = .false.
 logical :: use_tau          = .true.
 logical :: do_molecular_diffusion = .false.

 character(len=24) :: gust_scheme  = 'constant' ! valid schemes are:
                                                !   => 'constant'
                                                !   => 'beljaars'
 real              :: constant_gust = 1.0

 namelist /vert_turb_driver_nml/ do_shallow_conv,  &
                                 gust_scheme, constant_gust, use_tau, &
                                 do_molecular_diffusion

!-------------------- diagnostics fields -------------------------------

integer :: id_tke,    id_lscale, id_lscale_0, id_z_pbl, id_gust,  &
           id_diff_t, id_diff_m, id_diff_sc, id_z_full, id_z_half,&
	   id_uwnd,   id_vwnd

real :: missing_value = -999.

character(len=9) :: mod_name = 'vert_turb'

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine vert_turb_driver (is, js, Time, Time_next, dt, frac_land,   &
                             p_half, p_full, z_half, z_full,           &
                             u_star, b_star, rough,                    &
                             u, v, t, q, um, vm, tm, qm,               &
                             udt, vdt, tdt, qdt, diff_t, diff_m, gust, &
                             mask, kbot                                )

!-----------------------------------------------------------------------
integer,         intent(in)         :: is, js
type(time_type), intent(in)         :: Time, Time_next
   real,         intent(in)         :: dt
   real, intent(in), dimension(:,:) :: frac_land,   &
                                       u_star, b_star, rough
   real, intent(in), dimension(:,:,:) :: p_half, p_full, &
                                         z_half, z_full, &
                                         u, v, t, q, um, vm, tm, qm, &
                                         udt, vdt, tdt, qdt
   real, intent(out),   dimension(:,:,:) :: diff_t, diff_m
   real, intent(out),   dimension(:,:)   :: gust
   real, intent(in),optional, dimension(:,:,:) :: mask
integer, intent(in),optional, dimension(:,:) :: kbot
!-----------------------------------------------------------------------
real   , dimension(size(t,1),size(t,2),size(t,3))   :: ape, thv
logical, dimension(size(t,1),size(t,2),size(t,3)+1) :: lmask
real   , dimension(size(t,1),size(t,2),size(t,3)+1) :: el, diag3
real   , dimension(size(t,1),size(t,2))             :: el0, z_pbl
real   , dimension(size(diff_t,1),size(diff_t,2), &
                                  size(diff_t,3))   :: diff_sc
real   , dimension(size(t,1),size(t,2),size(t,3))   :: tt, qq, uu, vv
real    :: dt_tke
integer :: ie, je, nlev, sec, day
logical :: used
!-----------------------------------------------------------------------
!----------------------- vertical turbulence ---------------------------
!-----------------------------------------------------------------------

      if (do_init)  call error_mesg  &
                     ('vert_turb_driver in vert_turb_driver_mod',  &
                      'initialization has not been called', FATAL)

     nlev = size(p_full,3)
     ie = is + size(p_full,1) - 1
     je = js + size(p_full,2) - 1

!-----------------------------------------------------------------------
!---- set up state variable used by this module ----

      if (use_tau) then
      !-- variables at time tau
          uu = u
          vv = v
          tt = t
          qq = q
      else
      !-- variables at time tau+1
          uu = um + dt*udt
          vv = vm + dt*vdt
          tt = tm + dt*tdt
          qq = qm + dt*qdt
      endif

!-----------------------------------------------------------------------

!--------------------------------------------------------------------
!----------- compute molecular diffusion, if desired  ---------------

    if (do_molecular_diffusion) then
      call molecular_diff (tt, p_half, diff_m, diff_t)
    else
      diff_m = 0.0
      diff_t = 0.0
    endif

!---------------------------
!------------------- non-local K scheme --------------

    call diffusivity ( tt, qq, uu, vv, p_full, p_half, z_full, z_half,   &
                       u_star, b_star, z_pbl, diff_m, diff_t, &
                       kbot = kbot)

!---------------------------
!-----------------------------------------------------------------------
!------------------ shallow convection ???? ----------------------------

   if (do_shallow_conv) then
        call shallow_conv (tt, qq, p_full, p_half, diff_sc, kbot)
        diff_t = diff_t + diff_sc
   endif

!-----------------------------------------------------------------------
!------------- define gustiness ------------

     if ( trim(gust_scheme) == 'constant' ) then
          gust = constant_gust
     else if ( trim(gust_scheme) == 'beljaars' ) then
!    --- from Beljaars (1994) and Beljaars and Viterbo (1999) ---
          where (b_star > 0.)
             gust = (u_star*b_star*gust_zi)**(1./3.)
          elsewhere
             gust = 0.
          endwhere
     endif

!-----------------------------------------------------------------------
!------------------------ diagnostics section --------------------------

!------- boundary layer depth -------
      if ( id_z_pbl > 0 ) then
         used = send_data ( id_z_pbl, z_pbl, Time_next, is, js )
      endif

!------- gustiness -------
      if ( id_gust > 0 ) then
         used = send_data ( id_gust, gust, Time_next, is, js )
      endif


!------- output diffusion coefficients ---------

  if ( id_diff_t > 0 .or. id_diff_m > 0 .or. id_diff_sc > 0 ) then
!       --- set up local mask for fields without surface data ---
        if (present(mask)) then
            lmask(:,:,1:nlev) = mask(:,:,1:nlev) > 0.5
            lmask(:,:,nlev+1) = .false.
        else
            lmask(:,:,1:nlev) = .true.
            lmask(:,:,nlev+1) = .false.
        endif
!       -- dummy data at surface --
        diag3(:,:,nlev+1)=0.0
  endif

!------- diffusion coefficient for heat/moisture -------
   if ( id_diff_t > 0 ) then
      diag3(:,:,1:nlev) = diff_t(:,:,1:nlev)
      used = send_data ( id_diff_t, diag3, Time_next, is, js, 1, mask=lmask )
   endif

!------- diffusion coefficient for momentum -------
   if ( id_diff_m > 0 ) then
      diag3(:,:,1:nlev) = diff_m(:,:,1:nlev)
      used = send_data ( id_diff_m, diag3, Time_next, is, js, 1, mask=lmask )
   endif

!------- diffusion coefficient for shallow conv -------
 if (do_shallow_conv) then
   if ( id_diff_sc > 0 ) then
      diag3(:,:,1:nlev) = diff_sc(:,:,1:nlev)
      used = send_data ( id_diff_sc, diag3, Time_next, is, js, 1, mask=lmask)
   endif
 endif

!--- geopotential height relative to the surface on full and half levels ----

   if ( id_z_half > 0 ) then
      !--- set up local mask for fields with surface data ---
      if ( present(mask) ) then
         lmask(:,:,1)        = .true.
         lmask(:,:,2:nlev+1) = mask(:,:,1:nlev) > 0.5
      else
         lmask = .true.
      endif
      used = send_data ( id_z_half, z_half, Time_next, is, js, 1, mask=lmask )
   endif
   
   if ( id_z_full > 0 ) then
      used = send_data ( id_z_full, z_full, Time_next, is, js, 1, rmask=mask)
   endif
   
!--- zonal and meridional wind on mass grid -------

   if ( id_uwnd > 0 ) then
      used = send_data ( id_uwnd, uu, Time_next, is, js, 1, rmask=mask)
   endif
  
   if ( id_vwnd > 0 ) then
      used = send_data ( id_vwnd, vv, Time_next, is, js, 1, rmask=mask)
   endif
  
 
   
!-----------------------------------------------------------------------

end subroutine vert_turb_driver

!#######################################################################

subroutine vert_turb_driver_init (id, jd, kd, Time)

!-----------------------------------------------------------------------
   integer,         intent(in) :: id, jd, kd
   type(time_type), intent(in) :: Time
!-----------------------------------------------------------------------
   integer, dimension(3) :: full = (/1,2,3/), half = (/1,2,4/)
   integer :: ierr, unit, io

      if (.not.do_init)  &
          call error_mesg  &
                   ('vert_turb_driver_init in vert_turb_driver_mod',  &
                    'attempting to call initialization twice', FATAL)

!-----------------------------------------------------------------------
!--------------- read namelist ------------------

      if (file_exists('input.nml')) then
         ierr=1
         read  (input_nml_file, nml=vert_turb_driver_nml, iostat=io)
         ierr = check_nml_error (io, 'vert_turb_driver_nml')
      endif

!---------- output namelist --------------------------------------------

      call write_version_number(version, tag)
      if ( mpp_pe() == mpp_root_pe() ) write (stdlog(), nml=vert_turb_driver_nml)

!     --- check namelist option ---
      if ( trim(gust_scheme) /= 'constant' .and. &
           trim(gust_scheme) /= 'beljaars' ) call error_mesg &
         ('vert_turb_driver_mod', 'invalid value for namelist &
          &variable GUST_SCHEME', FATAL)

 
!-----------------------------------------------------------------------


      if (do_shallow_conv)  call shallow_conv_init (kd)

!-----------------------------------------------------------------------

   do_init=.false.

!-----------------------------------------------------------------------

end subroutine vert_turb_driver_init


!#######################################################################

subroutine vert_turb_driver_end

!-----------------------------------------------------------------------
!      if (do_mellor_yamada) call my25_turb_end
!-----------------------------------------------------------------------

end subroutine vert_turb_driver_end

!#######################################################################

end module vert_turb_driver_mod

