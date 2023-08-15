module simple_boundary_mod

 use constants_mod,      only: grav, rdgas, rvgas, RADIAN, kappa, radius
 use mpp_domains_mod,    only: mpp_update_domains
  use mpp_mod,            only: FATAL, mpp_error, mpp_pe, stdlog, &
                                 mpp_npes, mpp_get_current_pelist
 use time_manager_mod,   only: time_type, get_date, get_time
 use fv_grid_utils_mod,  only: g_sum
 use fv_mp_mod,          only: is_master

      implicit none
!-----------------------------------------------------------------------
      logical :: rf_initialized = .false.

      private
      public :: simple_boundary_tend, simple_boundary_init

  namelist /simple_boundary_nml/ sigma_b, p0, tau_f, mix_coeff_surf, mix_coeff_atmos

      real :: sigma_b = 0.7
      real :: p0 = 1.e5
      real :: tau_f = 86400
      real :: mix_coeff_atmos = 1.e6
      real :: mix_coeff_surf = 1.e6

contains

  subroutine simple_boundary_init

    integer :: f_unit, unit, ios, exists ! Namelist reading
    character(80) :: filename

    filename = "input.nml"
    inquire(file=filename,exist=exists)

!    master = gid == 0
    if (.not. exists) then  ! This will be replaced with fv_error wrapper
       if(is_master()) write(6,*) "file ",trim(filename)," doesn't exist"
       call mpp_error(FATAL,'FV core terminating')
    else

       open (newunit=f_unit,file=filename)
       ! Read initialisation namelist
       rewind (f_unit)
       read (f_unit,simple_boundary_nml ,iostat=ios)
       if (ios .gt. 0) then
          if(is_master()) write(6,*) 'simple_boundary_nml ERROR: reading ',trim(filename),', iostat=',ios
          call mpp_error(FATAL,'FV core terminating')
       endif
       unit = stdlog()
       write(unit, nml=simple_boundary_nml)

       close(f_unit)

    end if

  end subroutine simple_boundary_init
!-----------------------------------------------------------------------
 subroutine simple_boundary_tend(npx, npy, npz, is, ie, js, je, ng, nq,   &
                              pt, ts, ua, va, pe, delp, peln, pkz,   &
                              t_dt, ts_dt, u_dt, v_dt, agrid,  &
                              delz, hydrostatic,   &
                               master, Time)

      integer, INTENT(IN   ) :: npx, npy, npz
      integer, INTENT(IN   ) :: is, ie, js, je, ng, nq
      logical, intent(IN)    :: hydrostatic
      real   , INTENT(IN   ) ::  delz(is:ie,js:je,npz)

      real   , INTENT(INOUT) ::   ua(is-ng:ie+  ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) ::   va(is-ng:ie+  ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) ::   pt(is-ng:ie+  ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) ::   ts(is-ng:ie+  ng,js-ng:je+  ng)
      real   , INTENT(INOUT) :: delp(is-ng:ie+  ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) ::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
      real   , INTENT(INOUT) :: peln(is  :ie     ,1:npz+1,js   :je     )
      real   , INTENT(INOUT) ::  pkz(is  :ie     ,js   :je     ,1:npz)


! Tendencies:
      real, INTENT(INOUT):: u_dt(is:ie,js:je,npz)
      real, INTENT(INOUT):: v_dt(is:ie,js:je,npz)
      real, INTENT(INOUT):: t_dt(is:ie,js:je,npz)
      real, INTENT(INOUT):: ts_dt(is:ie,js:je)
      real   , INTENT(IN   ) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
      logical, INTENT(IN   ) :: master
      type(time_type), intent(in) :: Time



      real, allocatable, dimension(:,:,:)   :: sigma, D_r
      integer :: i,j,k
    

allocate(sigma   (is:ie,js:je,npz))
allocate(D_r   (is:ie,js:je,npz))

!$OMP parallel do default(none) shared(npz,is,ie,js,je,sigma,pe,p0,sigma_b,D_r,tau_f,u_dt,ua,v_dt,va)
do k = 1,npz
do i = is,ie
do j = js,je

sigma(i,j,k) = pe(i,k,j) / p0

if (sigma(i,j,k) < sigma_b) then
        D_r(i,j,k) = 0
else
        D_r(i,j,k) = (1/tau_f) * (sigma(i,j,k) - sigma_b)/(1-sigma_b)
end if

u_dt(i,j,k) = - ua(i,j,k) * D_r(i,j,k)
v_dt(i,j,k) = - va(i,j,k) * D_r(i,j,k)
end do
end do
end do

!$OMP parallel do default(none) shared(is,ie,js,je,npz,t_dt, ts, ts_dt,pt, mix_coeff_atmos, mix_coeff_surf)
 do j = js,je
    do i = is,ie
       t_dt(i,j,npz) = t_dt(i,j,npz) + (ts(i,j) - pt(i,j,npz))/mix_coeff_atmos
                ts_dt(i,j) = ts_dt(i,j) - (ts(i,j) - pt(i,j,npz))/mix_coeff_surf
            end do
        end do

 end subroutine simple_boundary_tend

 end module simple_boundary_mod
