module relax_mod

   use constants_mod,      only: grav, rdgas, rvgas, RADIAN,  radius
   use mpp_domains_mod,    only: mpp_update_domains
  use mpp_mod,            only: FATAL, mpp_error, mpp_pe, stdlog, &
                                 mpp_npes, mpp_get_current_pelist, get_unit
   use time_manager_mod,   only: time_type, get_date, get_time
   use diag_manager_mod,   only: send_data
   use fv_grid_utils_mod,  only: g_sum
   use fv_diagnostics_mod, only: prt_maxmin
   use fv_timing_mod,      only: timing_on, timing_off
   
   implicit none
!-----------------------------------------------------------------------
   logical :: rf_initialized = .false.

   private
   public :: relax_tend

   namelist /relax_nml/ T_surf, T_strat, dT_ep, dT_z, p0, kappa, sigma_b, tau_ru, tau_rd, tidally_locked

   real :: T_surf = 315.0
   real :: T_strat = 200.0
   real :: dT_ep = 60.0
   real :: dT_z = 10.0
   real :: p0 = 1.e5
   real :: kappa = 0.28571428571
   real :: sigma_b = 0.7
   real :: tau_ru = 345600
   real :: tau_rd = 3456000
   logical :: tidally_locked = .true.

contains

!-----------------------------------------------------------------------
   subroutine relax_tend(npx, npy, npz, is, ie, js, je, ng, nq,   &
      pt, pe, delp, peln, pkz,   &
      t_dt, agrid,  &
      delz, hydrostatic,   &
      master, Time)

      integer, INTENT(IN   ) :: npx, npy, npz
      integer, INTENT(IN   ) :: is, ie, js, je, ng, nq
      logical, intent(IN)    :: hydrostatic
      real   , INTENT(IN   ) ::  delz(is:ie,js:je,npz)
      
      real   , INTENT(INOUT) ::   pt(is-ng:ie+  ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) :: delp(is-ng:ie+  ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) ::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
      real   , INTENT(INOUT) :: peln(is  :ie     ,1:npz+1,js   :je     )
      real   , INTENT(INOUT) ::  pkz(is  :ie     ,js   :je     ,1:npz)


! Tendencies:
      real, INTENT(INOUT):: t_dt(is:ie,js:je,npz)
      real   , INTENT(IN   ) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
      logical, INTENT(IN   ) :: master
      type(time_type), intent(in) :: Time
      !real, INTENT(IN), optional:: time_total



      real, allocatable, dimension(:,:,:)   :: t0, sigma, D_n
      real, allocatable, dimension(:,:)   :: ts0
      integer :: i,j,k
 
    integer :: f_unit, unit, ios, exists ! Namelist reading
    character(80) :: filename

      allocate(t0   (is:ie,js:je,npz))
      allocate(sigma   (is:ie,js:je,npz))
      allocate(D_n   (is:ie,js:je,npz))
      allocate(ts0   (is:ie,js:je))


      ! Read in initialisation_nml, using block altered from control.f90
      filename = "input.nml"
      inquire(file=filename,exist=exists)

    !  master = gid == 0
      if (.not. exists) then  ! This will be replaced with fv_error wrapper
         if(master) write(6,*) "file ",trim(filename)," doesn't exist"
         call mpp_error(FATAL,'FV core terminating')
      else

         f_unit = get_unit()
         open (f_unit,file=filename)
         ! Read initialisation namelist
         rewind (f_unit)
         read (f_unit,relax_nml ,iostat=ios)
         if (ios .gt. 0) then
            if(master) write(6,*) 'relax_nml ERROR: reading ',trim(filename),', iostat=',ios
            call mpp_error(FATAL,'FV core terminating')
         endif
         unit = stdlog()
         write(unit, nml=relax_nml)

         close(f_unit)

      end if


      do i = is,ie
         do j = js,je
            do k = 1,npz
               sigma(i,j,k) = pe(i,k,j) / p0

               if (sigma(i,j,k) < sigma_b) then
                  D_n(i,j,k) = 1/tau_rd
               else
                  D_n(i,j,k) = 1/tau_rd + ((1/tau_ru)-(1/tau_rd)) &
                     *((sigma(i,j,k)-sigma_b)/(1-sigma_b)) &
                     *cos(agrid(i,j,2))**4
               end if



               if (tidally_locked == .false.) then
                  t0(i,j,k) = (sigma(i,j,k))**kappa * (T_surf - dT_ep * (sin(agrid(i,j,2)) )**2 &
                     - dT_z*(sigma(i,j,k)))*(cos(agrid(i,j,2) )**2)
               else
                  t0(i,j,k) = (sigma(i,j,k))**kappa * (T_surf - dT_ep * cos(agrid(i,j,1)) * cos(agrid(i,j,2)) &
                     - dT_z*(sigma(i,j,k)))*(cos(agrid(i,j,2) )**2)
               end if

               if (t0(i,j,k) < T_strat) then
                  t0(i,j,k) = T_strat
               end if
               t_dt(i,j,k) = (t0(i,j,k) - pt(i,j,k))*D_n(i,j,k)
            end do
         end do
      end do

   end subroutine relax_tend

end module relax_mod
