module ts_twostream_mod
   use, intrinsic :: iso_fortran_env
   implicit none

   !! Precision variables
   integer, parameter :: dp = REAL64

   !! Required constants
   real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
   real(dp), parameter :: sb = 5.670374419e-8_dp

   real(dp), parameter :: D = 1.66_dp  ! Diffusivity factor

   public :: ts_twostream
   private :: lw_grey_updown, sw_grey_updown_adding

contains

        subroutine ts_twostream(surf, thick_atmosphere, nlay, nlev, Ts, Tl, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, &
   & sw_a, sw_g, lw_a, lw_g, sw_a_surf, lw_a_surf, net_F, olr, net_Fs, surf_lw_down, surf_sw_down)
      implicit none

      !! Input variables\
      logical, intent(in) :: surf, thick_atmosphere
      integer, intent(in) :: nlay, nlev
      real(dp), intent(in) :: F0, mu_z, Tint, AB, sw_a_surf, lw_a_surf, Ts
      real(dp), dimension(nlay), intent(in) :: Tl
      real(dp), dimension(nlev), intent(in) :: tau_Ve, tau_IRe
      real(dp), dimension(nlay), intent(in) :: sw_a, sw_g, lw_a, lw_g

      !! Output variables
      real(dp), intent(out) :: olr, net_Fs, surf_lw_down, surf_sw_down
      real(dp), dimension(nlev), intent(out) :: net_F

      !! Work variables
      integer :: i
      real(dp) :: Finc
      real(dp), dimension(nlay) :: bl
      real(dp), dimension(nlev) :: sw_down, sw_up, lw_down, lw_up
      real(dp), dimension(nlev) :: lw_net, sw_net

      !! Shortwave flux calculation
      if (mu_z > 0.0_dp) then
         Finc = (1.0_dp - AB) * F0
         call sw_grey_updown_adding(nlay, nlev, Finc, tau_Ve(:), mu_z, sw_a, sw_g, sw_a_surf, sw_down(:), sw_up(:), Tint)
      else
         sw_down(:) = 0.0_dp
         sw_up(:) = 0.0_dp
      end if

      call lw_grey_updown(surf, thick_atmosphere, nlay, nlev, Ts, Tl, tau_IRe(:), lw_up(:), lw_down(:), Tint)

      !! Net fluxes at each level
      lw_net(:) = lw_up(:) - lw_down(:)
      sw_net(:) = sw_up(:) - sw_down(:)
      net_F(:) = lw_net(:) + sw_net(:)

      !! Net surface flux (for surface temperature evolution)
      net_Fs = sw_down(nlev) + lw_down(nlev) - lw_up(nlev)

      !! Output olr
      olr = lw_up(1)
      surf_lw_down = lw_down(nlev)
      surf_sw_down = sw_down(nlev)

   end subroutine ts_twostream

   subroutine sw_grey_updown_adding(nlay, nlev, Finc, tau_Ve, mu_z, w_in, g_in, w_surf, sw_down, sw_up, Tint)
      implicit none

      !! Input variables
      integer, intent(in) :: nlay, nlev
      real(dp), intent(in) :: Finc, mu_z, w_surf, Tint
      real(dp), dimension(nlev), intent(in) :: tau_Ve
      real(dp), dimension(nlay), intent(in) :: w_in, g_in

      !! Output variables
      real(dp), dimension(nlev), intent(out) :: sw_down, sw_up

      !! Work variables
      integer :: k
      real(dp) :: lamtau, e_lamtau, lim, arg, apg, amg
      real(dp), dimension(nlev) ::  w, g, f
      real(dp), dimension(nlev) :: tau_Ve_s
      real(dp), dimension(nlay) :: tau
      real(dp), dimension(nlay) :: tau_s, w_s, f_s, g_s
      real(dp), dimension(nlev) :: lam, u, N, gam, alp
      real(dp), dimension(nlev) :: R_b, T_b, R, T
      real(dp), dimension(nlev) :: Tf



      sw_down(:) = Finc * exp(-tau_Ve(:))
      sw_up(:) = 0.0_dp ! We assume no upward flux here even if surface albedo


   end subroutine sw_grey_updown_adding

   subroutine lw_grey_updown(surf, thick_atmosphere, nlay, nlev, Ts, Tl, tau_IRe, lw_up, lw_down, Tint)
      implicit none

      !! Input variables
      logical, intent(in) :: surf, thick_atmosphere
      integer, intent(in) :: nlay, nlev
      real(dp), intent(in) :: Ts, Tint
      real(dp), dimension(nlev), intent(in) :: tau_IRe
      real(dp), dimension(nlay), intent(in) :: Tl
      real(dp), dimension(nlay) :: bl, lw_dtrans

      !! Output variables
      real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

      !! Work variables and arrays
      real(dp) :: b_surf
      integer :: k
      real(dp), dimension(nlay) :: dtau, Tp, B0

      if (surf) then
        b_surf = sb*Ts**4
      else
        b_surf = sb*Tint**4
      end if

      bl = sb*Tl**4

      lw_dtrans(1) = exp(-(tau_IRe(1)))
      do k = 2, nlay
         lw_dtrans(k) = exp(-(tau_IRe(k)-tau_IRe(k-1)))
      end do

      lw_down(1) = 0.0


      if (thick_atmosphere) then
                do k = 1,nlay
                        lw_down(k+1) = 2.*bl(k) * (tau_IRe(k+1) - tau_IRe(k))/(2.+ (tau_IRe(k+1) - tau_IRe(k))) + &
                                lw_down(k) * (2.- (tau_IRe(k+1) - tau_IRe(k)))/(2.+ (tau_IRe(k+1) - tau_IRe(k)))
                end do
                lw_up(nlev) = b_surf + lw_down(nlev)
      else
                do k = 1,nlay
                        lw_down(k+1) = lw_down(k)*lw_dtrans(k) + bl(k)*(1.0 - lw_dtrans(k))
                end do
                lw_up(nlev) = b_surf
      end if


      if (thick_atmosphere) then
                do k = nlay,1,-1
                        lw_up(k) = 2.*bl(k) * -1.*(tau_IRe(k) - tau_IRe(k+1))/ (2.- (tau_IRe(k) - tau_IRe(k+1))) &
                            + lw_up(k+1)* (2. + (tau_IRe(k) - tau_IRe(k+1)))/(2. - (tau_IRe(k) - tau_IRe(k+1))) 
                end do
      else
                do k = nlay,1,-1
                        lw_up(k) = lw_up(k+1)*lw_dtrans(k) + bl(k)*(1.0 - lw_dtrans(k))
                end do
      end if

   end subroutine lw_grey_updown

end module ts_twostream_mod
