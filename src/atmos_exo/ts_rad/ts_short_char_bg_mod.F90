!!!
! Elspeth KH Lee - May 2021 : Initial version
!                - Oct 2021 : adding method & Bezier interpolation
! sw: Adding layer method with scattering
! lw: Two-stream method following the short characteristics method (e.g. Helios-r2: Kitzmann et al. 2018)
!     Uses the method of short characteristics (Olson & Kunasz 1987) with linear interpolants.
!     Pros: Very fast, accurate at high optical depths, very stable
!     Cons: No lw scattering
!!!

module ts_short_char_bg_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Gauss quadrature variables, cosine angle values (uarr] and weights (w)
  !! here you can comment in/out groups of mu values for testing
  !! make sure to make clean and recompile if you change these

  !! single angle diffusion factor approximation - typically 1/1.66
  !integer, parameter :: nmu = 1
  !real(dp), dimension(nmu), parameter :: uarr = (/1.0_dp/1.66_dp/)
  !real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)
  !real(dp), dimension(nmu), parameter :: wuarr = uarr * w


  !! Legendre quadrature for 2 nodes
  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.21132487_dp, 0.78867513_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! Lacis & Oinas (1991) 3 point numerical values - Does not work somehow, e-mail me if you know why :)
  ! integer, parameter :: nmu = 3
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.1_dp, 0.5_dp, 1.0_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = (/0.0433_dp, 0.5742_dp, 0.3825_dp/)

  !! Legendre quadrature for 4 nodes
  ! integer, parameter :: nmu = 4
  ! real(dp), dimension(nmu), parameter :: uarr = &
  !   & (/0.06943184_dp, 0.33000948_dp, 0.66999052_dp, 0.93056816_dp/)
  ! real(dp), dimension(nmu), parameter :: w = &
  !   & (/0.17392742_dp, 0.32607258_dp, 0.32607258_dp, 0.17392742_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! 5 point EGP quadrature values
  ! integer, parameter :: nmu = 5
  ! real(dp), dimension(nmu), parameter :: uarr = &
  !   &(/0.0985350858_dp, 0.3045357266_dp, 0.5620251898_dp, 0.8019865821_dp, 0.9601901429_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = &
  !   & (/0.0157479145_dp, 0.0739088701_dp, 0.1463869871_dp, 0.1671746381_dp, 0.0967815902_dp/)

  public :: ts_short_char_bg
  private :: bg_diffuse_updown_linear, bg_direct_beam_updown_adding, linear_log_interp, calc_cff, bezier_interp, BB_integrate

contains

  subroutine ts_short_char_bg(surf, Bezier, Planck_table, n_bands, nlay, nlev, Ts, Tl, pl, pe, tau_Ve, tau_IRe, tau_bg, mu_z, F0, &
                              Tint, AB, sw_a, sw_g, lw_a, lw_g, sw_a_surf, lw_a_surf, net_F, diffuse_up, diffuse_down, direct_up, &
                              direct_down, opr, net_Fs, cff, scff)
    implicit none

    !! Input variables
    logical, intent(in) :: surf, Bezier
    integer, intent(in) :: n_bands, nlay, nlev
    real(dp), intent(in) :: mu_z, Tint, AB, Ts, sw_a_surf, lw_a_surf
    real(dp), dimension(n_bands), intent(in) :: F0
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(n_bands+1,2999), intent(in) :: Planck_table
    real(dp), dimension(nlev), intent(in) :: tau_Ve, tau_IRe
    real(dp), dimension(n_bands,nlev), intent(in) :: tau_bg
    real(dp), dimension(nlay), intent(in) :: sw_a, sw_g, lw_a, lw_g

    !! Output variables
    real(dp), dimension(n_bands), intent(out) :: opr, net_Fs, scff
    real(dp), dimension(n_bands,nlev), intent(out) :: net_F
    real(dp), dimension(n_bands,nlev), intent(out) :: diffuse_up, diffuse_down, direct_up, direct_down
    real(dp), dimension(n_bands,nlay), intent(out) :: cff

    !! Work variables
    integer :: i, b, T_closest
    real(dp) :: Finc, be_int
    real(dp), dimension(nlev) :: Te, be
    real(dp), dimension(n_bands,nlev) :: be_bg, diffuse_net, direct_net
    real(dp), dimension(2999) :: T_range
    real(dp), dimension(n_bands,2) :: wn_edges

    !! Find temperature at layer edges through interpolation and extrapolation
    if (Bezier .eqv. .True.) then
      ! Perform interpolation using Bezier peicewise polynomial interpolation
      do i = 2, nlay-1
        call bezier_interp(pl(i-1:i+1), Tl(i-1:i+1), 3, pe(i), Te(i))
        !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
      end do
      call bezier_interp(pl(nlay-2:nlay), Tl(nlay-2:nlay), 3, pe(nlay), Te(nlay))
    else
      ! Perform interpolation using linear interpolation
      do i = 2, nlay
        call linear_log_interp(pe(i), pl(i-1), pl(i), Tl(i-1), Tl(i), Te(i))
        !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
      end do
    end if

    ! Edges are linearly interpolated
    Te(1) = 10.0_dp**(log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))

    if (surf .eqv. .True.) then
      Te(nlev) = Ts ! If surface, temperature at edge is the surface temperature
    else
      Te(nlev) = 10.0_dp**(log10(Tl(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))
    end if


    !! Shortwave flux calculation
    if (mu_z > 0.0_dp) then
      if (n_bands .eq. 1) then
        Finc = (1.0_dp - AB) * F0(1)
        call bg_direct_beam_updown_adding(nlay, nlev, Finc, tau_bg(2,:), mu_z, sw_a, sw_g, sw_a_surf, direct_down(1,:), direct_up(1,:))
      else if (n_bands .ge. 4) then
        do b = 1, n_bands
          Finc = (1.0_dp - AB) * F0(b)
          !print*, "b, (1.0_dp - AB), F0(b) = ", b, (1.0_dp - AB), F0(b)
          call bg_direct_beam_updown_adding(nlay, nlev, Finc, tau_bg(b,:), mu_z, sw_a, sw_g, sw_a_surf, direct_down(b,:), direct_up(b,:))
          !print*, "b, Finc = ", b, Finc
          !print*, "b, direct_down(b,:) = ", b, direct_down(b,:)
          !print*, "b, sum(direct_down(b,:)) = ", b, sum(direct_down(b,:))
        end do
      end if
      !print*, "======================================"
    else
      direct_down(:,:) = 0.0_dp
      direct_up(:,:)   = 0.0_dp
    end if

    !! Longwave two-stream flux calculation
    be(:) = (sb * Te(:)**4)/pi  ! Integrated planck function intensity at levels
    if (surf .eqv. .True.) then
      be_int = 0.0_dp
    else
      be_int = (sb * Tint**4)/pi ! Integrated planck function intensity for internal temperature
    end if

    if (n_bands .eq. 1) then
      call bg_diffuse_updown_linear(surf, nlay, nlev, be, be_int, tau_bg(1,:), lw_a_surf, diffuse_up(1,:), diffuse_down(1,:))

      !! Net fluxes at each level
      diffuse_net(1,:) = diffuse_up(1,:) - diffuse_down(1,:)
      direct_net(1,:)  = direct_up(1,:) - direct_down(1,:)
      net_F(1,:) = diffuse_net(1,:) + direct_net(1,:)

      !! Net surface flux (for surface temperature evolution)
      !! We have to define positive as downward (heating) and cooling (upward) in this case
      net_Fs(1) = direct_down(1,nlev) + diffuse_down(1,nlev) - diffuse_up(1,nlev)

      !! Output the olr
      opr(1) = diffuse_up(1,1)

      call calc_cff(surf, nlay, be, pe, tau_bg(1,:), cff(1,:), scff(1))

    else if (n_bands .ge. 4) then

      T_range = Planck_table(1,:)
      !! Fill
      do i = 1,nlev
        ! be cannot be zero. lw_grey_updown_linear divides by be(k) to get Am(k)
        T_closest = minloc(abs(T_range-(T_range/T_range)*Te(i)), DIM=1) ! returns the index of the closest element to Te(i) in TL
        be_bg(1,i) = max(Planck_table(2,T_closest), 1d-30)
        be_bg(2,i) = max(Planck_table(3,T_closest), 1d-30)
        be_bg(3,i) = max(Planck_table(4,T_closest), 1d-30)
        if (n_bands .eq. 4) then
          be_bg(4,i) = max(Planck_table(5,T_closest), 1d-30)
        elseif (n_bands .eq. 5) then
          be_bg(4,i) = max(Planck_table(5,T_closest), 1d-30)
          be_bg(5,i) = max(Planck_table(6,T_closest), 1d-30)
        elseif (n_bands .eq. 6) then
          be_bg(4,i) = max(Planck_table(5,T_closest), 1d-30)
          be_bg(5,i) = max(Planck_table(6,T_closest), 1d-30)
          be_bg(6,i) = max(Planck_table(7,T_closest), 1d-30)
        end if
      end do

      do b = 1, n_bands
        call bg_diffuse_updown_linear(surf, nlay, nlev, be_bg(b,:), be_int, tau_bg(b,:), lw_a_surf, diffuse_up(b,:), diffuse_down(b,:))

        !! Net fluxes at each level
        diffuse_net(b,:) = diffuse_up(b,:) - diffuse_down(b,:)
        direct_net(b,:)  = direct_up(b,:) - direct_down(b,:)
        net_F(b,:) = diffuse_net(b,:) + direct_net(b,:)

        !! Net surface flux (for surface temperature evolution)
        !! We have to define positive as downward (heating) and cooling (upward) in this case
        net_Fs(b) = direct_down(b,nlev) + diffuse_down(b,nlev) - diffuse_up(b,nlev)

        !! Output the olr
        opr(b) = diffuse_up(b,1)

        !! Output the contribution function
        call calc_cff(surf, nlay, be_bg(b,:), pe, tau_bg(b,:), cff(b,:), scff(b))

      end do
    end if

  end subroutine ts_short_char_bg


  subroutine bg_direct_beam_updown_adding(nlay, nlev, Finc, tau_Ve, mu_z, w_in, g_in, w_surf, sw_down, sw_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: Finc, mu_z, w_surf
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

    ! Design w and g to include surface property level
    w(1:nlay) = w_in(:)
    g(1:nlay) = g_in(:)

    w(nlev) = 0.0_dp
    g(nlev) = 0.0_dp

    ! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w(:) == 0.0_dp)) then
      sw_down(:) = Finc * mu_z * exp(-tau_Ve(:)/mu_z)
      sw_down(nlev) = sw_down(nlev) * (1.0_dp - w_surf) ! The surface flux for surface heating is the amount of flux absorbed by surface
      sw_up(:) = 0.0_dp ! We assume no upward flux here even if surface albedo
      return
    end if

    w(nlev) = w_surf
    g(nlev) = 0.0_dp

    ! Backscattering approximation
    f(:) = g(:)**2

    !! Do optical depth rescaling
    tau_Ve_s(1) = tau_Ve(1)
    do k = 1, nlay
      tau(k) = tau_Ve(k+1) - tau_Ve(k)
      tau_s(k) = tau(k) * (1.0_dp - w(k)*f(k))
      tau_Ve_s(k+1) = tau_Ve_s(k) + tau_s(k)
    end do

    do k = 1, nlev

      w_s(k) = w(k) * ((1.0_dp - f(k))/(1.0_dp - w(k)*f(k)))
      g_s(k) = (g(k) - f(k))/(1.0_dp - f(k))
      lam(k) = sqrt(3.0_dp*(1.0_dp - w_s(k))*(1.0_dp - w_s(k)*g_s(k)))
      gam(k) = 0.5_dp * w_s(k) * (1.0_dp + 3.0_dp*g_s(k)*(1.0_dp - w_s(k))*mu_z**2)/(1.0_dp - lam(k)**2*mu_z**2)
      alp(k) = 0.75_dp * w_s(k) * mu_z * (1.0_dp + g_s(k)*(1.0_dp - w_s(k)))/(1.0_dp - lam(k)**2*mu_z**2)
      u(k) = (3.0_dp/2.0_dp) * ((1.0_dp - w_s(k)*g_s(k))/lam(k))

      lamtau = min(lam(k)*tau_Ve_s(k),99.0_dp)
      e_lamtau = exp(-lamtau)

      N(k) = (u(k) + 1.0_dp)**2 * 1.0_dp/e_lamtau - (u(k) - 1.0_dp)**2  * e_lamtau

      R_b(k) = (u(k) + 1.0_dp)*(u(k) - 1.0_dp)*(1.0_dp/e_lamtau - e_lamtau)/N(k)
      T_b(k) = 4.0_dp * u(k)/N(k)

      arg = min(tau_Ve_s(k)/mu_z,99.0_dp)
      Tf(k) = exp(-arg)

      apg = alp(k) + gam(k)
      amg = alp(k) - gam(k)

      R(k) = amg*(T_b(k)*Tf(k) - 1.0_dp) + apg*R_b(k)

      T(k) = apg*T_b(k) + (amg*R_b(k) - (apg - 1.0_dp))*Tf(k)

      R(k) = max(R(k), 0.0_dp)
      T(k) = max(T(k), 0.0_dp)
      R_b(k) = max(R_b(k), 0.0_dp)
      T_b(k) = max(T_b(k), 0.0_dp)

    end do

    !! Calculate downward flux
    do k = 1, nlay
      sw_down(k) = Tf(k) + ((T(k) - Tf(k)) +  &
      & Tf(k)*R(k+1)*R_b(k))/(1.0_dp - R_b(k)*R_b(k+1))
    end do
    sw_down(nlev) = Tf(nlev)

    !! Calculate upward flux
    do k = 1, nlay
      sw_up(k) = (Tf(k)*R(k+1) + (T(k) - Tf(k))*R_b(k+1))/(1.0_dp - R_b(k)*R_b(k+1))
    end do
    sw_up(nlev) = sw_down(nlev) * w_surf

    !! Scale with the incident flux
    sw_down(:) = sw_down(:) * mu_z * Finc
    sw_up(:) = sw_up(:) * mu_z * Finc

  end subroutine bg_direct_beam_updown_adding

  subroutine bg_diffuse_updown_linear(surf, nlay, nlev, be, be_int, tau_IRe, lw_a_surf, lw_up, lw_down)
    implicit none

    !! Input variables
    logical, intent(in) :: surf
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_IRe
    real(dp), intent(in) :: be_int, lw_a_surf

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: k, m
    real(dp), dimension(nlay) :: dtau, edel
    real(dp) :: del, e0i, e1i, e1i_del
    real(dp), dimension(nlay) :: Am, Bm, Gp, Bp
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    !! Calculate dtau in each layer
    do k = 1, nlay
      dtau(k) = tau_IRe(k+1) - tau_IRe(k)
    end do

    ! Zero the total flux arrays
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do m = 1, nmu

      !! Prepare loop
      do k = 1, nlay
        ! Olson & Kunasz (1987) linear interpolant parameters
        del = dtau(k)/uarr(m)
        edel(k) = exp(-del)
        e0i = 1.0_dp - edel(k)
        e1i = del - e0i

        e1i_del = e1i/del ! The equivalent to the linear in tau term

        if (dtau(k) < 1.0e-6_dp) then
          ! If we are in very low optical depth regime, then use an isothermal approximation
          Am(k) = (0.5_dp*(be(k+1) + be(k)) * e0i)/be(k)
          Bm(k) = 0.0_dp
          Gp(k) = 0.0_dp
          Bp(k) = Am(k)
        else
          Am(k) = e0i - e1i_del ! Am(k) = Gp(k), just indexed differently
          Bm(k) = e1i_del ! Bm(k) = Bp(k), just indexed differently
          Gp(k) = Am(k)
          Bp(k) = Bm(k)
        end if
      end do

      !! Begin two-stream loops
      !! Perform downward loop first
      ! Top boundary condition - 0 flux downward from top boundary
      lw_down_g(1) = 0.0_dp
      do k = 1, nlay
        lw_down_g(k+1) = lw_down_g(k)*edel(k) + Am(k)*be(k) + Bm(k)*be(k+1) ! TS intensity
      end do

      !! Perform upward loop

      if (surf .eqv. .True.) then
        ! Surface boundary condition given by surface temperature + reflected longwave radiaiton
        lw_up_g(nlev) = lw_down_g(nlev)*lw_a_surf + be(nlev)
      else
        ! Lower boundary condition - internal heat definition Fint = F_down - F_up
        ! here the lw_a_surf is assumed to be = 1 as per the definition
        ! here we use the same condition but use intensity units to be consistent
        lw_up_g(nlev) = lw_down_g(nlev) + be_int
      end if

      do k = nlay, 1, -1
        lw_up_g(k) = lw_up_g(k+1)*edel(k) + Bp(k)*be(k) + Gp(k)*be(k+1) ! TS intensity
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      lw_down(:) = lw_down(:) + lw_down_g(:) * wuarr(m)
      lw_up(:) = lw_up(:) + lw_up_g(:) * wuarr(m)

    end do

    !! The flux is the intensity * 2pi
    lw_down(:) = twopi * lw_down(:)
    lw_up(:) = twopi * lw_up(:)

  end subroutine bg_diffuse_updown_linear

  subroutine calc_cff(surf, nlay, be, pe, tau, cff, scff )
    implicit none
    !! Input variables
    logical, intent(in) :: surf
    integer, intent(in) :: nlay
    real(dp), dimension(nlay+1), intent(in) :: be, pe, tau

    !! Output variables
    real(dp), dimension(nlay), intent(out) :: cff
    real(dp), intent(out) :: scff

    !! Work variables and arrays
    integer :: i
    real(dp) :: D_factor

    ! Contribution function to the planetary outgoing radiation
    D_factor = 1.66_dp ! constant diffusivity factor
    do i = 1, nlay
      !D_factor = 1.5_dp + 0.5_dp/(1_dp + 4_dp*(tau_struc_edges(k)+tau_struc_edges(k+1))/2_dp &
      !         + 10_dp*((tau_struc_edges(k)+tau_struc_edges(k+1))/2_dp)**2) ! diffusivity approximation (Thomas & Stammes 1999, section 11.2.5)
      cff(i) = ( 2.d0*pi*be(i)/D_factor ) * &
      & abs(exp(-D_factor*tau(i+1)) - exp(-D_factor*tau(i)) ) / &
      & abs( log10(pe(i+1))-log10(pe(i)) )
    end do
    !scff = 2_dp*pi*flux_up(nlay1)*exp(-D_factor*tau_struc_edges(nlay1))/D_factor
    scff = be(nlay+1)*exp(-D_factor*tau(nlay+1))

  end subroutine calc_cff

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / log10(x2/x1)

    yval = 10.0_dp**((ly1 * log10(x2/xval) + ly2 * log10(xval/x1)) * norm)

  end subroutine linear_log_interp

  subroutine bezier_interp(xi, yi, ni, x, y)
    implicit none

    integer, intent(in) :: ni
    real(dp), dimension(ni), intent(in) :: xi, yi
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y

    real(dp) :: xc, dx, dx1, dy, dy1, w, yc, t, wlim, wlim1

    !xc = (xi(1) + xi(2))/2.0_dp ! Control point (no needed here, implicitly included)
    dx = xi(2) - xi(1)
    dx1 = xi(3) - xi(2)
    dy = yi(2) - yi(1)
    dy1 = yi(3) - yi(2)

    if (x > xi(1) .and. x < xi(2)) then
      ! left hand side interpolation
      !print*,'left'
      w = dx1/(dx + dx1)
      !wlim = 1.0_dp + 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      !wlim1 = 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      !if (w < wlim .or. w > wlim1) then
      !  w = 1.0_dp
      !end if
      yc = yi(2) - dx/2.0_dp * (w*dy/dx + (1.0_dp - w)*dy1/dx1)
      t = (x - xi(1))/dx
      y = (1.0_dp - t)**2 * yi(1) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(2)
    else ! (x > xi(2) and x < xi(3)) then
      ! right hand side interpolation
      !print*,'right'
      w = dx/(dx + dx1)
      !wlim = 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      !wlim1 = 1.0_dp + 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      !if (w < wlim .or. w > wlim1) then
      !  w = 1.0_dp
      !end if
      yc = yi(2) + dx1/2.0_dp * (w*dy1/dx1 + (1.0_dp - w)*dy/dx)
      t = (x - xi(2))/(dx1)
      y = (1.0_dp - t)**2 * yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

  end subroutine bezier_interp

  subroutine BB_integrate(nlay1, Tint, wn_edges, bint)
    implicit none

    integer, intent(in) :: nlay1
    real(kind=dp), dimension(nlay1), intent(in) :: Tint
    real(kind=dp), dimension(2), intent(in) :: wn_edges
    real(kind=dp), dimension(nlay1), intent(out) :: bint

    integer :: i, w, j, intitera

    real(kind=dp), dimension(2) :: iB

    real(kind=dp) :: c1, x, x2, x3, itera, summ, n2, dn
    real(kind=dp) :: h, c, k, c2

    !! Code for integrating the blackbody function between two wavenumbers
    !! This is a method that uses a sum convergence.
    !! Taken from: spectralcalc.com/blackbody/inband_radiance.html

    h = 6.626075540D-34 
    c = 2.99792458D+08
    k = 1.38065812D-23
    c2 = c**2

    do i = 1, nlay1

      if (Tint(i) < 1e-6_dp) then
        bint(i) = 0.0_dp
        cycle
      end if

      do w = 1, 2

        c1 = (h * c) / k
        x = c1 * 100.0_dp * wn_edges(w)/ Tint(i)
        x2 = x**2
        x3 = x**3

        itera = 2.0_dp + 20.0_dp/x
        if (itera < 512) then
          itera = 512
        end if
        intitera = int(itera)

        summ = 0.0_dp
        do j = 1, intitera + 1
          dn = 1.0_dp/real(j,kind=dp)
          summ = summ + exp(-min(real(j,kind=dp)*x,300.0_dp))*&
          & (x3 + (3.0_dp * x2 + 6.0_dp*(x+dn)*dn)*dn)*dn
        end do

        n2 = 2.0_dp * h * c2
        iB(w) = n2 * (Tint(i)/c1)**4 * summ
      end do
      bint(i) = max(iB(2) - iB(1),0.0_dp)
    end do

  end subroutine BB_integrate

end module ts_short_char_bg_mod
