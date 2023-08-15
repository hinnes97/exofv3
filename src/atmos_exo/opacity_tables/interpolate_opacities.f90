module interpolate_opacities_mod
  use, intrinsic :: iso_fortran_env
  use linspace_mod, only: linspace
  use interp_mod, only: bilinear_interp
  use read_opacities_mod, only: nbPT, matrixIR,matrixW1,matrixW2,matrixSW,matrixUV,matrixVIS,matrixVIS1,matrixVIS2
  implicit none

  integer, parameter :: dp = REAL64

  contains

    subroutine interpolate_opacities(nlay,Tmid,pmid,n_bands,band_opacity)
      implicit none
      integer :: i
      real(dp) :: pL(nbPT), TL(nbPT)
      real(dp) :: p_closest,T_closest, p_closest_lower, p_closest_upper, T_closest_lower, T_closest_upper, dy
      real(dp), dimension(2) :: p2, T2
      real(dp), dimension(2,2) :: matrixIR_22,matrixW1_22,matrixW2_22,matrixSW_22,&
                                  matrixUV_22,matrixVIS_22,matrixVIS1_22,matrixVIS2_22
      integer, intent(in) :: nlay, n_bands
      real(dp), dimension(nlay), intent(in) :: pmid,Tmid
      real(dp), dimension(n_bands,nlay), intent(out) :: band_opacity
      real(kind=dp) :: SG_tweak_lw, SG_tweak_sw

      !----------------------------------- Reading opacity tables -----------------------------------
      ! Interpolated tables of size 100x100.
      ! Regions SW and IR: POKAZATEL
      ! Windows W1 and W2: MT_CKD 3.0
      ! P range: [1e-8 , 1e3] bar = [1e-3 , 1e8] Pa. T range: [150 , 2900] K.

      ! P and T arrays corresponding to the tables
      call linspace(pL,-3.0_dp,8.0_dp,nbPT)
      pL = 10.0_dp**pL ! numpy.logspace
      call linspace(TL,150.0_dp,2900.0_dp,nbPT)

      ! Fill
      do i = 1,nlay
        T_closest = minloc(abs(TL-(TL/TL)*Tmid(i)), DIM=1) ! returns the index of the closest element to Tmid(i) in TL

        if (TL(T_closest) < Tmid(i)) then ! lower neighbor
          T_closest_lower = T_closest
          T_closest_upper = T_closest+1
        else if (TL(T_closest) > Tmid(i)) then ! upper neighbor
          T_closest_upper = T_closest
          T_closest_lower = T_closest-1
        else if (TL(T_closest) == Tmid(i)) then
        end if

        p_closest = minloc(abs(pL-(pL/pL)*pmid(i)), DIM=1)

        if (pL(p_closest) < pmid(i)) then ! lower neighbor
          p_closest_lower = p_closest
          p_closest_upper = p_closest+1
        else if (pL(p_closest) > pmid(i)) then ! upper neighbor
          p_closest_upper = p_closest
          p_closest_lower = p_closest-1
        else if (pL(p_closest) == pmid(i)) then
        end if

        T2            = (/ TL(T_closest_lower), TL(T_closest_upper) /)
        p2            = (/ pL(p_closest_lower), pL(p_closest_upper) /)

        matrixIR_22   = matrixIR( p_closest_lower:p_closest_upper,T_closest_lower:T_closest_upper )
        matrixW1_22   = matrixW1( p_closest_lower:p_closest_upper,T_closest_lower:T_closest_upper )
        matrixW2_22   = matrixW2( p_closest_lower:p_closest_upper,T_closest_lower:T_closest_upper )

        ! Interpolate tables further to avoid vertical jumps in opacity and in net flux
        call bilinear_interp(pmid(i),Tmid(i),p2,T2,matrixIR_22,band_opacity(1,i))
        call bilinear_interp(pmid(i),Tmid(i),p2,T2,matrixW1_22,band_opacity(2,i))
        call bilinear_interp(pmid(i),Tmid(i),p2,T2,matrixW2_22,band_opacity(3,i))

        ! The interpolation works only if Tmid(i) stays within the P/T range of the tables
        if ( (Tmid(i) .gt. maxval(TL)) .or. (Tmid(i) .lt. minval(TL))&
        .or. (pmid(i) .gt. maxval(pL)) .or. (pmid(i) .lt. minval(pL)) ) then
          band_opacity(1,i) = matrixIR(p_closest,T_closest)
          band_opacity(2,i) = matrixIR(p_closest,T_closest)
          band_opacity(3,i) = matrixIR(p_closest,T_closest)
        end if

        ! Tweaking factors
        SG_tweak_lw = 0d0
        !band_opacity(1,i) = SG_tweak_lw !1e-3 !band_opacity(1,i) !*0.0
        !band_opacity(2,i) = SG_tweak_lw !1e-3 !band_opacity(2,i)!*0.0 ! 10bar: 1, 260bar: 7d-2
        !band_opacity(3,i) = SG_tweak_lw !1e-3 !band_opacity(3,i)!*0.0 ! 10bar: 1, 260bar: 7d-3

        !band_opacity(1,i) = band_opacity(1,i)*1e-5
        !band_opacity(2,i) = band_opacity(2,i)*1e-5
        !band_opacity(3,i) = band_opacity(3,i)*1e-5

        if (n_bands .eq. 4) then
          matrixSW_22 = matrixSW( p_closest_lower:p_closest_upper,T_closest_lower:T_closest_upper )

          call bilinear_interp(pmid(i),Tmid(i),p2,T2,matrixSW_22,band_opacity(4,i))

          if ( (Tmid(i) .gt. maxval(TL)) .or. (Tmid(i) .lt. minval(TL))&
          .or. (pmid(i) .gt. maxval(pL)) .or. (pmid(i) .lt. minval(pL)) ) then
            band_opacity(4,i) = matrixSW(p_closest,T_closest)
          end if

          band_opacity(4,i) = band_opacity(4,i)*3d-2 ! 10bar: 3d-2 , 260bar: 1d-3

        elseif (n_bands .eq. 5) then
          matrixUV_22  = matrixUV ( p_closest_lower:p_closest_upper,T_closest_lower:T_closest_upper )
          matrixVIS_22 = matrixVIS( p_closest_lower:p_closest_upper,T_closest_lower:T_closest_upper )
           
          call bilinear_interp(pmid(i),Tmid(i),p2,T2,matrixUV_22,band_opacity(4,i))
          call bilinear_interp(pmid(i),Tmid(i),p2,T2,matrixVIS_22,band_opacity(5,i))

          if ( (Tmid(i) .gt. maxval(TL)) .or. (Tmid(i) .lt. minval(TL))&
          .or. (pmid(i) .gt. maxval(pL)) .or. (pmid(i) .lt. minval(pL)) ) then
            band_opacity(4,i) = matrixUV(p_closest,T_closest)
            band_opacity(5,i) = matrixVIS(p_closest,T_closest)
          end if

          !band_opacity(4,i) = band_opacity(4,i)*0.0
          !band_opacity(5,i) = band_opacity(5,i)*0.0
        elseif (n_bands .eq. 6) then
          matrixUV_22   = matrixUV  ( p_closest_lower:p_closest_upper,T_closest_lower:T_closest_upper )
          matrixVIS1_22 = matrixVIS1( p_closest_lower:p_closest_upper,T_closest_lower:T_closest_upper )
          matrixVIS2_22 = matrixVIS2( p_closest_lower:p_closest_upper,T_closest_lower:T_closest_upper )

          call bilinear_interp(pmid(i),Tmid(i),p2,T2,matrixUV_22,band_opacity(4,i))
          call bilinear_interp(pmid(i),Tmid(i),p2,T2,matrixVIS1_22,band_opacity(5,i))
          call bilinear_interp(pmid(i),Tmid(i),p2,T2,matrixVIS2_22,band_opacity(6,i))

          if ( (Tmid(i) .gt. maxval(TL)) .or. (Tmid(i) .lt. minval(TL))&
          .or. (pmid(i) .gt. maxval(pL)) .or. (pmid(i) .lt. minval(pL)) ) then
            band_opacity(4,i) = matrixUV(p_closest,T_closest)
            band_opacity(5,i) = matrixVIS1(p_closest,T_closest)
            band_opacity(6,i) = matrixVIS2(p_closest,T_closest)
          end if

          SG_tweak_sw = 0d0
          !band_opacity(4,i) = SG_tweak_sw !6e-4
          !band_opacity(5,i) = SG_tweak_sw !6e-4
          !band_opacity(6,i) = SG_tweak_sw !6e-4

          !band_opacity(4,i) = band_opacity(4,i)*1e-5
          !band_opacity(5,i) = band_opacity(5,i)*1e-5
          !band_opacity(6,i) = band_opacity(6,i)*1e-5

        end if

      end do
  
    end subroutine interpolate_opacities

end module interpolate_opacities_mod 
