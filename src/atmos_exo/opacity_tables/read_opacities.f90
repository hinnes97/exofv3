module read_opacities_mod 
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = REAL64
  integer, parameter :: nbPT = 100 ! number of PT points on which the opacities are computed
  real(dp), dimension(nbPT,nbPT) :: matrixIR,matrixW1,matrixW2,matrixSW,matrixUV,matrixVIS,matrixVIS1,matrixVIS2

  contains

    subroutine read_opacities(n_bands)
      implicit none
      integer, intent(in) :: n_bands
      character(len=150) :: OpacityfileIR,OpacityfileW1,OpacityfileW2,OpacityfileSW,&
                            OpacityfileUV,OpacityfileVIS,OpacityfileVIS1,OpacityfileVIS2,&
                            path_to_opacities
      logical :: UseRosselandMean, UsePlanckMean, UseRossPlanckmean
      integer :: i
      
      !----------------------------------- Reading opacity tables -----------------------------------
      ! Interpolated tables of size 100x100.
      ! Regions SW and IR: POKAZATEL
      ! Windows W1 and W2: MT_CKD 3.0
      ! P range: [1e-8 , 1e3] bar = [1e-3 , 1e8] Pa. T range: [150 , 2900] K.

      ! Choose the Rosseland mean, the Planck mean, or the mean of the two.
      UseRosselandMean  = .False.
      UsePlanckMean     = .False.
      UseRossPlanckmean = .True.

      path_to_opacities = '../opacities'
      if (UseRosselandMean .eqv. .true.) then
        OpacityfileIR   = trim(path_to_opacities)//'/POKAZATEL/Rosseland_means/regionIR_Rosselandmean.txt'
        OpacityfileW1   = trim(path_to_opacities)//'/MTCKD/Rosseland_means/regionW1_Rosselandmean.txt'
        OpacityfileW2   = trim(path_to_opacities)//'/MTCKD/Rosseland_means/regionW2_Rosselandmean.txt'
        if (n_bands .eq. 4) then
          OpacityfileSW  = trim(path_to_opacities)//'/POKAZATEL/Rosseland_means/regionSW_Rosselandmean.txt'
        elseif (n_bands .eq. 5) then
          OpacityfileUV  = trim(path_to_opacities)//'/POKAZATEL/Rosseland_means/regionUV_Rosselandmean.txt'
          OpacityfileVIS = trim(path_to_opacities)//'/POKAZATEL/Rosseland_means/regionVIS_Rosselandmean.txt'
        elseif (n_bands .eq. 6) then
          OpacityfileUV   = trim(path_to_opacities)//'/POKAZATEL/Rosseland_means/regionUV_Rosselandmean.txt'
          OpacityfileVIS1 = trim(path_to_opacities)//'/POKAZATEL/Rosseland_means/regionVIS1_Rosselandmean.txt'
          OpacityfileVIS2 = trim(path_to_opacities)//'/POKAZATEL/Rosseland_means/regionVIS2_Rosselandmean.txt'
        end if
      elseif (UsePlanckmean .eqv. .true.) then
        OpacityfileIR   = trim(path_to_opacities)//'/POKAZATEL/Planck_means/regionIR_Planckmean.txt'
        OpacityfileW1   = trim(path_to_opacities)//'/MTCKD/Planck_means/regionW1_Planckmean.txt'
        OpacityfileW2   = trim(path_to_opacities)//'/MTCKD/Planck_means/regionW2_Planckmean.txt'
        if (n_bands .eq. 4) then
          OpacityfileSW  = trim(path_to_opacities)//'/POKAZATEL/Planck_means/regionSW_Planckmean.txt'
        elseif (n_bands .eq. 5) then
          OpacityfileUV  = trim(path_to_opacities)//'/POKAZATEL/Planck_means/regionUV_Planckmean.txt'
          OpacityfileVIS = trim(path_to_opacities)//'/POKAZATEL/Planck_means/regionVIS_Planckmean.txt'
        elseif (n_bands .eq. 6) then
          OpacityfileUV   = trim(path_to_opacities)//'/POKAZATEL/Planck_means/regionUV_Planckmean.txt'
          OpacityfileVIS1 = trim(path_to_opacities)//'/POKAZATEL/Planck_means/regionVIS1_Planckmean.txt'
          OpacityfileVIS2 = trim(path_to_opacities)//'/POKAZATEL/Planck_means/regionVIS2_Planckmean.txt'
        end if
      elseif (UseRossPlanckmean .eqv. .true.) then
        OpacityfileIR   = trim(path_to_opacities)//'/POKAZATEL/RossPlanck_means/regionIR_RossPlanckmean.txt'
        OpacityfileW1   = trim(path_to_opacities)//'/MTCKD/RossPlanck_means/regionW1_RossPlanckmean.txt'
        OpacityfileW2   = trim(path_to_opacities)//'/MTCKD/RossPlanck_means/regionW2_RossPlanckmean.txt'
        if (n_bands .eq. 4) then
          OpacityfileSW  = trim(path_to_opacities)//'/POKAZATEL/RossPlanck_means/regionSW_RossPlanckmean.txt'
        elseif (n_bands .eq. 5) then
          OpacityfileUV  = trim(path_to_opacities)//'/POKAZATEL/RossPlanck_means/regionUV_RossPlanckmean.txt'
          OpacityfileVIS = trim(path_to_opacities)//'/POKAZATEL/RossPlanck_means/regionVIS_RossPlanckmean.txt'
        elseif (n_bands .eq. 6) then
          OpacityfileUV   = trim(path_to_opacities)//'/POKAZATEL/RossPlanck_means/regionUV_RossPlanckmean.txt'
          OpacityfileVIS1 = trim(path_to_opacities)//'/POKAZATEL/RossPlanck_means/regionVIS1_RossPlanckmean.txt'
          OpacityfileVIS2 = trim(path_to_opacities)//'/POKAZATEL/RossPlanck_means/regionVIS2_RossPlanckmean.txt'
        end if
      end if

      open(11,file=OpacityfileIR)
      do i=1,nbPT
        read(11,*) matrixIR(i,:) ! matrixIR(i,j), i is the pressure index, j is the temperature index
      end do
      close(11)
      !print*, "matrixIR(55,:)", matrixIR(55,:)

      open(11,file=OpacityfileW1)
      do i=1,nbPT
        read(11,*) matrixW1(i,:)
      end do
      close(11)
      !print*, "matrixW1(55,:)", matrixW1(55,:)

      open(11,file=OpacityfileW2)
      do i=1,nbPT
        read(11,*) matrixW2(i,:)
      end do
      close(11)

      if (n_bands .eq. 4) then

        open(11,file=OpacityfileSW)
        do i=1,nbPT
          read(11,*) matrixSW(i,:)
        end do
        close(11)

      elseif (n_bands .eq. 5) then

        open(11,file=OpacityfileUV)
        do i=1,nbPT
          read(11,*) matrixUV(i,:)
        end do
        close(11)

        open(11,file=OpacityfileVIS)
        do i=1,nbPT
          read(11,*) matrixVIS(i,:)
        end do
        close(11)

      elseif (n_bands .eq. 6) then

        open(11,file=OpacityfileUV)
        do i=1,nbPT
          read(11,*) matrixUV(i,:)
        end do
        close(11)

        open(11,file=OpacityfileVIS1)
        do i=1,nbPT
          read(11,*) matrixVIS1(i,:)
        end do
        close(11)

        open(11,file=OpacityfileVIS2)
        do i=1,nbPT
          read(11,*) matrixVIS2(i,:)
        end do
        close(11)

      end if

    end subroutine read_opacities

end module read_opacities_mod
