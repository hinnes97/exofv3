module rad_coupler_mod


   use constants_mod,        only: grav, rdgas, rvgas, RADIAN, kappa, radius, cp_air, rho_cp
   use mpp_domains_mod,      only: mpp_update_domains
   use mpp_mod,              only: FATAL, mpp_error, mpp_pe, stdlog, &
                                   mpp_npes, mpp_get_current_pelist, get_unit
   use time_manager_mod,     only: time_type, get_date, get_time
   use diag_manager_mod,     only: send_data, register_diag_field
   use fv_grid_utils_mod,    only: g_sum
   use fv_diagnostics_mod,   only: prt_maxmin
   use fv_timing_mod,        only: timing_on, timing_off
   use ts_short_char_mod,    only: ts_short_char
   use ts_short_char_mod_Bezier, only: ts_short_char_bezier
   use ts_short_char_bg_mod, only: ts_short_char_bg
   use ts_isothermal_mod,    only: ts_isothermal
   use ts_twostream_mod,    only: ts_twostream
   use read_opacities_mod,   only: read_opacities
   use read_planck_mod,      only: read_planck, Planck_table
   use interpolate_opacities_mod,   only: interpolate_opacities
   use spectral_partitioner_mod,    only: spectral_partition
   use netcdf_reg_mod,     only: netcdf_reg_RT
   implicit none

   logical :: rf_initialized = .false.
   private
   public :: rad_coupler, rad_coupler_init, rad_coupler_end

   CHARACTER(LEN=32) :: rad_scheme = 'None'
   real    :: I0 = 1000.0
   real    :: Tint = 200.0
   real    :: tau_IRe0 = 1.0
   real    :: tau_Ve0 = 1.0
   real    :: n_IR = 2.0
   real    :: n_V = 1.0
   real    :: f1 = 1.0
   real    :: fsw = 1.0
   real    :: flw = 1.0
   real    :: kappa_lw = 1.e-2
   real    :: kappa_sw = 1.e-4
   logical :: thick_atmosphere = .false.
   CHARACTER(LEN=32) :: star = 'Sun'

   integer :: n_bands = 1
   real, dimension(:), allocatable :: f_star, I_star

   namelist /rad_coupler_nml/ I0, Tint, rad_scheme, tau_IRe0, tau_Ve0, n_IR, n_V, f1, star, thick_atmosphere, &
                              kappa_sw, kappa_lw, flw, fsw
   namelist /n_bands_nml/ n_bands

contains


   subroutine rad_coupler_init(axes,Time)
      !use n_bands_mod, only: n_bands_init
      integer :: f_unit, unit, ios ! Namelist reading
      character(80) :: filename
      integer, intent(in), dimension(4) :: axes
      type(time_type), intent(in)       :: Time
      filename = "input.nml"


      f_unit = get_unit()
      open (f_unit,file=filename)
      rewind (f_unit)
      read (f_unit,rad_coupler_nml ,iostat=ios)
      unit = stdlog()
      write(unit, nml=rad_coupler_nml)

      close(f_unit)

      call netcdf_reg_RT(axes,Time,rad_scheme)

      select case(rad_scheme)
       case('ts_short_char_bg')

         f_unit = get_unit()
         open (f_unit,file=filename)
         rewind (f_unit)
         read (f_unit,n_bands_nml ,iostat=ios)
         unit = stdlog()
         write(unit, nml=n_bands_nml)

         close(f_unit)

         allocate(Planck_table(n_bands+1,2999), f_star(n_bands), I_star(n_bands))
         call spectral_partition(n_bands,star,f_star,I_star)
         ! Read tables before the timestepping
         call read_planck(n_bands)
         call read_opacities(n_bands)
      end select
         

   end subroutine rad_coupler_init

   subroutine rad_coupler_end()

     deallocate(Planck_table,f_star,I_star)

   end subroutine rad_coupler_end


   subroutine rad_coupler(surface_on, tidally_locked, is, ie, js, je, npz, ng, ts, pt, pl, pe, agrid, net_F, olr, opr_IR, &
                          opr_W1, opr_W2, opr_UV, opr_VIS1, opr_VIS2, net_Fs, surf_lw_down, surf_sw_down, cff, scff, direct_down)


      logical, intent(IN) :: surface_on, tidally_locked
      integer, intent(IN) :: npz, ng
      integer, intent(IN) :: is, ie, js, je

      real   , intent(IN) ::   pt(is:ie,js:je,npz)
      real   , INTENT(IN) ::   ts(is:ie,js:je)
      real   , intent(IN) ::   pl(is:ie,js:je,1:npz)
      real   , intent(IN) ::   pe(is:ie,js:je,1:npz+1)
      real   , intent(IN   ) :: agrid(is:ie,js:je, 2)

      real, intent(OUT)   :: net_F(is:ie,js:je,1:npz+1), cff(is:ie,js:je,1:npz), direct_down(is:ie,js:je,1:npz+1)
      real, intent(OUT)   :: olr(is:ie,js:je), opr_IR(is:ie,js:je), opr_W1(is:ie,js:je), opr_W2(is:ie,js:je),&
                             opr_UV(is:ie,js:je), opr_VIS1(is:ie,js:je), opr_VIS2(is:ie,js:je), scff(is:ie,js:je)
      real, intent(OUT)   :: net_Fs(is:ie,js:je)
      real, intent(OUT)   :: surf_lw_down(is:ie,js:je), surf_sw_down(is:ie,js:je)

      integer :: i,j,k,b
      real :: PI=4.D0*DATAN(1.D0)

      real, dimension(npz+1) :: tau_Ve, tau_IRe, net_F_1D, sw_down_1D, mu_z_arr
      real, dimension(npz) :: sw_a, sw_g, lw_a, lw_g
      real :: mu_z, AB, sw_a_surf, lw_a_surf, olr_1D, net_Fs_1D, I_1D
      real :: surf_lw_down_1D, surf_sw_down_1D
      real, dimension(n_bands) :: I_1D_b, opr_1D, net_Fs_bg_1D, scff_1D
      real, dimension(n_bands,npz+1) :: diffuse_up_1D, diffuse_down_1D, direct_up_1D, direct_down_1D, net_F_bg_1D
      real, dimension(n_bands,npz) :: cff_1D

      real,dimension(is:ie,js:je,1:npz) :: tau_Ve_3D,tau_IRe_3D
      real,dimension(is:ie,js:je) :: I_2D

      real, dimension(is:ie,js:je,n_bands,npz) :: aa, gg, kRoss
      real, dimension(is:ie,js:je,n_bands,npz+1) :: tau_bg

      !-----------------------------------------------------------------------------------------------
      !                                          RADIATION
      !-----------------------------------------------------------------------------------------------

      ! Choice of scheme (SELECT CASE ?): user input should be in a namelist
      ! Choices should include:
      !   - ts_short_char : would be the standard scheme, replacing radiation.f90
      !   - ts_disort_scatter
      !   - ts_Heng
      !   - ts_isothermal_2
      !   - ts_isothermal
      !   - ts_Lewis_scatter
      !   - ts_Mendonca
      !   - ts_Toon
      !   - ts_Toon_scatter
      !   - ts_short_char_bg:     band-grey (4, 5, or 6 bands) version of ts_short_char 
      !   - ts_disort_scatter_bg: band-grey (4, 5, or 6 bands) version of ts_disort_scatter with parameterised clouds
      !   - SOCRATES

      mu_z = 0.577
      AB = 0.0
      sw_a(1:npz) = 0.0
      sw_g(1:npz) = 0.0
      lw_a(1:npz) = 0.0
      lw_g(1:npz) = 0.0
      sw_a_surf = 0.0
      lw_a_surf = 0.0



      select case(rad_scheme)
       case('ts_short_char_bg')

         do i = is,ie
            do j = js,je

               do k = 1, 3
                 aa(i,j,k,1:npz) = lw_a
                 gg(i,j,k,1:npz) = lw_g
               end do
               do k = 4, 6
                 aa(i,j,k,1:npz) = sw_a
                 gg(i,j,k,1:npz) = sw_g
               end do

               call interpolate_opacities(npz,pt(i,j,1:npz),pl(i,j,1:npz),n_bands,kRoss(i,j,:,:))

               do k = 1, npz+1
                 tau_IRe(k) = f1*tau_IRe0 * (pe(i,j,k) / pe(i,j,npz+1)) + (1-f1)*tau_IRe0 * (pe(i,j,k) / pe(i,j,npz+1)) **n_IR
                 tau_Ve(k) = tau_Ve0 * (pe(i,j,k) / pe(i,j,npz+1)) ** n_V
               end do

               tau_bg(i,j,:,1) = 0.0

               if (n_bands .eq. 1) then
                  
                 do k = 1, npz
                   !kRoss(1,k) = tau_IRe0*10.0/pe(i,j,npz+1) ! to be exactly equivalent to the tauinf*p/ps approach
                   tau_bg(i,j,1,k+1) = tau_bg(i,j,1,k) + (kRoss(i,j,1,k) * abs(pe(i,j,k+1) - pe(i,j,k))) / grav
                   !kRoss(2,k) = tau_Ve0*10.0/pe(i,j,npz+1)
                   tau_bg(i,j,2,k+1) = tau_bg(i,j,2,k) + (kRoss(i,j,4,k) * abs(pe(i,j,k+1) - pe(i,j,k))) / grav
                 end do

                 if (tidally_locked) then
                    I_1D_b(1) = I0
                    mu_z = cos(agrid(i,j,1))*cos(agrid(i,j,2))
                 else
                   I_1D_b(1) = I0 * cos(agrid(i,j,2))
                 end if

               else if (n_bands .ge. 4) then

                 do k = 1, npz
                   !kRoss(i,j,:,k) = tau_IRe0*10.0/pe(i,j,npz+1)
                   tau_bg(i,j,:,k+1) = tau_bg(i,j,:,k) + (kRoss(i,j,:,k) * abs(pe(i,j,k+1) - pe(i,j,k))) / grav
                 end do

                 do b = 1, n_bands
                   if (tidally_locked) then
                      I_1D_b(b) = I0 * f_star(b)
                      mu_z = cos(agrid(i,j,1))*cos(agrid(i,j,2))
                   else
                     I_1D_b(b) = I0 * f_star(b) * cos(agrid(i,j,2))
                   end if
                 end do

               end if


               call ts_short_char_bg(surface_on, .true., Planck_table, n_bands, npz, npz+1, ts(i,j), pt(i,j,1:npz), pl(i,j,1:npz), &
                                     pe(i,j,1:npz+1), tau_Ve, tau_IRe, tau_bg(i,j,:,:), mu_z, I_1D_b, Tint, AB, sw_a, sw_g, lw_a, &
                                     lw_g, sw_a_surf, lw_a_surf, net_F_bg_1D, diffuse_up_1D, diffuse_down_1D, direct_up_1D, &
                                     direct_down_1D, opr_1D, net_Fs_bg_1D, cff_1D, scff_1D)

               if (n_bands .eq. 1) then
                 net_F(i,j,1:npz+1) = net_F_bg_1D(1,:)
                 net_Fs(i,j) = net_Fs_bg_1D(1)
                 olr(i,j) = opr_1D(1)
                 cff(i,j,1:npz) = cff_1D(1,1:npz) 
                 scff(i,j) = scff_1D(1)
                 direct_down(i,j,1:npz+1) = direct_down_1D(1,:)


               else if (n_bands .ge. 4) then
                 net_F(i,j,1:npz+1) = sum(net_F_bg_1D(:,:),dim=1)
                 net_Fs(i,j) = sum(net_Fs_bg_1D(:))
                 olr(i,j) = sum(opr_1D(:))
                 opr_IR(i,j) = opr_1D(1)
                 opr_W1(i,j) = opr_1D(2)
                 opr_W2(i,j) = opr_1D(3)
                 opr_UV(i,j) = opr_1D(4)
                 opr_VIS1(i,j) = opr_1D(5)
                 opr_VIS2(i,j) = opr_1D(6)
                 cff(i,j,1:npz) = sum(cff_1D(:,1:npz),dim=1)
                 scff(i,j) = sum(scff_1D(:))
                 direct_down(i,j,1:npz+1) = sum(direct_down_1D(:,:),dim=1)


               end if

            end do
         end do

      case('ts_short_char')
!$OMP parallel do default (none) firstprivate(  &
!$OMP                           tau_IRe,tau_Ve,mu_z,sw_a,sw_g,lw_a,lw_g,sw_a_surf,lw_a_surf,&
!$OMP                                net_F_1D, olr_1D, net_Fs_1D, surf_lw_down_1D, surf_sw_down_1D,&
!$OMP                                sw_down_1D) &
!$OMP     shared(js,je,is,ie,npz,tau_IRe0, tau_Ve0, n_V, n_IR, tidally_locked, I_1D, I0, surface_on, &
!$OMP             pe, f1, ts, pt, pl, Tint, AB, &
!$OMP              surf_lw_down, surf_sw_down, direct_down, &
!$OMP                                   agrid,olr, net_Fs, net_F)    
         do j = js,je
            do i = is,ie

               do k = 1, npz+1
                  tau_IRe(k) = f1*tau_IRe0 * (pe(i,j,k) / pe(i,j,npz+1)) &
                       + (1-f1)*tau_IRe0 * (pe(i,j,k) / pe(i,j,npz+1)) **n_IR
                  tau_Ve(k) = tau_Ve0 * (pe(i,j,k) / pe(i,j,npz+1)) ** n_V
               end do

               if (tidally_locked) then
                  I_1D = I0
                  mu_z = -cos(agrid(i,j,1))*cos(agrid(i,j,2))
                  if (I_1D < 0) then
                     I_1D = 0
                  end if
               else
                  I_1D = I0 * cos(agrid(i,j,2))
               end if

               call ts_short_char(surface_on, .true.,npz, npz+1, ts(i,j), pt(i,j,1:npz), pl(i,j,1:npz), & 
                    pe(i,j,1:npz+1), tau_Ve, tau_IRe, mu_z, I_1D, Tint, AB, sw_a, sw_g, lw_a, lw_g,     &
                    sw_a_surf, lw_a_surf, net_F_1D, olr_1D, net_Fs_1D,surf_lw_down_1D,surf_sw_down_1D,sw_down_1D)

               net_F(i,j,:) = net_F_1D
               net_Fs(i,j) = net_Fs_1D
               olr(i,j) = olr_1D
               surf_lw_down(i,j) = surf_lw_down_1D
               surf_sw_down(i,j) = surf_sw_down_1D
               direct_down(i,j,1:npz+1)  = sw_down_1D

            end do
         end do
      case('ts_short_char_bezier')
!$OMP parallel do default (none) firstprivate(  &
!$OMP                           tau_IRe,tau_Ve,mu_z,mu_z_arr,sw_a,sw_g,lw_a,lw_g,sw_a_surf,lw_a_surf,&
!$OMP                                net_F_1D, olr_1D, net_Fs_1D, surf_lw_down_1D, surf_sw_down_1D,&
!$OMP                                sw_down_1D) &
!$OMP     shared(js,je,is,ie,npz,tau_IRe0, tau_Ve0, n_V, n_IR, tidally_locked, I_1D, I0, surface_on, &
!$OMP             pe, f1, ts, pt, pl, Tint, AB, &
!$OMP              surf_lw_down, surf_sw_down, direct_down, &
!$OMP                                   agrid,olr, net_Fs, net_F, &
!$OMP                   kappa_lw, kappa_sw, fsw, flw, pi)    
         do j = js,je
            do i = is,ie
               tau_IRe(1) = kappa_lw*(flw + (1.-flw)*0.5*pe(i,j,1)/pe(i,j,npz+1))*pe(i,j,1)/grav
               tau_Ve(1)  = kappa_sw*(fsw + (1.-fsw)*0.5*pe(i,j,1)/pe(i,j,npz+1))*pe(i,j,1)/grav
               do k = 2, npz+1
!                  tau_IRe(k) = f1*tau_IRe0 * (pe(i,j,k) / pe(i,j,npz+1)) &
!                       + (1-f1)*tau_IRe0 * (pe(i,j,k) / pe(i,j,npz+1)) **n_IR
!                  tau_Ve(k) = tau_Ve0 * (pe(i,j,k) / pe(i,j,npz+1)) ** n_V

                  tau_IRe(k) = tau_IRe(k-1) + kappa_lw*(flw + (1.-flw)*(pl(i,j,k-1)/pe(i,j,npz+1)))*(pe(i,j,k) - pe(i,j,k-1))/grav
                  tau_Ve(k)  =  tau_Ve(k-1) + kappa_sw*(fsw + (1.-fsw)*(pl(i,j,k-1)/pe(i,j,npz+1)))*(pe(i,j,k) - pe(i,j,k-1))/grav
               end do

               if (tidally_locked) then
                  I_1D = I0
                  mu_z = -cos(agrid(i,j,1))*cos(agrid(i,j,2))
                  mu_z_arr = mu_z
                  if (I_1D < 0) then
                     I_1D = 0
                  end if
               else
                  ! Divide by pi here so that \int_{-pi/2}^{pi/2} 2*pi*a^2 cos^2(theta) d(theta) = pi*a^2
                  I_1D = I0/pi
                  mu_z_arr(:) =  cos(agrid(i,j,2))
               end if

               call ts_short_char_Bezier(.true., surface_on, npz, npz+1, ts(i,j), pt(i,j,1:npz), pl(i,j,1:npz), & 
                    pe(i,j,1:npz+1), tau_Ve, tau_IRe, mu_z_arr, I_1D, Tint, AB, sw_a, sw_g,     &
                    sw_a_surf, lw_a_surf, net_F_1D, olr_1D, net_Fs_1D,surf_lw_down_1D,surf_sw_down_1D,sw_down_1D)

               net_F(i,j,:) = net_F_1D
               net_Fs(i,j) = net_Fs_1D
               olr(i,j) = olr_1D
               surf_lw_down(i,j) = surf_lw_down_1D
               surf_sw_down(i,j) = surf_sw_down_1D
               direct_down(i,j,1:npz+1)  = sw_down_1D

            end do
         end do

       case('ts_isothermal')
         do i = is,ie
            do j = js,je

               do k = 1, npz+1
                  tau_IRe(k) = f1*tau_IRe0 * (pe(i,j,k) / pe(i,j,npz+1)) &
                          + (1-f1)*tau_IRe0 * (pe(i,j,k) / pe(i,j,npz+1)) **n_IR
                  tau_Ve(k) = tau_Ve0 * (pe(i,j,k) / pe(i,j,npz+1)) ** n_V
               end do

               if (tidally_locked) then
                  I_1D = I0
                  mu_z = -cos(agrid(i,j,1))*cos(agrid(i,j,2))
                  if (I_1D < 0) then
                     I_1D = 0
                  end if
               else
                  I_1D = I0 * cos(agrid(i,j,2))
               end if

               call ts_isothermal(surface_on, npz, npz+1, ts(i,j), pt(i,j,1:npz), tau_Ve, tau_IRe, mu_z, &
                  I_1D, Tint, AB, sw_a, sw_g, lw_a, lw_g, sw_a_surf, lw_a_surf, net_F_1D, olr_1D,        &
                  net_Fs_1D,surf_lw_down_1D,surf_sw_down_1D)

               net_F(i,j,1:npz+1) = net_F_1D
               net_Fs(i,j) = net_Fs_1D
               olr(i,j) = olr_1D
               surf_lw_down(i,j) = surf_lw_down_1D
               surf_sw_down(i,j) = surf_sw_down_1D

            end do
         end do

       case('ts_twostream')

         do i = is,ie
            do j = js,je

               do k = 1, npz+1
                  tau_IRe(k) = f1*tau_IRe0 * (pe(i,j,k) / pe(i,j,npz+1)) &
                          + (1-f1)*tau_IRe0 * (pe(i,j,k) / pe(i,j,npz+1)) **n_IR
                  tau_Ve(k) = tau_Ve0 * (pe(i,j,k) / pe(i,j,npz+1)) ** n_V
               end do


               if (tidally_locked) then
                  I_1D = I0 * cos(agrid(i,j,2)) * cos(agrid(i,j,1)+PI)
                  if (I_1D < 0) then
                     I_1D = 0
                  end if
               else
                  I_1D = I0 * cos(agrid(i,j,2))
               end if

               call ts_twostream(surface_on, thick_atmosphere, npz, npz+1, ts(i,j), pt(i,j,1:npz), tau_Ve, tau_IRe, mu_z, &
                  I_1D, Tint, AB, sw_a, sw_g, lw_a, lw_g, sw_a_surf, lw_a_surf, net_F_1D, olr_1D,       &
                  net_Fs_1D,surf_lw_down_1D,surf_sw_down_1D)

               net_F(i,j,1:npz+1) = net_F_1D
               net_Fs(i,j) = net_Fs_1D
               olr(i,j) = olr_1D
               surf_lw_down(i,j) = surf_lw_down_1D
               surf_sw_down(i,j) = surf_sw_down_1D

            end do
         end do

       case('None')

         net_F = 0.0
         net_Fs = 0.0
         olr = 0.0

      end select


      !-----------------------------------------------------------------------------------------------
      !                                      SEND DATA TO NETCDF
      !-----------------------------------------------------------------------------------------------

      ! Uses the send_data scheme with the template associated to output X:

      !if ( id_X > 0 ) then
      !  used = send_data ( id_X, X, Time)
      !endif

      !-----------------------------------------------------------------------------------------------
      !                                         DEALLOCATIONS
      !-----------------------------------------------------------------------------------------------

      ! Deallocate here any array allocated above.


   end subroutine rad_coupler


end module rad_coupler_mod

