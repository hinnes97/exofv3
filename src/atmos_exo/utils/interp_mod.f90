module interp_mod
  use, intrinsic :: iso_fortran_env
  !use fv_control_mod, only: master
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  contains

    ! Perform linear interpolation in log10 space
    subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
      implicit none
  
      real(dp), intent(in) :: xval, y1, y2, x1, x2
      real(dp) :: lxval, ly1, ly2, lx1, lx2
      real(dp), intent(out) :: yval
      real(dp) :: norm
  
      lxval = log10(xval)
      lx1 = log10(x1); lx2 = log10(x2)
      ly1 = log10(y1); ly2 = log10(y2)
  
      norm = 1.0_dp / (lx2 - lx1)
  
      yval = 10.0_dp**((ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm)
  
    end subroutine linear_log_interp
  
    ! Perform linear interpolation in linear space
    subroutine linear_interp(xval, x1, x2, y1, y2, yval)
      implicit none
  
      real(dp), intent(in) :: xval, y1, y2, x1, x2
      real(dp), intent(out) :: yval
  
      yval = (y1 * (x2 - xval) + y2 * (xval - x1))/(x2 - x1)
  
    end subroutine linear_interp

    ! Perform bilinear interpolation in linear space
    subroutine bilinear_interp(xval, yval, x12, y12, matrix, P)
      implicit none

      real(dp), intent(in) :: xval, yval
      real(dp), dimension(2), intent(in) :: x12, y12
      real(dp), dimension(2,2), intent(in) :: matrix
      real(dp), intent(out) :: P
      real(dp) :: yval1, yval2

      yval1 = (matrix(1,1) * (x12(2) - xval) + matrix(2,1) * (xval - x12(1)))/(x12(2) - x12(1))
      yval2 = (matrix(1,2) * (x12(2) - xval) + matrix(2,2) * (xval - x12(1)))/(x12(2) - x12(1))

      P = ( 1.0_dp/((x12(2)-x12(1))*(y12(2)-y12(1))) ) * ( matrix(1,1)*(x12(2)-xval)*(y12(2)-yval) + &
                                                           matrix(2,1)*(xval-x12(1))*(y12(2)-yval) + &
                                                           matrix(1,2)*(x12(2)-xval)*(yval-y12(1)) + &
                                                           matrix(2,2)*(xval-x12(1))*(yval-y12(1)) )
    end subroutine bilinear_interp

    subroutine get_pk_edge(npz, is, ie, js, je, peln, kap_loc, pk_edge, kap_edge)
      integer, intent(in) :: npz, is, ie, js, je
      real, dimension(is:ie,1:npz+1,js:je), intent(in) :: peln
      real, dimension(is:ie,js:je,1:npz), intent(in) :: kap_loc
      real, dimension(is:ie,js:je,1:npz+1), intent(out) :: pk_edge
      real, dimension(is:ie,js:je,1:npz+1), optional, intent(out) :: kap_edge
      
      real :: pav(is:ie,js:je,1:npz)
      real :: kap_e(is:ie,js:je,1:npz+1)
      integer :: i,j,k

      call get_kappa_edge(npz, is,ie, js,je, peln,kap_loc, kap_e)
      if (present(kap_edge)) kap_edge(is:ie,js:je,1:npz+1) = kap_e(is:ie,js:je,1:npz+1)
      do k=1,npz+1
         do j=js,je
            do i=is,ie
               pk_edge(i,j,k) = exp(kap_e(i,j,k)*peln(i,k,j))
            enddo
         enddo
      enddo
      
    end subroutine get_pk_edge


    subroutine get_kappa_edge(npz,is,ie,js,je,peln,kap_loc,kap_edge)
      integer, intent(in) :: npz, is, ie, js, je
      real, dimension(is:ie,1:npz+1,js:je), intent(in) :: peln
      real, dimension(is:ie,js:je,1:npz), intent(in) :: kap_loc
      real, dimension(is:ie,js:je,1:npz+1), intent(out) :: kap_edge

      real :: pav(is:ie,js:je,1:npz)
      integer :: i,j,k
! Edge values
      do j=js,je
         do i=is,ie
            kap_edge(i,j,1) = kap_loc(i,j,1)
            kap_edge(i,j,npz+1) = kap_loc(i,j,npz)
         enddo
      enddo
   

! pav
      do k=1,npz
         do j=js,je
            do i=is,ie
               pav(i,j,k) = 0.5*(peln(i,k,j) + peln(i,k+1,j))
            enddo
         enddo
      enddo

! Middle values
      do k=2,npz
         do j=js,je
            do i=is,ie
               kap_edge(i,j,k) = kap_loc(i,j,k-1) + &
                    (peln(i,k,j) - pav(i,j,k-1))/(pav(i,j,k) - pav(i,j,k-1))*(kap_loc(i,j,k) - kap_loc(i,j,k-1))
               !pk_edge(i,j,k) = exp(kap_loc(i,j,k)*peln(i,k,j))
            enddo
         enddo
      enddo

    end subroutine get_kappa_edge
    
    subroutine get_kappa_edge_delp(npz,is,ie,js,je,ptop, delp,kap_loc,kap_edge, peln_out)
      integer, intent(in) :: npz, is, ie, js, je
      real, dimension(is:ie,js:je,npz), intent(in) :: delp
      real, intent(in) :: ptop
      real, dimension(is:ie,js:je,1:npz), intent(in) :: kap_loc
      real, dimension(is:ie,js:je,1:npz+1), intent(out) :: kap_edge
      real, dimension(is:ie,js:je,1:npz+1), intent(out), optional :: peln_out

      real :: pav(is:ie,js:je,1:npz)
      real :: peln(is:ie,js:je,1:npz+1)
      real :: pe(is:ie,js:je,1:npz+1)
      integer :: i,j,k
! Edge values
      do j=js,je
         do i=is,ie
            kap_edge(i,j,1) = kap_loc(i,j,1)
            kap_edge(i,j,npz+1) = kap_loc(i,j,npz)
         enddo
      enddo
   

! pav
      pe(is:ie,js:je,1) = ptop
      peln(is:ie,js:je,1) = log(ptop)

      do k=1,npz
         do j=js,je
            do i=is,ie
               pe(i,j,k+1) = pe(i,j,k) + delp(i,j,k)
               peln(i,j,k+1) = log(pe(i,j,k+1))
            enddo
         enddo
      enddo
         
      do k=1,npz
         do j=js,je
            do i=is,ie
               pav(i,j,k) = 0.5*(peln(i,j,k) + peln(i,j,k+1))
            enddo
         enddo
      enddo

! Middle values
      do k=2,npz
         do j=js,je
            do i=is,ie
               kap_edge(i,j,k) = kap_loc(i,j,k-1) + &
                    (peln(i,k,j) - pav(i,j,k-1))/(pav(i,j,k) - pav(i,j,k-1))*(kap_loc(i,j,k) - kap_loc(i,j,k-1))
!pk_edge(i,j,k) = exp(kap_loc(i,j,k)*peln(i,k,j))
            enddo
         enddo
      enddo

      if (present(peln_out)) peln_out(is:ie,js:je,1:npz+1) = peln
         
!      do k=1,npz
!         do j=js+1,je-1
!            do i=is+1,ie-1
!               if (kap_edge(i,j,k) .ne. kap_edge(i,j,k)) then
!                  write(*,*) i,j,k
!                  write(*,*) 'nan in kap edge'
!                  write(*,*) peln(i,j,k), pav(i,j,k), pav(i,j,k-1), kap_loc(i,j,k-1), kap_loc(i,j,k), kap_edge(i,j,k)
!                  write(*,*) delp(i,j,:)
!                  stop
!               endif
!            enddo
!         enddo
!      enddo
      

    end subroutine get_kappa_edge_delp

    subroutine get_pk_edge_delp(npz,is,ie,js,je,ptop, delp,kap_loc,pk_edge,kap_edge)
      integer, intent(in) :: npz, is, ie, js, je
      real, dimension(is:ie,js:je,npz), intent(in) :: delp
      real, intent(in) :: ptop
      real, dimension(is:ie,js:je,1:npz), intent(in) :: kap_loc
      real, dimension(is:ie,js:je,1:npz+1), intent(out) :: pk_edge
      real, dimension(is:ie,js:je,1:npz+1), optional, intent(out) :: kap_edge

      real :: kap_e(is:ie,js:je,1:npz+1)
      real :: peln(is:ie,js:je,1:npz+1)
      integer :: i,j,k

      
      call get_kappa_edge_delp(npz,is,ie,js,je,ptop,delp,kap_loc,kap_e, peln)

      if (present(kap_edge)) kap_edge(is:ie,js:je,1:npz+1) = kap_e(is:ie,js:je,1:npz+1)
      do k=1,npz+1
         do j=js,je
            do i=is,ie
               pk_edge(i,j,k) = exp(kap_e(i,j,k)*peln(i,k,j))
            enddo
         enddo
      enddo

      
    end subroutine get_pk_edge_delp
    
    subroutine get_pk_edge_2d(npz, is, ie, peln, kap_loc, pk_edge)
      integer, intent(in) :: npz,is,ie
      real, intent(in), dimension(is:ie,1:npz+1) :: peln
      real, intent(in), dimension(is:ie,1:npz)   :: kap_loc
      real, intent(out), dimension(is:ie,1:npz+1) :: pk_edge

      integer :: i,k, test(2)
      real, dimension(is:ie,1:npz+1) :: kap_edge
      real, dimension(is:ie,1:npz)   :: pav

      do k=1,npz
         do i=is,ie
            pav(i,k) = 0.5*(peln(i,k) + peln(i,k+1))
         enddo
      enddo
      do i=is,ie
         kap_edge(i,1) = kap_loc(i,1)
         kap_edge(i,npz+1) = kap_loc(i,npz)
      enddo
      do k=2,npz
         do i=is,ie
            kap_edge(i,k) = kap_loc(i,k-1) + &
                 (peln(i,k) - pav(i,k-1))/(pav(i,k) - pav(i,k-1))*(kap_loc(i,k) - kap_loc(i,k-1))
!pk_edge(i,j,k) = exp(kap_loc(i,j,k)*peln(i,k,j))
            pk_edge(i,k) = exp(kap_edge(i,k)*peln(i,k))
         enddo
      enddo
      
      
    end subroutine get_pk_edge_2d
    
end module interp_mod
