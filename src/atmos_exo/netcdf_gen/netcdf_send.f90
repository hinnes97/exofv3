module netcdf_send_mod 
  use diag_manager_mod, only: send_data
  use time_manager_mod, only: time_type, operator(+), operator(-), operator(/=)

  implicit none

  contains

    subroutine send_data_to_netcdf(is,ie,js,je,npz,n,Time,id_X,X)
      implicit none
      integer, intent(in) :: is,ie,js,je,npz,n
      type(time_type), intent(in)       :: Time
      integer, intent(in)               :: id_X
      real, dimension(is:ie,js:je,n), intent(in) :: X
      logical :: used

      if ( id_X > 0 ) then
        if (n .eq. npz+1) then 
          used = send_data ( id_X, X(is:ie,js:je,2:n), Time)
        else if (n .eq. npz) then
          used = send_data ( id_X, X(is:ie,js:je,1:n), Time)
        end if
      end if
      
    end subroutine send_data_to_netcdf

    subroutine send_2D_data_to_netcdf(is,ie,js,je,Time,id_X,X_2D)
      implicit none
      integer, intent(in) :: is,ie,js,je
      type(time_type), intent(in)       :: Time
      integer, intent(in)               :: id_X
      real, dimension(is:ie,js:je), optional, intent(in) :: X_2D
      logical :: used

      if ( id_X > 0 ) then
        used = send_data ( id_X, X_2D(is:ie,js:je), Time)
      end if

    end subroutine send_2D_data_to_netcdf

end module netcdf_send_mod
