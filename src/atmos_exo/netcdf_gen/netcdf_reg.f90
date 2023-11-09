module netcdf_reg_mod 
  use diag_manager_mod, only: register_diag_field
  use time_manager_mod, only: time_type, operator(+), operator(-), operator(/=)

  implicit none

  real :: missing_value = -999.
  integer :: id_net_F, id_net_Fs, id_olr, id_opr_IR, id_opr_W1, id_opr_W2, id_opr_UV, id_opr_VIS1, id_opr_VIS2,&
             id_cff, id_scff, id_direct_down, id_t_dt_rad, id_t_dt_conv, id_t_dt_conv_moist, id_height, id_isr,&
             id_surf_sw_down, id_surf_lw_down, id_flux_t, id_t_dt_ding, id_t_dt_lsc, id_t_dt_rainout, &
             id_q_dt_ding, id_q_dt_lsc, id_q_dt_rainout, id_sw_dn, id_lw_dn, id_lw_up

  !integer, dimension(n_bands) :: id_net_F_bg

  contains

    subroutine register_diagnostic_field(axes,Time,mod_name,quantity,label,units,id_X, half)
      
      integer, intent(in), dimension(4) :: axes
      type(time_type), intent(in)       :: Time
      character(len=*), intent(in)     :: mod_name, units
      character(len=*), intent(in)    :: quantity, label
      !logical, intent(in)               :: 3D ! X(lat,lon,z) or X(lat,lon) ?
      integer, intent(out)               :: id_X
      logical, intent(in),optional                :: half

!if (3D .eqv. .true.) then
      if (present(half) .and. half) then
         id_X = register_diag_field(mod_name, quantity, (/axes(1), axes(2), axes(4)/), Time, label, units, &
              missing_value = missing_value)
      else
         id_X = register_diag_field ( mod_name, quantity, axes(1:3), Time, label, units, missing_value=missing_value )
      endif
      !else
      !  id_X = register_diag_field ( mod_name, quantity, axes(1:2), Time, label, units, missing_value=missing_value )
      !end if

    end subroutine register_diagnostic_field

    subroutine register_diagnostic_field_2D(axes,Time,mod_name,quantity,label,units,id_X)
      
      integer, intent(in), dimension(4) :: axes
      type(time_type), intent(in)       :: Time
      character(len=*), intent(in)     :: mod_name, units
      character(len=*), intent(in)    :: quantity, label
      integer, intent(out)               :: id_X

      id_X = register_diag_field ( mod_name, quantity, axes(1:2), Time, label, units, missing_value=missing_value )

    end subroutine register_diagnostic_field_2D

    subroutine netcdf_reg_RT(axes,Time,scheme)!,n_bands)
      
      integer, intent(in), dimension(4) :: axes
      type(time_type), intent(in)       :: Time
      CHARACTER(LEN=*), intent(in) :: scheme
      !integer, optional, intent(in) :: n_bands
      !integer :: b 
      !character (len=20) :: b_string

      select case(scheme)
       case('ts_short_char_bg')

         call register_diagnostic_field(axes,Time,"dynamics",'net_F',&
                                        &'net flux','watts/m2',id_net_F)
         call register_diagnostic_field_2D(axes,Time,"dynamics",'olr',&
                                        &'outgoing radiation','watts/m2',id_olr)
         call register_diagnostic_field_2D(axes,Time,"dynamics",'opr_IR',&
                                        &'outgoing radiation','watts/m2',id_opr_IR)
         call register_diagnostic_field_2D(axes,Time,"dynamics",'opr_W1',&
                                        &'outgoing radiation','watts/m2',id_opr_W1)
         call register_diagnostic_field_2D(axes,Time,"dynamics",'opr_W2',&
                                        &'outgoing radiation','watts/m2',id_opr_W2)
         call register_diagnostic_field_2D(axes,Time,"dynamics",'opr_UV',&
                                        &'outgoing radiation','watts/m2',id_opr_UV)
         call register_diagnostic_field_2D(axes,Time,"dynamics",'opr_VIS1',&
                                        &'outgoing radiation','watts/m2',id_opr_VIS1)
         call register_diagnostic_field_2D(axes,Time,"dynamics",'opr_VIS2',&
                                        &'outgoing radiation','watts/m2',id_opr_VIS2)
         call register_diagnostic_field(axes,Time,"dynamics",'cff',&
                                        &'flux contribution function','watts/m2/dex',id_cff)
         call register_diagnostic_field_2D(axes,Time,"dynamics",'scff',&
                                        &'surface flux contribution function','watts/m2',id_scff)
         call register_diagnostic_field(axes,Time,"dynamics",'direct_down',&
                                        &'downward direct flux','watts/m2',id_direct_down)
         call register_diagnostic_field(axes,Time,"dynamics",'t_dt_conv_moist',&
                                        &'temperature tendency due to moist convection','K/s',id_t_dt_conv_moist)

       case('ts_short_char')

         call register_diagnostic_field(axes,Time,"dynamics",'net_F',&
                                        &'net flux','watts/m2',id_net_F, .true.)
         call register_diagnostic_field_2D(axes,Time,"dynamics",'olr',&
                                        &'outgoing radiation','watts/m2',id_olr)
         call register_diagnostic_field(axes,Time,"dynamics",'direct_down',&
              &'downward direct flux','watts/m2',id_direct_down)
         

      case('ts_short_char_bezier')

         call register_diagnostic_field(axes,Time,"dynamics",'net_F',&
                                        &'net flux','watts/m2',id_net_F, .true.)
         call register_diagnostic_field_2D(axes,Time,"dynamics",'olr',&
                                        &'outgoing radiation','watts/m2',id_olr)
         call register_diagnostic_field(axes,Time,"dynamics", "sw_dn", &
         & "downward sw_flux", "w/m^2", id_sw_dn, .true.)
         call register_diagnostic_field(axes,Time,"dynamics", "lw_dn", &
         & "downward lw_flux", "w/m^2", id_lw_dn, .true.)
         call register_diagnostic_field(axes,Time,"dynamics", "lw_up", &
         & "upward lw_flux", "w/m^2", id_lw_up, .true.)
         
       case('ts_isothermal')

         call register_diagnostic_field(axes,Time,"dynamics",'net_F',&
                                        &'net flux','watts/m2',id_net_F)
         call register_diagnostic_field_2D(axes,Time,"dynamics",'olr',&
                                        &'outgoing radiation','watts/m2',id_olr)

       case('ts_twostream')

         call register_diagnostic_field(axes,Time,scheme,'net_F',&
                                        &'net flux','watts/m2',id_net_F)
         call register_diagnostic_field_2D(axes,Time,scheme,'olr',&
                                        &'outgoing radiation','watts/m2',id_olr)

       case('None')
         print*, 'No radiation scheme selected.'

       case default
         print*, 'Please select a radiation scheme, or enter "None".'
         stop

      end select

         call register_diagnostic_field(axes,Time,"dynamics",'t_dt_rad',&
                                        &'temperature tendency due to radiation','K/s',id_t_dt_rad)
         call register_diagnostic_field(axes,Time,"dynamics",'t_dt_conv',&
                                        &'temperature tendency due to dry convection','K/s',id_t_dt_conv)
         call register_diagnostic_field(axes,Time,'dynamics','height',&
                                        &'height','m',id_height)
         call register_diagnostic_field_2D(axes,Time,'dynamics','isr',&
                                        &'isr','m',id_isr)
         call register_diagnostic_field_2D(axes,Time,'dynamics','surf_sw_down',&
                                        &'surf_sw_down','m',id_surf_sw_down)
         call register_diagnostic_field_2D(axes,Time,'dynamics','surf_lw_down',&
                                        &'surf_lw_down','m',id_surf_lw_down)
         call register_diagnostic_field_2D(axes,Time,'dynamics','flux_t',&
                                        &'flux_t','m',id_flux_t)
    end subroutine netcdf_reg_RT

    subroutine netcdf_reg_ding(axes,Time)

      integer, intent(in), dimension(4) :: axes
      type(time_type), intent(in)       :: Time
      
      call register_diagnostic_field(axes, Time, "dynamics", 't_dt_ding', &
           'temperature tendency due to ding convection', 'K/s', id_t_dt_ding)
      call register_diagnostic_field(axes, Time, "dynamics", 't_dt_lsc', &
           'temperature tendency due to large scale condensation', 'K/s', id_t_dt_lsc)
      call register_diagnostic_field(axes, Time, "dynamics", 't_dt_rainout', &
           'temperature tendency due to rainout (re-evaporation)', 'K/s', id_t_dt_rainout)
      call register_diagnostic_field(axes, Time, "dynamics", 'q_dt_ding', &
           'vapour tendency due to ding convection', 'kg/kg/s', id_q_dt_ding)
      call register_diagnostic_field(axes, Time, "dynamics", 'q_dt_lsc', &
           'vapour tendency due to large scale condensation', 'kg/kg/s', id_q_dt_lsc)
      call register_diagnostic_field(axes, Time, "dynamics", 'q_dt_rainout', &
           'vapour tendency due to rainout (re-evaporation)', 'kg/kg/s', id_q_dt_rainout)

    end subroutine netcdf_reg_ding
end module netcdf_reg_mod
