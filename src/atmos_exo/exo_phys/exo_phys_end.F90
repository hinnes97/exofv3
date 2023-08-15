module exo_phys_end_mod

  !use radiation_mod,       only: radiation_end
  use rad_coupler_mod,    only: rad_coupler_end

  implicit none

  private
  public :: Exo_End

contains

  subroutine Exo_End

    call rad_coupler_end()

  end subroutine Exo_End

end module exo_phys_end_mod

