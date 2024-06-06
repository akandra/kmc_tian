module temperature_laws
!-----------------------------------------------------------------------------
!             Temperature dependence law subroutines
!-----------------------------------------------------------------------------

  use constants

  implicit none

contains

!-----------------------------------------------------------------------------
  real(dp) function rct_law(rct_law_id, temperature, parameters)
!-----------------------------------------------------------------------------
    integer,  intent(in) :: rct_law_id
    real(dp), intent(in) :: temperature
    real(dp), dimension(:), intent(in) :: parameters

    select case (rct_law_id)
      case (Arrhenius_id)
        rct_law = arrhenius(temperature, parameters)
      case (extArrhenius_id)
        rct_law = extArrhenius(temperature, parameters)
      case default
        stop " rct_law function error: Should never occur"
    end select

  end function rct_law

!-----------------------------------------------------------------------------
  real(dp) function arrhenius(temperature, parameters)
!-----------------------------------------------------------------------------
    real(dp), intent(in) :: temperature
    real(dp), dimension(:), intent(in) :: parameters
    real(dp) :: prefactor, act_energy

    prefactor  = parameters(1)
    act_energy = parameters(2)

    arrhenius = prefactor*exp(-act_energy/(kB*temperature))

  end function arrhenius

!-----------------------------------------------------------------------------
  real(dp) function extArrhenius(temperature, parameters)
!-----------------------------------------------------------------------------
    real(dp), intent(in) :: temperature
    real(dp), dimension(:), intent(in) :: parameters
    real(dp) :: prefactor, act_energy, power, kT

      prefactor  = parameters(1)
      act_energy = parameters(2)
      power = parameters(3)
      kT = kB*temperature

      extArrhenius = prefactor*exp(-act_energy/kT)/(kT**power)

  end function extArrhenius

end module temperature_laws

