module rate_constant_laws

!-----------------------------------------------------------------------------
! Rate Constant Temperature (rct) and Interaction Correction (rcic) laws subroutines
!-----------------------------------------------------------------------------

  use constants

  implicit none

  type :: int_law_pars
    integer                              :: id   ! law id
    real(dp), dimension(n_max_rcic_pars) :: pars ! parameter list
  end type

contains

!-----------------------------------------------------------------------------
! Rate Constant Temperature (rct) laws
!-----------------------------------------------------------------------------

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


!-----------------------------------------------------------------------------
! Rate Constant Interaction Correction (rcic) laws
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  real(dp) function rcic_law(rcic_info, V_A, V_B)
!-----------------------------------------------------------------------------
    type(int_law_pars), intent(in) :: rcic_info
    real(dp), intent(in) :: V_A, V_B

    select case (rcic_info%id)
      case (rcic_linear_id)
        rcic_law = rcic_linear(V_A, V_B, rcic_info%pars)
      case default
        stop " rcic_law function error: Should never occur"
    end select

  end function rcic_law

!-----------------------------------------------------------------------------
  real(dp) function rcic_linear(V_A, V_B, parameters)
!-----------------------------------------------------------------------------
    real(dp), intent(in) :: V_A, V_B
    real(dp), dimension(:), intent(in) :: parameters
    real(dp) :: a, b

    a = parameters(1)
    b = parameters(2)

    rcic_linear = a*V_A + b*V_B

  end function rcic_linear



end module rate_constant_laws

