module reaction_class

  use constants
  use control_parameters_class
  use mc_lat_class
  use energy_parameters_class
  use energy_mod
  use open_file
  use utilities
  use rates_hopping_class
  use rates_desorption_class

  implicit none

  private
  public    :: reaction_init, reaction_type


  type :: reaction_type

    type(   hopping_type) :: hopping
    type(desorption_type) :: desorption

    real(dp) :: beta                          ! inverse thermodynamic temperature
    integer  :: n_ads_total

  contains
    procedure :: construct
!    procedure :: print

  end type


contains
!-------------------------------------1-----------------------------------------
  function reaction_init(c_pars, lat, e_pars)
!------------------------------------------------------------------------------
    type(reaction_type) reaction_init

    type(control_parameters), intent(in)    :: c_pars ! Warning: check where c_pars gets redefined
    type(mc_lat)            , intent(in)    :: lat
    type(energy_parameters) , intent(in)    :: e_pars


    ! Create a rate structure
    reaction_init%hopping    =    hopping_init(c_pars, lat, e_pars)
    reaction_init%desorption = desorption_init(c_pars, lat, e_pars)
!    call reaction_init%hopping%print(c_pars)
!    call reaction_init%desorption%print(c_pars)

    reaction_init%beta = 1.0_dp/(kB*c_pars%temperature)
    reaction_init%n_ads_total = lat%n_ads_tot()

  end function reaction_init

!-----------------------------------------------------------------------------
  subroutine construct(this, lat, e_pars)
!-----------------------------------------------------------------------------
    class(reaction_type),     intent(inout) :: this
    class(mc_lat),            intent(inout) :: lat
    class(energy_parameters), intent(in)    :: e_pars

    integer :: i

    ! Hopping
    if (this%hopping%is_defined) then
      do i=1,this%n_ads_total
        call this%hopping%construct(i, lat, e_pars, this%beta)
      end do
    end if

    ! Desorption
    if (this%desorption%is_defined) then
      do i=1,this%n_ads_total
        call this%desorption%construct(i, lat, e_pars, this%beta)
      end do
    end if

  end subroutine construct


end module reaction_class
