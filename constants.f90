module constants

  use, intrinsic :: iso_fortran_env

  implicit none
  public

  integer, parameter :: dp = REAL64

  ! Physical and mathematical constants
  real(dp), parameter :: big_bang       = 0.0_dp ! in fortnights
  real(dp), parameter :: sqrt2          = sqrt(2.0_dp)
  real(dp), parameter :: isqrt2         = 1.0_dp/sqrt2
  real(dp), parameter :: sqrt3          = sqrt(3.0_dp)
  real(dp), parameter :: pi             = acos(-1.0_dp)
  real(dp), parameter :: kB             = 8.61733238496e-5_dp       ! eV / K
  real(dp), parameter :: hbar           = 0.6582119514467406e-15_dp ! eV * s

  !-------- Site type info ---------!

  ! Number of site types
  integer, parameter :: n_max_site_types = 3
  ! List of site types
  character(len=10), dimension(n_max_site_types), parameter :: &
            site_names = [character(10)::'terrace','step','corner']
  ! Site type ids
  integer, parameter :: terrace_site = 1
  integer, parameter :: step_site    = 2
  integer, parameter :: corner_site  = 3

  !-------- Adsorption sites info ---------!

  ! Maximum number of adsorption sites in the unit cell
  integer, parameter :: n_max_ads_sites = 6
  ! List of adsorption sites for the hexagonal lattice
  character(len=3), dimension(n_max_ads_sites), parameter :: &
          ads_site_names = [character(3)::'top','fcc','hcp','br1','br2','br3']
  ! Adsorption site ids
  integer, parameter :: top_id = 1
  integer, parameter :: fcc_id = 2
  integer, parameter :: hcp_id = 3
  integer, parameter :: br1_id = 4
  integer, parameter :: br2_id = 5
  integer, parameter :: br3_id = 6

  !-------- Interaction energy laws ---------!

  ! List of interaction laws
  character(len=10), dimension(2), parameter :: &
          int_law_names = [character(10)::'linear','sqrt']
  ! Interaction energy law ids
  integer, parameter :: linear_id = 1
  integer, parameter :: sqrt_id   = 2
  ! Maximum number of interaction parameters
  ! defined by the number of shells
  integer, parameter :: n_shells = 3

  !-------- Reactions info ---------!

  ! Number of reaction types
  integer, parameter :: n_reaction_types = 2
  ! List of reaction types
  character(len=10), dimension(n_reaction_types), parameter :: &
            reaction_names = [character(10)::'hopping','desorption']
  ! Reaction ids
  integer, parameter :: hopping_id    = 1
  integer, parameter :: desorption_id = 2

  !-------- Temperature dependences info ---------!

  ! Number of temperature dependence equations
  integer, parameter :: n_laws = 2
  ! List of laws
  character(len=20), dimension(n_laws), parameter :: &
            law_names = [character(20)::'Arrhenius',&
                         'extArrhenius']
  ! Law ids
  integer, parameter :: Arrhenius_id   = 1
  integer, parameter :: extArrhenius_id = 2


  ! file units
  integer, parameter :: inp_unit     = 1
  integer, parameter :: outcfg_unit  = 10
  integer, parameter :: outeng_unit  = 11
  integer, parameter :: outhst_unit  = 12

  ! Internal program constants
  integer, parameter :: randseed(13) = [7,5,3,11,9,1,17,2,9,6,4,5,8]
  integer, parameter :: max_string_length    = 1000
  real(dp), parameter :: tolerance           = 1.0e-9_dp

  ! Conversion constants to program units
  !
  ! Program basic units
  !           Temperature   : Kelvin
  !           Time          : s
  !           Energy        : eV
  real(dp), parameter :: Kelvin2eV        = kB


end module constants
