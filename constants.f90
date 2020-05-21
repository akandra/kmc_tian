module constants

  use, intrinsic :: iso_fortran_env

  implicit none

  integer, parameter :: dp = REAL64

  ! Physical and mathematical constants
  real(dp), parameter :: big_bang       = 0.0_dp ! in fortnights
  real(dp), parameter :: sqrt2          = sqrt(2.0_dp)
  real(dp), parameter :: isqrt2         = 1.0_dp/sqrt2
  real(dp), parameter :: sqrt3          = sqrt(3.0_dp)
  real(dp), parameter :: pi             = acos(-1.0_dp)
  real(dp), parameter :: kB             = 8.61733238496e-5_dp       ! eV / K
  real(dp), parameter :: hbar           = 0.6582119514467406e-15_dp ! eV * s

  ! Number of site types
  integer, parameter :: n_site_types = 3
  ! Surface site types
  character(len=10), dimension(n_site_types), parameter :: &
            site_names = ['terrace','step','corner']
  integer, parameter :: terrace_site = 1
  integer, parameter :: step_site    = 2
  integer, parameter :: corner_site  = 3

  ! Maximum number of adsorption sites in the unit cell
  integer, parameter :: n_ads_sites = 6
  ! adsorption sites for the hexagonal lattice
  character(len=3), dimension(n_ads_sites), parameter :: &
          ads_site_names = ['top','fcc','hcp','br1','br2','br3']
  integer, parameter :: top_id = 1
  integer, parameter :: fcc_id = 2
  integer, parameter :: hcp_id = 3
  integer, parameter :: br1_id = 4
  integer, parameter :: br2_id = 5
  integer, parameter :: br3_id = 6

  ! List of interaction laws
  character(len=10), dimension(2), parameter :: &
          int_law_names = ['linear','sqrt']
  ! Maximum number of interaction parameters
  integer, parameter :: n_shells = 3

  ! file units
  integer, parameter :: inp_unit     = 1
  integer, parameter :: outcfg_unit  = 10
  integer, parameter :: outeng_unit  = 11
  integer, parameter :: outhst_unit  = 12

  ! Internal program constants
  integer, parameter :: randseed(13)      = [7,5,3,11,9,1,17,2,9,6,4,5,8]
  integer, parameter :: max_string_length = 1000
  real(dp), parameter :: tolerance        = 1.0e-9_dp

  ! Conversion constants to program units
  !
  ! Program basic units
  !           Temperature   : Kelvin
  !           Time          : s
  !           Energy        : eV
  real(dp), parameter :: Kelvin2eV        = kB


end module constants
