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

  ! Surface site types
  integer(int8), parameter :: terrace_site = 1
  integer(int8), parameter :: step_site    = 2
  integer(int8), parameter :: corner_site  = 3

  ! file units
  integer, parameter :: inp_unit     = 1
  integer, parameter :: outcfg_unit  = 10
  integer, parameter :: outeng_unit  = 11
  integer, parameter :: outhst_unit  = 12


  ! Internal program constants
  integer, parameter :: randseed(13)      = [7,5,3,11,9,1,17,2,9,6,4,5,8]
  integer, parameter :: max_string_length = 1000
  real(dp), parameter :: tolerance        = 1.0e-9_dp

  ! Defaults
  integer,   parameter  :: default_int           = 0
  real(dp),  parameter  :: default_real          = 0.0_dp
  character, parameter  :: default_string        = ""
  logical,   parameter  :: default_bool          = .true.


  ! Conversion constants to program units
  !
  ! Program basic units
  !           Temperature   : Kelvin
  !           Time          : s
  !           Energy        : eV
  real(dp), parameter :: Kelvin2eV        = kB


end module constants
