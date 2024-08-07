module constants

  use, intrinsic :: iso_fortran_env

  implicit none
  public

  integer, parameter :: dp = REAL64
  integer, parameter :: i8 = int8

  ! Physical and mathematical constants
  real(dp), parameter :: big_bang       = 0.0_dp
  real(dp), parameter :: sqrt2          = sqrt(2.0_dp)
  real(dp), parameter :: isqrt2         = 1.0_dp/sqrt2
  real(dp), parameter :: sqrt3          = sqrt(3.0_dp)
  real(dp), parameter :: pi             = acos(-1.0_dp)
  real(dp), parameter :: kB             = 8.61733238496e-5_dp       ! eV / K
  real(dp), parameter :: hbar           = 0.6582119514467406e-15_dp ! eV * s

  ! Version number taken from a preprocessor variable
  ! Linux
  character(len=40), parameter :: version = __GIT_VERSION__
  ! Windows
  ! character(len=40), parameter :: version =  "windows v1"

  !-------- Site type info ---------!

  ! Number of site types
  integer, parameter :: n_max_lat_site_types = 13
  ! List of site types
  character(len=10), dimension(n_max_lat_site_types), parameter :: &
            lat_site_names = [character(10):: 'terrace', 'step', 'corner', 'tc1', 'tc2', 'ts2', 'ts1',&
                                                         'stepA','cornerA','tc1A','tc2A','ts2A','ts1A']
  character(len=10), parameter :: same_lst_mark = '*'


  ! Site type ids
  integer, parameter :: terrace_site = 1
  integer, parameter :: step_site    = 2
  integer, parameter :: corner_site  = 3
  integer, parameter :: tc1_site     = 4
  integer, parameter :: tc2_site     = 5
  integer, parameter :: ts2_site     = 6
  integer, parameter :: ts1_site     = 7
  integer, parameter :: stepA_site   = 8
  integer, parameter :: cornerA_site = 9
  integer, parameter :: tc1A_site    = 10
  integer, parameter :: tc2A_site    = 11
  integer, parameter :: ts2A_site    = 12
  integer, parameter :: ts1A_site    = 13

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
  integer, parameter :: n_reaction_types = 5
  ! List of reaction types
  character(len=20), dimension(n_reaction_types), parameter :: &
            reaction_names = [character(20)::'hopping','desorption','dissociation','association','bimolecular']
  ! Reaction ids
  integer, parameter :: hopping_id      = 1
  integer, parameter :: desorption_id   = 2
  integer, parameter :: dissociation_id = 3
  integer, parameter :: association_id  = 4
  integer, parameter :: bimolecular_id  = 5

  !-------- Hopping Rate Constants' Interaction Correction Laws ---------!

  ! List of interaction laws
  character(len=10), dimension(1), parameter :: &
          rcic_law_names = [character(10)::'linear']
  ! law ids and number of parameters
  integer, parameter :: rcic_linear_id = 1
  integer, parameter :: rcic_linear_npars = 2
  ! Maximum number of parameters
  integer, parameter :: n_max_rcic_pars = 2


  !-------- Check mark
  character(len=2), parameter :: check_mark = '->'

  !-------- Section end keyword
  character(len=10), parameter :: section_end = ''


  !-------- Temperature dependences info ---------!

  ! Number of temperature dependence equations
  integer, parameter :: n_rct_laws = 2
  ! List of laws
  character(len=20), dimension(n_rct_laws), parameter :: &
            rct_law_names = &
                [character(20)::'Arrhenius','extArrhenius']
  ! Law ids
  integer, parameter :: Arrhenius_id    = 1
  integer, parameter :: extArrhenius_id = 2
  ! Maximum number of temperature law parameters
  integer, parameter :: n_max_rct_pars = 3


  ! file units
  integer, parameter :: inp_unit     = 1
  integer, parameter :: outcfg_unit  = 10
  integer, parameter :: outeng_unit  = 11
  integer, parameter :: outhst_unit  = 12
  integer, parameter :: outcnt_unit  = 13
  integer, parameter :: outrdf_unit  = 14
  integer, parameter :: outrav_unit  = 15

  ! Internal program constants
  integer, parameter :: randseed(13) = [7,5,3,11,9,1,17,2,9,6,4,5,8]
  integer, parameter :: max_string_length    = 1000
  real(dp), parameter :: tolerance           = 1.0e-9_dp
  real(dp), parameter :: exp_arg_too_big     = 100.0_dp

  ! Conversion constants to program units
  !
  ! Program basic units
  !           Temperature   : Kelvin
  !           Time          : s
  !           Energy        : eV
  real(dp), parameter :: Kelvin2eV        = kB


  ! Global variables for debugging
  !   1 prints out rct and rcic laws parameters
  logical :: debug(10)


end module constants
