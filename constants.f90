module constants

  use, intrinsic:: iso_fortran_env
  use version_info 
    
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
  ! character(len=40), parameter :: version = __GIT_VERSION__
  ! Windows
  !character(len=40), parameter :: version =  "windows v1"
  
  character(len=40), parameter :: version =  VERSION_STRING
  
  
    
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

  !-------- Neighbors info ---------!

  ! Maximum number of members in a neighbor shell
  integer, parameter :: max_n_neighbors = 10

  !-------- Interaction energy laws ---------!

  ! List of interaction laws
  character(len=10), dimension(2), parameter :: &
          int_law_names = [character(10)::'linear']
  ! Interaction energy law ids
  integer, parameter :: linear_id = 1

  !-------- Reactions info ---------!

  ! Number of reaction types
  integer, parameter :: n_reaction_types = 5
  ! List of reaction types
  character(len=20), dimension(n_reaction_types), parameter :: &
            reaction_names = [character(20)::'hopping', 'desorption','dissociation','association','bimolecular']
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
  ! Don't use '' as section_end
  ! it will confuse the code in its current state
  character(len=10), parameter :: section_end = 'end'


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


  ! Stopping trigger name
  character(len=20) :: essential_name = 'essential'

  ! Global variables for debugging
  !   1 prints out rct and rcic laws parameters
  logical :: debug(10)
 
!
! Lattice definition
!
! hexagonal: fcc with step density 1/6 and lower
!
    ! Adsorption sites on the unit hex cell
    !  T---------B1---------.           1 top       (T)
    !   \  .              .  .          2 fcc       (F)
    !    \     S         .    .         3 hcp       (H)
    !     \        .   .       .        4 bridge 1  (B1)
    !      B3        B2         .       5 bridge 2  (B2)
    !       \       .     .      .      6 bridge 3  (B3)
    !        \    .          H    .
    !         \ .                  .
    !          \. . . . . . . . . . .

    ! lattice vectors for hex lattice
  real(dp), parameter ::  hex_lat_vec_1(2) = [ cos(0.0_dp),    -sin(0.0_dp)   ]
  real(dp), parameter ::  hex_lat_vec_2(2) = [ cos(pi/3.0_dp), -sin(pi/3.0_dp)]

    ! Shell-wise number of neighbors for hex lattice with A- and B-type steps
  integer, parameter :: hex_n_nn_terrace(3) = [6,6, 6]
  integer, parameter :: hex_n_nn_step(3)    = [5,5, 7]
  integer, parameter :: hex_n_nn_corner(3)  = [4,4,10]
  integer, parameter :: hex_n_nn_tc1(3)     = [5,5, 5]
  integer, parameter :: hex_n_nn_ts1(3)     = [6,7, 5]
!  integer, parameter :: hex_n_nn_ts1(3)  = [6,7,5] ! the nn for the 3rd shell is changed from 4 to 5 on 2025-08-13

  ! terrace, ts2, or tc2

  ! NN list for the hexagonal 111-structure (row,col)
  !  tc1   tc2   ts2    ts1   s
  !  c     tc1   tc2    ts2   ts1
  !  11    12    13***  14**  15***
  !
  !     21    22**   23*   24*   25**
  !
  !        31*** 32*   33    34*   35***
  !
  !           41**  42*   43*   44**  45
  !
  !              51*** 52**  53*** 54    55

  integer, parameter :: hex_shell_list_terrace(3,6,2) = &
    reshape([ &
        
      ! Nearest-neigbour (1st) shell (d = 1))
      [ (/ 0, 1/), (/ 1, 0/), (/ 1,-1/), (/ 0,-1/), (/-1, 0/), (/-1, 1/) ], & 
       
      ! Next-Nearest-neigbour (2nd) shell  (d = sqrt(3))
      [ (/ 1, 1/), (/ 2,-1/), (/ 1,-2/), (/-1,-1/), (/-2, 1/), (/-1, 2/) ], &
        
      ! Next-Next-Nearest-neigbour (3rd) shell  (d = 2)
      [ (/ 0, 2/), (/ 2, 0/), (/ 2,-2/), (/ 0,-2/), (/-2, 0/), (/-2, 2/) ]  &
    ], shape=[3,6,2], order=[3,2,1])

    ! step

  ! NN list for the hexagonal structure (row,col)

  ! ts2   ts1    s      c     tc1
  !                   04***
  !
  !  11    12    13***  14**  15
  !
  !     21    22**   23*   24*   25***
  !
  !        31*** 32*   33    34**   35
  !
  !           41**  42*   43*   44*** 45
  !
  !              51*** 52**  53*** 54    55

  integer, parameter :: hex_shell_list_step(3,7,2) = &
    reshape([ &
        
      ! Nearest-neigbour (1st) shell)
      [ (/-1, 1/), & ! step -> corner (d=2/3)
        (/ 1, 0/), & ! step -> step   (d=1)
        (/-1, 0/), & ! step -> step   (d=1)
        (/ 1,-1/), & ! step -> ts1    (d=1)
        (/ 0,-1/), & ! step -> ts1    (d=1)
        (/ 0, 0/), & ! dummy entry
        (/ 0, 0/)  & ! dummy entry
       ], &

      ! Next-Nearest-neigbour (2nd) shell
      [ (/ 0, 1/), & ! step -> corner (d=sqrt(13)/3)
        (/-2, 1/), & ! step -> corner (d=sqrt(13)/3)
        (/ 2,-1/), & ! step -> ts1    (d=sqrt(3))
        (/-1,-1/), & ! step -> ts1    (d=sqrt(3))
        (/ 1,-2/), & ! step -> ts2    (d=sqrt(3))
        (/ 0, 0/), & ! dummy entry
        (/ 0, 0/)  & ! dummy entry
       ], &

      ! Next-Next-Nearest-neigbour (3rd) shell
      [ (/-1, 2/), & ! step -> tc1  (d=4/3)
        (/ 2, 0/), & ! step -> step   (d=2)
        (/-2, 0/), & ! step -> step   (d=2)
        (/ 2,-2/), & ! step -> ts2    (d=2)
        (/ 0,-2/), & ! step -> ts2    (d=2)
        (/-3, 1/), & ! step -> corner (d>2)
        (/ 1, 1/)  & ! step -> corner (d>2)
       ] &
    ], shape=[3,7,2], order=[3,2,1])


  ! corner
  ! NN list for the hexagonal 111-structure (row,col)

  ! ts1    s     c      tc1   tc2
  !  11    12    13***  14    15***
  !
  !     21    22    23*   24**  25***
  !
  !        31    32**  33    34*   35***
  !
  !           41*** 42*   43*   44**  45
  !
  !              51*** 52**  53*** 54    55

  integer, parameter :: hex_shell_list_corner(3,10,2) = &
    reshape([ &

      ! Nearest-neigbour (1st) shell
      [ (/ 1,-1/), & ! corner -> step   (d=2/3)
        (/ 0, 1/), & ! corner -> tc1    (d=2/3)
        (/ 1, 0/), & ! corner -> corner (d=1)
        (/-1, 0/), & ! corner -> corner (d=1)
        (/ 0, 0/), & ! dummy entry
        (/ 0, 0/), & ! dummy entry
        (/ 0, 0/), & ! dummy entry
        (/ 0, 0/), & ! dummy entry
        (/ 0, 0/), & ! dummy entry
        (/ 0, 0/)  & ! dummy entry
      ], &

      ! Next-Nearest-neigbour (2nd) shell
      [ (/-1, 1/), & ! corner -> tc1  (d=sqrt(3)-2/3)
        (/ 1, 1/), & ! corner -> tc1  (d=sqrt(3)-2/3)
        (/ 2,-1/), & ! corner -> step (d=sqrt(13)/3)
        (/ 0,-1/), & ! corner -> step (d=sqrt(13)/3)
        (/ 0, 0/), & ! dummy entry
        (/ 0, 0/), & ! dummy entry
        (/ 0, 0/), & ! dummy entry
        (/ 0, 0/), & ! dummy entry
        (/ 0, 0/), & ! dummy entry
        (/ 0, 0/)  & ! dummy entry
      ], &

      ! Next-Next-Nearest-neigbour (3rd) shell
      [ (/ 0, 2/), & ! corner -> tc2    (d=sqrt( (sqrt(3)-1/3)^2 + 1/4 ) )
        (/ 2, 0/), & ! corner -> corner (d=2)
        (/ 2,-2/), & ! corner -> ts1    (d=sqrt( (5/3)^2 + 1/4 ) )
        (/ 1,-2/), & ! corner -> ts2    (d=sqrt(3))
        (/-2, 0/), & ! corner -> corner (d=2)
        (/-1, 2/), & ! corner -> ts1    (d=sqrt( (5/3)^2 + 1/4 ) )
        (/-2, 2/), & ! corner -> tc2    (d= )                      ---------------------------------
        (/ 0,-2/), & ! corner -> ts1    (d= )                      ---------------------------------
        (/-1,-1/), & ! corner -> step   (d= )                      ---------------------------------
        (/ 3,-1/)  & ! corner -> step   (d= )                      ---------------------------------
      ] &
    ], shape=[3,10,2], order=[3,2,1])

    ! ts1 (terrace adjacent to step)

    ! NN list for the hexagonal structure (row,col)

    !  tc2   ts2   ts1    s     c
    !  11    12    13***  14**  15**
    !
    !     21    22**   23*   24*   25**
    !
    !        31*** 32*   33    34*   35***
    !
    !           41**  42*   43*   44**  45
    !
    !              51*** 52**  53*** 54    55

  integer, parameter :: hex_shell_list_ts1(3,7,2) = &
    reshape([ &

      ! Nearest-neigbour (1st) shell
        [ (/ 0, 1/), & ! ts1 -> step (d=1)
          (/ 1, 0/), & ! ts1 -> ts1  (d=1)
          (/ 1,-1/), & ! ts1 -> ts2  (d=1)
          (/ 0,-1/), & ! ts1 -> ts2  (d=1)
          (/-1, 0/), & ! ts1 -> ts1  (d=1)
          (/-1, 1/), & ! ts1 -> step (d=1)
          (/ 0, 0/)  & ! dummy entry
        ], &

        ! Next-Nearest-neigbour (2nd) shell
          [ (/ 2,-1/), & ! ts1 -> ts2 (d=sqrt(3))
            (/ 1,-2/), & ! ts1 -> tc2 (d=sqrt(3))
            (/-1,-1/), & ! ts1 -> ts2 (d=sqrt(3))
            (/ 1, 1/), & ! ts1 -> s   (d=sqrt(3))
            (/-2, 1/), & ! ts1 -> s   (d=sqrt(3))
            (/-2, 2/), & ! ts1 -> c   (d=sqrt( (5/3)^2 + 1/4 ) )
            (/-1, 2/)  & ! ts1 -> c   (d=sqrt( (5/3)^2 + 1/4 ) )
          ], &

        ! Next-Next-Nearest-neigbour (3rd) shell
          [ (/ 0, 2/), & ! ts1 -> c  (d=sqrt( (sqrt(3)-1/3)^2 + 1/4 ) )
            (/ 2, 0/), & ! ts1 -> ts1 (d=2)
            (/ 2,-2/), & ! ts1 -> tc2 (d=2)
            (/ 0,-2/), & ! ts1 -> tc2 (d=2)
            (/-2, 0/), & ! ts1 -> ts1 (d=2) -------------------------------------------------------- 
            (/ 0, 0/), & ! dummy entry
            (/ 0, 0/)  & ! dummy entry
          ] &
    ], shape=[3,7,2], order=[3,2,1])

    ! tc1 (terrace adjacent to corner)

    ! NN list for the hexagonal structure (row,col)

    !  s     c     tc1    tc2   ts2
    !  11    12    13***  14**  15***
    !
    !     21   22**   23*   24*   25**
    !
    !        31    32*   33    34*   35***
    !
    !           41*** 42**  43*   44**  45
    !
    !              51    52    53*** 54    55

    integer, parameter :: hex_shell_list_tc1(3,5,2) = &
    reshape([ &

        ! Nearest-neigbour (1st) shell
          [ (/ 0,-1/), & ! tc1 -> c   (d=2/3)
            (/ 1, 0/), & ! tc1 -> tc1 (d=1)
            (/-1, 0/), & ! tc1 -> tc1 (d=1)
            (/ 0, 1/), & ! tc1 -> tc2 (d=1)
            (/-1, 1/)  & ! tc1 -> tc2 (d=1)
          ], &

        ! Next-Nearest-neigbour (2nd) shell
          [ (/-1,-1/), & ! tc1 -> c   (d = sqrt(3)-2/3)
            (/ 1,-1/), & ! tc1 -> c   (d = sqrt(3)-2/3)
            (/ 1, 1/), & ! tc1 -> tc2 (d = sqrt(3))
            (/-2, 1/), & ! tc1 -> tc2 (d = sqrt(3))
            (/-1, 2/)  & ! tc1 -> ts2 (d = sqrt(3))
          ], &
          
        ! Next-Next-Nearest-neigbour (3rd) shell  (d = 2)
          [ (/ 0, 2/), & ! tc1 -> ts2 (d=2)
            (/ 2, 0/), & ! tc1 -> tc1 (d=2)
            (/-2, 0/), & ! tc1 -> tc1 (d=2)
            (/-2, 2/), & ! tc1 -> ts2 (d=2)
            (/ 1,-2/)  & ! tc1 -> s   (d=sqrt(13)/3)
          ] &
    ], shape=[3,5,2], order=[3,2,1])


end module constants
