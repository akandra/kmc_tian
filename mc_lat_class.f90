module mc_lat_class

  use constants

  implicit none
  !public ! for performance


! module global variables
  real(8) :: temperature  ! temperature in K
  real(8) :: coverage     ! coverage in ML
  integer :: nsteps       ! number of Metropolis MC steps
  integer :: nsave        ! period for for conf. output

!  type, public :: adsorbate
!
!    integer :: i
!    integer :: j
!    integer :: site
!    character(len=4) :: name
!
!  end type adsorbate


  type, public :: mc_lat

    integer :: n_row       ! number of rows    in 2D lattice
    integer :: n_col       ! number of columns in 2D lattice
    integer :: n_ads_sites ! number of columns in 2D lattice

    integer, dimension(:,:,:), allocatable  :: occupations  !  n_row x n_col x n_ads_sites
    integer, dimension(:,:  ), allocatable  :: site_type   ! n_row x n_col

    integer :: n_nn      ! number of nearest neighbors (6 for hex lattice)
    integer, dimension(:,:), allocatable  :: nn_list

    integer :: n_ads     ! number of adsorbates



    integer, dimension(:,:), allocatable  :: ads_list

    contains
      procedure :: print_ocs  => mc_lat_print_ocs

  end type mc_lat

  interface mc_lat
    module procedure :: mc_lat_init
  end interface


contains

  function mc_lat_init(rows, cols)

    integer, intent(in) :: rows, cols
    type(mc_lat) mc_lat_init

    mc_lat_init%n_row = rows
    mc_lat_init%n_col = cols

    ! Adsorption sites on the unit hex cell
    !  T.........B1..........           1 top       (T)
    !   .  .              .  .          2 fcc       (F)
    !    .     F         .    .         3 hcp       (H)
    !     .        .   .       .        4 bridge 1  (B1)
    !      B3        B2         .       5 bridge 2  (B2)
    !       .       .     .      .      6 bridge 3  (B3)
    !        .    .          H    .
    !         . .                  .
    !          ......................
    mc_lat_init%n_ads_sites = 6

    allocate(mc_lat_init%occupations(rows,cols))
    mc_lat_init%occupations(i,j) = 0

    mc_lat_init%n_ads = 0

    mc_lat_init%n_nn = 6

    ! NN list for the hexagonal structure
    !  11    12*   13*   14
    !
    !     21*   22*   23*   24
    !
    !        31*   32*   33    34
    !
    !           41    42    43    44

    allocate(mc_lat_init%nn_list(mc_lat_init%n_nn,2))
    mc_lat_init%nn_list(1,:) = (/ 0, 1/)
    mc_lat_init%nn_list(2,:) = (/ 1, 0/)
    mc_lat_init%nn_list(3,:) = (/ 1,-1/)
    mc_lat_init%nn_list(4,:) = (/ 0,-1/)
    mc_lat_init%nn_list(5,:) = (/-1, 0/)
    mc_lat_init%nn_list(6,:) = (/-1, 1/)


  end function

!------------------------------------------------------------------------------
!  subroutine mc_lat_print_ocs
!  print the mc_lat occupations matrix
!
!------------------------------------------------------------------------------
  subroutine mc_lat_print_ocs (this)
    class(mc_lat), intent(in) :: this
    integer                   :: i,j


    print '(/A)','occupations'
    do i=1,this%n_row
      write(6,'(100i4)') (this%occupations(i,j), j=1,this%n_col)
    end do
    print *
  end subroutine mc_lat_print_ocs




end module mc_lat_class


