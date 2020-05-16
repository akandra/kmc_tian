module mc_lat_class

  use constants

  implicit none
  !public ! for performance


! module global variables
  real(8) :: temperature  ! temperature in K
  real(8) :: coverage     ! coverage in ML
  integer :: nsteps       ! number of Metropolis MC steps
  integer :: nsave        ! period for for conf. output

  type, public :: adsorbate

    integer :: row
    integer :: col
    integer :: site
    integer :: id

  end type adsorbate


  type, public :: mc_lat

    integer :: n_row       ! number of rows    in 2D lattice
    integer :: n_col       ! number of columns in 2D lattice
    integer :: n_ads_sites ! number of adsorbtion site in the unit cell

    integer, dimension(:,:,:), allocatable  :: occupations  !  n_row x n_col x n_ads_sites
    integer, dimension(:,:  ), allocatable  :: site_type   ! n_row x n_col

    integer :: n_nn      ! number of nearest neighbors (6 for hex lattice)
    integer, dimension(:,:), allocatable  :: nn_list

    integer :: n_ads     ! number of adsorbates
    type(adsorbate), dimension(:), allocatable  :: ads_list

    contains
      procedure :: print_ocs  => mc_lat_print_ocs
      procedure :: print_ads  => mc_lat_print_ads

  end type mc_lat

  interface mc_lat
    module procedure :: mc_lat_init
  end interface


contains

  function mc_lat_init(rows, cols, nads)

    integer, intent(in) :: rows, cols, nads
    type(mc_lat) mc_lat_init

    mc_lat_init%n_row = rows
    mc_lat_init%n_col = cols

    ! Adsorption sites on the unit hex c ell
    !  T. . . . .B1 . . . . .           1 top       (T)
    !   .  .              .  .          2 fcc       (F)
    !    .     F         .    .         3 hcp       (H)
    !     .        .   .       .        4 bridge 1  (B1)
    !      B3        B2         .       5 bridge 2  (B2)
    !       .       .     .      .      6 bridge 3  (B3)
    !        .    .          H    .
    !         . .                  .
    !          .. . . . . . . . . . .
    mc_lat_init%n_ads_sites = 6

    allocate(mc_lat_init%occupations(rows,cols,mc_lat_init%n_ads_sites))
    mc_lat_init%occupations = 0

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

    mc_lat_init%n_ads = nads

    allocate(mc_lat_init%ads_list(nads))
    mc_lat_init%ads_list = adsorbate(0,0,0,0)

  end function

!------------------------------------------------------------------------------
!  subroutine mc_lat_print_ocs
!  print the mc_lat occupations matrix
!
!------------------------------------------------------------------------------
  subroutine mc_lat_print_ocs (this)
    class(mc_lat), intent(in) :: this
    integer                   :: i,j,k

    do k=1,this%n_ads_sites
      print '(A3,A)', ads_site_names(k), ' occupations:'
      do i=1,this%n_row
        write(6,'(100i4)') (this%occupations(i,j,k), j=1,this%n_col)
      end do
    end do
    print *
  end subroutine mc_lat_print_ocs

!------------------------------------------------------------------------------
!  subroutine mc_lat_print_ads
!  print the mc_lat occupations matrix
!
!------------------------------------------------------------------------------
  subroutine mc_lat_print_ads (this)
    class(mc_lat), intent(in) :: this
    integer                   :: i

    print '(A)', 'adsorbate list:'
    write(6,'(4i4)') this%ads_list
    print *

  end subroutine mc_lat_print_ads



end module mc_lat_class


