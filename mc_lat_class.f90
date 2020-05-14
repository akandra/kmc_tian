module mc_lat_class

  use constants

  implicit none
  !public ! for performance


! module global variables
  real(8) :: temperature  ! temperature in K
  real(8) :: coverage     ! coverage in ML
  real(8) :: eps          ! O-O interaction energy in eV
  integer :: nsteps       ! number of Metropolis MC steps
  integer :: nsave        ! period for for conf. output

  type, public :: adsorbate

    integer :: i
    integer :: j
    integer :: site

  end type adsorbate

  type, public :: mc_lat
    integer       :: nlat         ! size of 2D lattice (nlat x nlat)
    integer       :: nads         ! number of adsorbates
    integer       :: nnn          ! number of nearest neighbors (6 for hex lattice)


    integer, dimension(:,:), allocatable  :: occupations  ! (nlat x nlat) x nads_sites
    integer, dimension(:,:  ), allocatable  :: site_type    ! (nlat x nlat)

    integer, dimension(:,:), allocatable  :: clusters           ! (nlat x nlat)

    integer, dimension(:,:), allocatable  :: ads_list, nn_list

    contains
      procedure :: print_ocs  => mc_lat_print_ocs
      procedure :: print_clus => mc_lat_print_clus

  end type mc_lat

  interface mc_lat
    module procedure :: mc_lat_init
  end interface


contains

  function mc_lat_init(n)

    integer, intent(in) :: n
    type(mc_lat) mc_lat_init

    mc_lat_init%nlat = n
    allocate(mc_lat_init%occupations(n,n))
    allocate(mc_lat_init%clusters(n,n))
    mc_lat_init%occupations = 0
    mc_lat_init%clusters    = 0
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
    do i=1,this%nlat
      write(6,'(100i4)') (this%occupations(i,j), j=1,this%nlat)
    end do
    print *
  end subroutine mc_lat_print_ocs


!------------------------------------------------------------------------------
!  subroutine mc_lat_print_clus
!    print the clusters matrix
!------------------------------------------------------------------------------

  subroutine mc_lat_print_clus (this)
    class(mc_lat), intent(in) :: this
    integer                   :: i,j

    print '(/A)','clusters'
    do i=1,this%nlat
      write(6,'(100i4)') (this%clusters(i,j), j=1,this%nlat)
    end do
    print *
  end subroutine mc_lat_print_clus
!------------------------------------------------------------------------------



end module mc_lat_class


