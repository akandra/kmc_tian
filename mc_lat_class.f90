module mc_lat_class

  use constants
  use control_parameters_class

  implicit none
  !public ! for performance

  type, public :: adsorbate

    integer :: row
    integer :: col
    integer :: site
    integer :: id

  end type adsorbate


  type, public :: mc_lat

    integer :: n_rows       ! number of rows    in 2D lattice
    integer :: n_cols       ! number of columns in 2D lattice
    integer :: n_ads_sites ! number of adsorbtion site in the unit cell

    integer, dimension(:,:), allocatable  :: occupations  !  n_rows x n_cols
    integer, dimension(:,:), allocatable  :: site_type    !  n_rows x n_cols

    integer :: n_nn      ! number of nearest neighbors (6 for hex lattice)
    integer, dimension(:,:), allocatable  :: nn_list

    integer, dimension(:), allocatable :: n_ads  ! initial number of adsorbates
    type(adsorbate), dimension(:), allocatable  :: ads_list

    contains

      procedure :: print_ocs  => mc_lat_print_ocs
      procedure :: print_st   => mc_lat_print_st
      procedure :: print_ads  => mc_lat_print_ads

  end type mc_lat

  interface mc_lat

    module procedure :: mc_lat_init

  end interface


contains

  function mc_lat_init(control_pars)

    type(control_parameters), intent(in) :: control_pars
    type(mc_lat) mc_lat_init

    integer :: n_ads_total, i, j
    integer :: counter, s_counter, current_species


    mc_lat_init%n_rows = control_pars%n_rows
    mc_lat_init%n_cols = control_pars%n_cols

    ! Adsorption sites on the unit hex cell
    !  T. . . . .B1 . . . . .           1 top       (T)
    !   .  .              .  .          2 fcc       (F)
    !    .     F         .    .         3 hcp       (H)
    !     .        .   .       .        4 bridge 1  (B1)
    !      B3        B2         .       5 bridge 2  (B2)
    !       .       .     .      .      6 bridge 3  (B3)
    !        .    .          H    .
    !         . .                  .
    !          .. . . . . . . . . . .

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


    mc_lat_init%n_ads_sites = n_ads_sites

    allocate(mc_lat_init%occupations( control_pars%n_rows,&
                                      control_pars%n_cols) )
    allocate(mc_lat_init%site_type  ( control_pars%n_rows,&
                                      control_pars%n_cols) )
    allocate(mc_lat_init%n_ads(control_pars%n_species))
    ! Warning: the allocation assumes 1 ads. per unit cell
    allocate(mc_lat_init%ads_list( control_pars%n_rows*&
                                   control_pars%n_cols) )

    mc_lat_init%occupations = 0
    mc_lat_init%site_type   = 0
    mc_lat_init%n_ads       = control_pars%n_ads
    mc_lat_init%ads_list    = adsorbate(0,0,0,0)

    if (control_pars%cfg_file_name=='none') then

      n_ads_total = 0
      do i=1,size(mc_lat_init%n_ads)
        n_ads_total = n_ads_total + mc_lat_init%n_ads(i)
      end do

      counter = 0
      s_counter = 0
      current_species = 1
      loop1: do j=1,mc_lat_init%n_cols
      do i=1,mc_lat_init%n_rows
        counter   = counter + 1
        s_counter = s_counter+1
        if (s_counter>mc_lat_init%n_ads(current_species)) then
          s_counter = 1
          current_species = current_species + 1
        end if
        mc_lat_init%occupations(i,j) = counter
        ! Warning: arbitrary choice for the ads. site (top_id)!
        mc_lat_init%ads_list(counter) = adsorbate(i,j,top_id,current_species)
        if (counter == n_ads_total) exit loop1
      end do
      end do loop1

    else

!        call open_for_read(inp_unit,trim(cfg_fname))
!        read(inp_unit,*) nlat_old, nads_old
!        if (nlat_old /= nlat .OR. nads_old /= nads ) &
!            stop 'Error: inconsistent input and configuration files'
!        read(inp_unit,cfg_fmt) (occupations(i,:), i=1,nlat)
!        close(inp_unit)
!        print*
!        print*,'    Dear Sir,'
!        print*,'I would like to humbly notify you that the initial configuration is read from:'
!        print*, trim(cfg_fname)
!        print*,'Remaining your humble servant, K.M.C. Code'
!        print*

    end if

    mc_lat_init%site_type = terrace_site
    ! Define where steps and corners are
    if ( control_pars%step_period > 0) then
    do j=1,mc_lat_init%n_cols,control_pars%step_period
        mc_lat_init%site_type(:,j)   = step_site
        mc_lat_init%site_type(:,j+1) = corner_site
    end do

end if


  end function

!------------------------------------------------------------------------------
!  subroutine mc_lat_print_ocs
!  print the mc_lat occupations matrix
!
!------------------------------------------------------------------------------
  subroutine mc_lat_print_ocs (this)
    class(mc_lat), intent(in) :: this
    integer                   :: i,j

    print '(A)',' occupations:'
    do i=1,this%n_rows
      write(6,'(100i4)') (this%occupations(i,j), j=1,this%n_cols)
    end do
    print *
  end subroutine mc_lat_print_ocs

!------------------------------------------------------------------------------
!  subroutine mc_lat_print_st
!  print the mc_lat site_type matrix
!
!------------------------------------------------------------------------------
  subroutine mc_lat_print_st (this)
    class(mc_lat), intent(in) :: this
    integer                   :: i,j

    print '(A)',' site type matrix:'
    do i=1,this%n_rows
      write(6,'(100i4)') (this%site_type(i,j), j=1,this%n_cols)
    end do
    print *
  end subroutine mc_lat_print_st

!------------------------------------------------------------------------------
!  subroutine mc_lat_print_ads
!  print the mc_lat occupations matrix
!
!------------------------------------------------------------------------------
  subroutine mc_lat_print_ads (this)
    class(mc_lat), intent(in) :: this
    integer                   :: i, n_ads_total

    print '(A)', 'adsorbate list:'
    n_ads_total = 0
    do i=1,size(this%n_ads)
      n_ads_total = n_ads_total + this%n_ads(i)
    end do

    do i=1,n_ads_total
      write(6,'(5i4)') i, this%ads_list(i)
    end do

    print *

  end subroutine mc_lat_print_ads



end module mc_lat_class


