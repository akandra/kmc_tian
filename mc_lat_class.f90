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
    ! shellwise number of neighbors
    integer, dimension(n_shells) :: n_nn
    integer, dimension(:,:,:), allocatable  :: shell_list

    integer, dimension(:), allocatable :: n_ads  ! initial number of adsorbates
    type(adsorbate), dimension(:), allocatable  :: ads_list

    contains

      procedure :: print_ocs  => mc_lat_print_ocs
      procedure :: print_st   => mc_lat_print_st
      procedure :: print_ads  => mc_lat_print_ads
      procedure :: hop        => mc_lat_hop_with_pbc
      procedure :: n_ads_tot  => mc_lat_n_ads_total

  end type mc_lat

  interface mc_lat

    module procedure :: mc_lat_init

  end interface


contains

  function mc_lat_init(control_pars) result(lat)

    type(control_parameters), intent(inout) :: control_pars
    type(mc_lat) lat

    integer :: i, j
    integer :: counter, s_counter, current_species
    integer :: n_rows_in, n_cols_in, step_period_in
    character(len=max_string_length) :: line
    character(len=max_string_length) :: ads_names_in(100)
    integer :: nw, n_species_in
    integer, allocatable :: n_ads_in(:)

    lat%n_rows = control_pars%n_rows
    lat%n_cols = control_pars%n_cols

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



    ! Shellwise number of neighbors for hex lattice
    lat%n_nn = [6,6,6]

    ! NN list for the hexagonal structure
    !  11    12*   13*   14
    !
    !     21*   22*   23*   24
    !
    !        31*   32*   33    34
    !
    !           41    42    43    44

    allocate(lat%shell_list(n_shells,maxval(lat%n_nn),2))
    ! Nearest-neigbour (1st) shell (d = 1))
    lat%shell_list(1,1,:) = (/ 0, 1/)
    lat%shell_list(1,2,:) = (/ 1, 0/)
    lat%shell_list(1,3,:) = (/ 1,-1/)
    lat%shell_list(1,4,:) = (/ 0,-1/)
    lat%shell_list(1,5,:) = (/-1, 0/)
    lat%shell_list(1,6,:) = (/-1, 1/)
    ! Next-Nearest-neigbour (2nd) shell  (d = sqrt(3))
    lat%shell_list(2,1,:) = (/ 1, 1/)
    lat%shell_list(2,2,:) = (/ 2,-1/)
    lat%shell_list(2,3,:) = (/ 1,-2/)
    lat%shell_list(2,4,:) = (/-1,-1/)
    lat%shell_list(2,5,:) = (/-2, 1/)
    lat%shell_list(2,6,:) = (/-1, 2/)
    ! Next-Next-Nearest-neigbour (3rd) shell  (d = 2)
    lat%shell_list(3,1,:) = (/ 0, 2/)
    lat%shell_list(3,2,:) = (/ 2, 0/)
    lat%shell_list(3,3,:) = (/ 2,-2/)
    lat%shell_list(3,4,:) = (/ 0,-2/)
    lat%shell_list(3,5,:) = (/-2, 0/)
    lat%shell_list(3,6,:) = (/-2, 2/)

    ! Check the distances to the neibours
!    print*, sqrt( &
!                 ( lat%shell_list(2,:,1)*cos(pi/3.0_dp) &
!                  +lat%shell_list(2,:,2)               )**2    &
!                +                                                    &
!                 ( lat%shell_list(2,:,1)*sin(pi/3.0_dp))**2    &
!        )
!    stop

    lat%n_ads_sites = n_ads_sites

    ! Allocate arrays
    allocate(lat%occupations(control_pars%n_rows,&
                             control_pars%n_cols))
    allocate(lat%site_type  (control_pars%n_rows,&
                             control_pars%n_cols))
    allocate(lat%n_ads(control_pars%n_species))
    ! Warning: the allocation assumes 1 ads. per unit cell
    allocate(lat%ads_list(control_pars%n_rows*control_pars%n_cols))
    ! Initialize arrays
    lat%occupations = 0
    lat%site_type   = 0
    lat%n_ads       = control_pars%n_ads
    lat%ads_list    = adsorbate(0,0,0,0)

    if (control_pars%cfg_file_name=='none') then
      ! Populate the lattice
      counter = 0
      s_counter = 0
      current_species = 1
      loop1: do j=1,lat%n_cols
      do i=1,lat%n_rows
        counter   = counter + 1
        s_counter = s_counter+1
        if (s_counter>lat%n_ads(current_species)) then
          s_counter = 1
          current_species = current_species + 1
        end if
        lat%occupations(i,j) = counter
        ! Warning: arbitrary choice for the ads. site (top_id)!
        lat%ads_list(counter) = adsorbate(i,j,top_id,current_species)
        if (counter == lat%n_ads_tot()) exit loop1
      end do
      end do loop1

    else
      ! Read info from the configuration file
      call open_for_read(inp_unit,trim(control_pars%cfg_file_name))
      read(inp_unit,'(A)') line
      read(inp_unit,*) n_rows_in, n_cols_in, step_period_in
      read(inp_unit,'(A)') line
      call split_string (line, ads_names_in, n_species_in)
      allocate(n_ads_in(n_species_in))
      read(inp_unit,*) n_ads_in

      ! Check consistency with control parameters
      if (n_rows_in /= control_pars%n_rows .OR.&
          n_cols_in /= control_pars%n_cols .OR.&
          step_period_in /= control_pars%step_period) then
        print*, 'Error! '
        print*, ' inconsistency in input and configuration files:'
        print*, ' lattice definition'
        stop
      end if
      if (n_species_in /= control_pars%n_species .OR.&
          any(n_ads_in /= control_pars%n_ads)) then
        print*, 'Error! '
        print*, ' inconsistency in input and configuration files:'
        print*, ' adsorbates definition'
        stop
      end if

      ! Read the configuration
      read(inp_unit,'(A)') line
      do i=1,lat%n_ads_tot()
        read(inp_unit,*) j,lat%ads_list(i)
        lat%occupations(lat%ads_list(i)%row,&
                        lat%ads_list(i)%col) = i
      end do

      close(inp_unit)

    end if

    lat%site_type = terrace_site
    ! Define where steps and corners are
    if ( control_pars%step_period > 0) then
      do j=1,lat%n_cols,control_pars%step_period
          lat%site_type(:,j)   = step_site
          lat%site_type(:,j+1) = corner_site
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
  subroutine mc_lat_print_ads (this, out_unit)
    class(mc_lat), intent(in) :: this
    integer, optional         :: out_unit

    integer                   :: i, n_ads_total
    integer                   :: my_unit

    if (present(out_unit)) then
      my_unit = out_unit
    else
      my_unit = output_unit
      print '(A)', 'adsorbate list:'
    end if

    n_ads_total = 0
    do i=1,size(this%n_ads)
      n_ads_total = n_ads_total + this%n_ads(i)
    end do

    do i=1,n_ads_total
      write(my_unit,'(5i10)') i, this%ads_list(i)
    end do

  end subroutine mc_lat_print_ads

!------------------------------------------------------------------------------
!  subroutine pbc
!  applies pbc
!
!------------------------------------------------------------------------------
  subroutine mc_lat_hop_with_pbc(this,i,ihop, row, col)

    class(mc_lat), intent(in) :: this
    integer, intent(in)  :: i, ihop
    integer, intent(out) :: row, col

    row = modulo(this%ads_list(i)%row &
               + this%shell_list(1,ihop,1) - 1, this%n_rows) + 1
    col = modulo(this%ads_list(i)%col &
               + this%shell_list(1,ihop,2) - 1, this%n_cols) + 1

  end subroutine

  integer function mc_lat_n_ads_total(this) result(n_ads_total)

    class(mc_lat), intent(in) :: this

    integer :: i

      n_ads_total = 0
      do i=1,size(this%n_ads)
        n_ads_total = n_ads_total + this%n_ads(i)
      end do

  end function


end module mc_lat_class


