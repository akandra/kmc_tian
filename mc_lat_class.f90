module mc_lat_class

  use constants
  use utilities
  use control_parameters_class
  use energy_parameters_class

  implicit none

  private
  public :: mc_lat_init

  type, private :: adsorbate  ! dja changed from public 2020-06-02

    integer :: row  ! lattice row     occupied by an adsorbate
    integer :: col  ! lattice column  occupied by an adsorbate
    integer :: ast  ! adsorption site occupied by an adsorbate
    integer :: id   ! adsorbate species

  end type adsorbate

  type, private :: v_list

    integer, dimension(:),allocatable :: list

  end type


  type, public :: mc_lat

    integer :: n_rows       ! number of rows    in 2D lattice
    integer :: n_cols       ! number of columns in 2D lattice
    integer :: n_max_ads_sites ! number of adsorbtion sites in the unit cell

    integer, dimension(:,:), allocatable  :: occupations  !  n_rows x n_cols
    ! array of lattice site types
    integer, dimension(:,:), allocatable  :: lst          !  n_rows x n_cols
    ! lattice vectors
    real(dp), dimension(2) :: lat_vec_1, lat_vec_2
    ! shellwise number of neighbors
    integer, dimension(n_shells) :: n_nn
    integer, dimension(:,:,:), allocatable  :: shell_list
    ! List of additional nn directions to scan after reaction
    integer, dimension(:,:), allocatable :: nn_new
    ! number of adsorbates per species
    integer, dimension(:), allocatable :: n_ads
    ! initial and current adsorbates lists
    type(adsorbate), dimension(:), allocatable  :: ads_list_ini, ads_list
    ! Available adsorbtion sites list (n_species x n_max_lat_site_types x (n_avail_ads_sites))
    type(v_list), dimension(:,:), allocatable :: avail_ads_sites

    contains

      procedure :: print_ocs  => mc_lat_print_ocs
      procedure :: print_st   => mc_lat_print_st
      procedure :: print_ads  => mc_lat_print_ads
      procedure :: hop        => mc_lat_hop_with_pbc
      procedure :: neighbor   => mc_lat_neighbor_with_pbc
      procedure :: distance   => mc_lat_distance_with_pbc
      procedure :: n_ads_tot  => mc_lat_n_ads_total
      procedure :: hoshen_kopelman
      procedure :: cluster_size => count_cluster_sizes
      procedure :: rdf => radial_distribution_function
      procedure :: restore    => mc_lat_restore
      procedure, nopass, public :: mc_lat_init

  end type mc_lat


contains

  function mc_lat_init(control_pars, energy_pars) result(lat)

    type(control_parameters), intent(inout) :: control_pars
    type(energy_parameters),     intent(in) :: energy_pars
    type(mc_lat)                            :: lat

    integer :: i, j, k, m, ios
    integer :: counter, s_counter, current_species
    integer :: n_rows_in, n_cols_in, step_period_in
    character(len=max_string_length) :: line
    character(len=max_string_length) :: ads_names_in(100)
    integer :: nw, n_species_in, n_site_types, ads_site
    integer, allocatable :: n_ads_in(:)

    lat%n_rows = control_pars%n_rows
    lat%n_cols = control_pars%n_cols

    ! Adsorption sites on the unit hex cell
    !  T---------B1---------.           1 top       (T)
    !   \  .              .  .          2 fcc       (F)
    !    \     F         .    .         3 hcp       (H)
    !     \        .   .       .        4 bridge 1  (B1)
    !      B3        B2         .       5 bridge 2  (B2)
    !       \       .     .      .      6 bridge 3  (B3)
    !        \    .          H    .
    !         \ .                  .
    !          \. . . . . . . . . . .

    ! lattice vectors for hex lattice
    lat%lat_vec_1 = [ 1.0_dp, 0.0_dp]
    lat%lat_vec_2 = [ 0.5_dp,-cos(pi/3)]

    ! Shellwise number of neighbors for hex lattice
    lat%n_nn = [6,6,6]

    ! NN list for the hexagonal structure
    !  11    12    13***  14**  15***
    !
    !     21    22**   23*   24*   25**
    !
    !        31*** 32*   33    34*   35***
    !
    !           41**  42*   43*   44**  45
    !
    !              51*** 52**  53*** 54    55

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

    ! Check the distances to the neighbours
!    print*, sqrt( &
!                 ( lat%shell_list(2,:,1)*cos(pi/3.0_dp) &
!                  +lat%shell_list(2,:,2)               )**2    &
!                +                                                    &
!                 ( lat%shell_list(2,:,1)*sin(pi/3.0_dp))**2    &
!        )
!    stop

    ! List of additional nn directions for hex lattice to scan after a reaction
    allocate(lat%nn_new( lat%n_nn(1), lat%n_nn(1)/2) )
    do m=1,lat%n_nn(1)
    do i=1,lat%n_nn(1)/2
      lat%nn_new(m,i) = modulo( m + i - lat%n_nn(1)/2, lat%n_nn(1) ) + 1
    end do
    end do

    lat%n_max_ads_sites = n_max_ads_sites

    ! Allocate arrays
    allocate(lat%occupations(control_pars%n_rows,&
                             control_pars%n_cols))
    allocate(lat%lst  (control_pars%n_rows,&
                             control_pars%n_cols))
    allocate(lat%n_ads(control_pars%n_species))
    ! Warning: the allocation assumes 1 ads. per unit cell
    allocate(lat%ads_list(    control_pars%n_rows*control_pars%n_cols))
    allocate(lat%ads_list_ini(control_pars%n_rows*control_pars%n_cols))
    ! Initialize arrays
    lat%occupations  = 0
    lat%n_ads        = control_pars%n_ads
    lat%ads_list     = adsorbate(0,0,0,0)
    lat%ads_list_ini = adsorbate(0,0,0,0)

    lat%lst = terrace_site
    ! Define where steps and corners are
    if ( control_pars%step_period > 0) then
      do j=1,lat%n_cols,control_pars%step_period
          lat%lst(:,j)   = step_site
          lat%lst(:,j+1) = corner_site
      end do
      n_site_types = n_max_lat_site_types
    else
      n_site_types = 1
    end if

    allocate(lat%avail_ads_sites(control_pars%n_species,n_site_types))

    do i=1,control_pars%n_species
    do j=1,n_site_types
      counter = 0
      do k=1,n_max_ads_sites
        if (energy_pars%ads_energy(i,j,k) < energy_pars%undefined_energy) then
          counter = counter + 1
        end if
      end do
      allocate(lat%avail_ads_sites(i,j)%list(counter))
      counter = 0
      do k=1,n_max_ads_sites
        if (energy_pars%ads_energy(i,j,k) < energy_pars%undefined_energy) then
          counter = counter + 1
          lat%avail_ads_sites(i,j)%list(counter) = k
        end if
      end do
    end do
    end do

!!----------------------------------------------------------------------------------
!! debug printout for construction of avail_ads_sites
!print *, 'debug printout for construction of avail_ads_sites mc_lat_class line 191'
!do i=1,control_pars%n_species
!do j=1,n_site_types
!  write(*, '(A, A5, A, A10)',advance='no')' species = ',control_pars%ads_names(i),'site type =',lat_site_names(j)
!  print '(A,A5,1x,A,1x,A,1x,A,1x,A,1x)','  list = ',(ads_site_names(lat%avail_ads_sites(i,j)%list(k)), &
!                    k=1,size(lat%avail_ads_sites(i,j)%list))
!  end do
!end do
!----------------------------------------------------------------------------------

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
        ! Warning: arbitrary choice for the ads. site (the 1st available)!
        ads_site = lat%avail_ads_sites(current_species,lat%lst(i,j))%list(1)
        lat%ads_list(counter) = adsorbate(i,j,ads_site,current_species)
        if (counter == lat%n_ads_tot()) exit loop1
      end do
      end do loop1

    else
      ! Read info from the configuration file
      call open_for_read(inp_unit,trim(control_pars%cfg_file_name))
      read(inp_unit,'(A)') line
      read(inp_unit,*) n_rows_in, n_cols_in, step_period_in
      ! Check consistency with control parameters
      if (n_rows_in /= control_pars%n_rows .OR.&
          n_cols_in /= control_pars%n_cols .OR.&
          step_period_in /= control_pars%step_period) then
        print*, 'Error! '
        print*, ' inconsistency in input and configuration files:'
        print*, ' lattice definition'
        stop
      end if

      read(inp_unit,'(A)') line
      call split_string (line, ads_names_in, n_species_in)
      allocate(n_ads_in(n_species_in))

      ios = 0
      do while (ios == 0)
        read(inp_unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        read(inp_unit,*) n_ads_in
        ! Check consistency with control parameters
        if (n_species_in /= control_pars%n_species &
        !.OR. any(n_ads_in /= control_pars%n_ads) --- Get rid of that line when our wisdom get improved
          ) then
          print*, 'Error! '
          print*, ' inconsistency in input and configuration files:'
          print*, ' adsorbates definition'
          stop
        end if
        ! Read the configuration
        do i=1,lat%n_ads_tot()
          read(inp_unit, *, iostat=ios) j,lat%ads_list(i)
          if (ios /= 0) then
            stop "Error in control file: premature end of the final configuration"
          end if
        end do
      end do ! ios

      close(inp_unit)

      ! Fill the lattice with adsorbates
      do i=1,lat%n_ads_tot()
        lat%occupations(lat%ads_list(i)%row, lat%ads_list(i)%col) = i
      end do

    end if

    ! Save initial adsorbate list
    lat%ads_list_ini = lat%ads_list

  end function mc_lat_init

! ---------------------------------------------------------------------
! Subroutine restoring mc_lat structure from initial adsorbate list
! ---------------------------------------------------------------------
  subroutine mc_lat_restore(this)

    class(mc_lat), intent(inout) :: this

    integer, dimension(size(this%n_ads)) :: n_ads_counter
    integer :: i, ads_id

    ! Restore the ads list
    this%ads_list = this%ads_list_ini
    ! Restore species-specific numbers of adsorbates
    n_ads_counter = 0
    do i=1,size(this%ads_list)
      ads_id = this%ads_list(i)%id
      if (ads_id>0) n_ads_counter(ads_id) = n_ads_counter(ads_id) + 1
    end do
    this%n_ads = n_ads_counter
    ! Clean up and restore occupations
    this%occupations = 0
    do i=1,this%n_ads_tot()
      this%occupations(this%ads_list(i)%row, this%ads_list(i)%col) = i
    end do

  end subroutine

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
      write(6,'(100i4)') (this%lst(i,j), j=1,this%n_cols)
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
  subroutine mc_lat_hop_with_pbc(this,i,ihop, row, col, ads_site)

    class(mc_lat), intent(in) :: this
    integer, intent(in)  :: i, ihop
    integer, intent(out) :: row, col, ads_site
    integer :: list_size, id, site_type


    row = modulo(this%ads_list(i)%row &
               + this%shell_list(1,ihop,1) - 1, this%n_rows) + 1
    col = modulo(this%ads_list(i)%col &
               + this%shell_list(1,ihop,2) - 1, this%n_cols) + 1

    id = this%ads_list(i)%id
    site_type = this%lst(row,col)
    list_size = size( this%avail_ads_sites(id,site_type)%list )
    if ( list_size > 1 ) then
      ads_site = this%avail_ads_sites(id,site_type)%list(irand(list_size))
    else
      ads_site = this%avail_ads_sites(id,site_type)%list(1)
    end if

  end subroutine

! ---------------------------------------------------------------------
! Subroutine finding  position of neighbor accounting for PBCs
! ---------------------------------------------------------------------
  subroutine mc_lat_neighbor_with_pbc(this,i,ihop, row, col, shell)

    class(mc_lat), intent(in) :: this
    integer, intent(in)  :: i, ihop
    integer, intent(out) :: row, col
    integer, intent(in), optional :: shell

    integer :: ishell

    ishell = 1
    if (present(shell)) ishell = shell
    row = modulo(this%ads_list(i)%row &
               + this%shell_list(ishell,ihop,1) - 1, this%n_rows) + 1
    col = modulo(this%ads_list(i)%col &
               + this%shell_list(ishell,ihop,2) - 1, this%n_cols) + 1

  end subroutine

! ---------------------------------------------------------------------
! Subroutine applying PBCs to calculate distance
! ---------------------------------------------------------------------
  subroutine mc_lat_distance_with_pbc(this, ads1, ads2, r)

    class(mc_lat), intent(in) :: this
    integer, intent(in)   :: ads1, ads2
    real(dp), intent(out) :: r

    integer :: d_row, d_col

    d_row = this%ads_list(ads1)%row - this%ads_list(ads2)%row
    d_col = this%ads_list(ads1)%col - this%ads_list(ads2)%col

    itemp = (2*d_row)/this%n_rows
    d_row = d_row - itemp*this%n_rows/2
    itemp = (2*d_col)/this%n_cols
    d_col = d_col - itemp*this%n_cols/2

    r = sqrt( (this%lat_vec_1(1)*d_col + this%lat_vec_2(1)*d_col)**2 &
             +(this%lat_vec_1(2)*d_row + this%lat_vec_2(2)*d_row)**2  )

  end subroutine

! ---------------------------------------------------------------------
! Function returning the total number of adsorbates
! ---------------------------------------------------------------------
  integer function mc_lat_n_ads_total(this) result(n_ads_total)

    class(mc_lat), intent(in) :: this

    integer :: i

      n_ads_total = 0
      do i=1,size(this%n_ads)
        n_ads_total = n_ads_total + this%n_ads(i)
      end do

  end function

! ---------------------------------------------------------------------
! Subroutine finding clusters with Hoshen-Kopelman algorithm
! ---------------------------------------------------------------------
  subroutine hoshen_kopelman(this, species, cluster_label, largest_label)

  implicit none

  class (mc_lat), intent(in) :: this
  integer, intent(in) :: species
  integer, dimension(:,:), intent(out) :: cluster_label
  integer, intent(out) :: largest_label

  type nn_position
    integer :: shell,m
  end type nn_position

  integer :: i, j, m, n, row, col, row_new,col_new, shell
  integer :: itemp, ip, jp, i_ads, i_ads_nn, n_neighbors
  integer :: ads_site

  integer, dimension(n_shells) :: nnn2
  integer, dimension(this%n_ads(species)) :: labels
  integer, dimension(n_shells, maxval(this%n_nn)/2) :: scanned_nn_occs, row_nn, col_nn
  type(nn_position), dimension(sum(this%n_nn)/2) :: scanned_nn_list


  nnn2 = this%n_nn/2
  scanned_nn_occs = 0
  cluster_label   = 0
  do i=1,this%n_ads(species)
      labels(i) = i
  end do

  largest_label = 0
  do row=1,this%n_rows
  do col=1,this%n_cols

    i_ads = this%occupations(row,col)

    if ( i_ads > 0 ) then
      if ( this%ads_list(i_ads)%id == species ) then

        do shell=1,n_shells
        do m=1,nnn2(shell)

          ! Take previously scanned neighbours
          call this%neighbor(i_ads,nnn2(shell)+m,row_nn(shell,m),col_nn(shell,m), shell)
          i_ads_nn = this%occupations(row_nn(shell,m), col_nn(shell,m))

          scanned_nn_occs(shell,m) = 0
          if ( i_ads_nn > 0 ) then
            if ( this%ads_list(i_ads_nn)%id == species ) scanned_nn_occs(shell,m) = 1
          end if

        end do
        end do

        ! Switch off periodic boundary conditions
        ! Warning: works only for hexagonal lattice, because of Corona virus

        ! 1st shell
        if (row==1) scanned_nn_occs(1,2:3) = 0
        if (col==1) scanned_nn_occs(1,1) = 0
        if (col==this%n_cols) scanned_nn_occs(1,3) = 0

        ! 2nd shell
        if (row==1) scanned_nn_occs(2,:) = 0
        if (row==2) scanned_nn_occs(2,2) = 0
        if (col==1) scanned_nn_occs(2,1) = 0
        if (col==this%n_cols) scanned_nn_occs(2,2:3) = 0
        if (col==this%n_cols-1) scanned_nn_occs(2,3) = 0

        ! 3rd shell
        if (row==1 .or. row==2) scanned_nn_occs(3,2:3) = 0
        if (col==1 .or. col==2) scanned_nn_occs(3,1) = 0
        if (col==this%n_cols .or. col==this%n_cols-1) scanned_nn_occs(3,3) = 0

        ! create list of neighbors
        n_neighbors = 0
        scanned_nn_list = nn_position(0,0)
        do shell=1,n_shells
        do m=1,nnn2(shell)
          if ( scanned_nn_occs(shell,m) == 1 ) then
            n_neighbors = n_neighbors + 1
            scanned_nn_list(n_neighbors)%shell = shell
            scanned_nn_list(n_neighbors)%m     = m
          end if
        end do
        end do

       select case (n_neighbors)

          case (0)
!            print*
!            print*,"Case 0"
            largest_label = largest_label + 1
            cluster_label(row,col) = largest_label

          case (1)
!            print*
!            print*,"Case 1"
            shell = scanned_nn_list(1)%shell
            m     = scanned_nn_list(1)%m
            cluster_label(row,col) = lfind( cluster_label(row_nn(shell,m),col_nn(shell,m)), labels)

          case default
            shell = scanned_nn_list(1)%shell
            m     = scanned_nn_list(1)%m
            itemp = cluster_label(row_nn(shell,m),col_nn(shell,m))
            do i=2,n_neighbors
              shell = scanned_nn_list(i)%shell
              n     = scanned_nn_list(i)%m
              call lunion(itemp,cluster_label(row_nn(shell,n),col_nn(shell,n)),labels)
            end do
            cluster_label(row,col) = lfind(itemp, labels)

        end select

      end if
    end if

  end do
  end do

!  print*
!  call this%print_ocs
!  print*
!  print*,'cluster labels before applying PBC:'
!  write(*,'(20i4)') transpose(cluster_label)
!  print*,'cluster labels list:'
!  do i=1,largest_label
!      print('(i3,a5,i3)'), i, ' --> ' , labels(i)
!  end do

  ! Apply PBC
  do row=1,this%n_rows
  do col=1,2
    ip = cluster_label(row,col)
    if (ip > 0) then
      do shell=1,n_shells
      do m=1,this%n_nn(shell)
        call this%neighbor(this%occupations(row,col),m,row_new,col_new,shell)
        jp = cluster_label(row_new,col_new)
        if (jp > 0) call lunion(ip,jp,labels)
      end do
      end do
    end if
  end do
  end do

  do col=3,this%n_cols-2
  do row=1,2
    ip = cluster_label(row,col)
    if (ip > 0) then
      do shell=1, n_shells
      do m=1,this%n_nn(shell)
        call this%neighbor(this%occupations(row,col),m,row_new,col_new,shell)
        jp = cluster_label(row_new,col_new)
        if (jp > 0) call lunion(ip,jp,labels)
      end do
      end do
    end if
  end do
  end do

!  print*
!  print*,'cluster labels list after PBC:'
!  do i=1,largest_label
!      print('(i3,a5,i3)'), i, ' --> ' , labels(i)
!  end do

  ! Going down to roots
  do row=1,this%n_rows
  do col=1,this%n_cols
    if (cluster_label(row,col) > 0) &
       cluster_label(row,col) = lfind(cluster_label(row,col), labels)
  end do
  end do

!  print*
!  print*,'cluster labels after PBC and rooting:'
!  write(*,'(20i4)') transpose(cluster_label)
!  pause

  contains

    integer function lfind(x, labels)

      integer, intent(inout) :: x
      integer, dimension(:), intent(inout) :: labels

      integer :: y, z

      y = x
        do while (labels(y) /= y)
            y = labels(y)
        end do

    !    do while (labels(x) /= x)
    !        z = labels(x)
    !        labels(x) = y
    !        x = z
    !    end do
    !
        lfind = y

    end function

    subroutine lunion(x, y, labels)

      integer, intent(inout) :: x, y
      integer, dimension(:), intent(inout) :: labels

        if (x > y) then
            labels( lfind(x, labels) ) = lfind(y, labels)
        else if ( x < y) then
            labels( lfind(y, labels) ) = lfind(x, labels)
        endif

    end subroutine

  end subroutine

! ---------------------------------------------------------------------
! Subroutine counting cluster sizes
! ---------------------------------------------------------------------
  subroutine count_cluster_sizes(this, species, cluster_label, cluster_sizes)

  implicit none

  class (mc_lat), intent(in) :: this
  integer, intent(in) :: species
  integer, dimension(:,:), intent(in) :: cluster_label
  integer, dimension(:),  intent(out) :: cluster_sizes
  integer :: i, i_label

    cluster_sizes   = 0
    do i=1, this%n_ads_tot()
      if (this%ads_list(i)%id == species) then
        i_label = cluster_label(this%ads_list(i)%row,this%ads_list(i)%col)
        cluster_sizes(i_label) = cluster_sizes(i_label) + 1
      end if
    end do

  end subroutine

! ---------------------------------------------------------------------
! Subroutine calculating radial distribution function
! ---------------------------------------------------------------------
  subroutine radial_distribution_function(this, rdf)

  implicit none

  class (mc_lat), intent(in) :: this
  integer, dimension(:),  intent(out) :: rdf

  integer :: ads1, row1, col1
  integer :: ads2, row2, col2
  integer :: n_ads_total

  n_ads_total = this%n_ads_tot()
  do ads1 = 1, n_ads_total-1
  do ads2 = ads1+1, n_ads_total

    row1 = this%ads_list(ads1)%row
    col1 = this%ads_list(ads1)%col

    row2 = this%ads_list(ads2)%row
    col2 = this%ads_list(ads2)%col

    d_row = row2 - row1
    d_col = col2 - col1

  end do
  end do

  end subroutine

end module mc_lat_class




