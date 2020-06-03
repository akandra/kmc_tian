module mc_lat_class

  use constants
  use utilities
  use control_parameters_class

  implicit none

  private
  public :: mc_lat_init

  type, private :: adsorbate  ! dja changed from public 2020-06-02

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
      procedure :: hoshen_kopelman
      procedure :: cluster_size => count_cluster_sizes
      procedure, nopass, public :: mc_lat_init

  end type mc_lat

  interface mc_lat

    module procedure :: mc_lat_init

  end interface


contains

  function mc_lat_init(control_pars) result(lat)

    type(control_parameters), intent(inout) :: control_pars
    type(mc_lat)               :: lat

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

! ---------------------------------------------------------------------
! Subroutine finding clusters with Hoshen-Kopelman algorithm
! ---------------------------------------------------------------------
  subroutine hoshen_kopelman(this, species, cluster_label, largest_label)

  implicit none

  class (mc_lat), intent(in) :: this
  integer, intent(in) :: species
  integer, dimension(:,:), intent(out) :: cluster_label
  integer, intent(out) :: largest_label

  integer :: i, j, m, n, nnn2, row, col, row_new,col_new
  integer :: itemp, ip, jp, i_ads, i_ads_nn

  integer, dimension(this%n_ads(species)) :: labels
  integer, dimension(this%n_nn(1)/2) :: scanned_nn_occs, row_nn, col_nn

  nnn2 = this%n_nn(1)/2
  scanned_nn_occs = 0
  cluster_label   = 0
  do i=1,this%n_ads(species)
      labels(i) = i
  end do

  largest_label = 0
  do row=1,this%n_rows
  do col=1,this%n_cols

    i_ads = this%occupations(row,col)

    if ( i_ads > 0 .and. this%ads_list(i_ads)%id == species) then

      do m=1,nnn2

        ! Take previously scan neighbours
        call this%hop(i_ads,nnn2+m,row_nn(m),col_nn(m))
        i_ads_nn = this%occupations(row_nn(m), col_nn(m))
        if ( i_ads_nn > 0 .and. this%ads_list(i_ads_nn)%id == species ) then
            scanned_nn_occs(m) = 1
        else
            scanned_nn_occs(m) = 0
        end if

        ! Switch off periodic boundary conditions
        ! Warning: works only for hexagonal lattice,
        !          maybe, because of Corona virus
        if (row==1) scanned_nn_occs(2:3) = 0
        if (col==1) scanned_nn_occs(1) = 0
        if (col==this%n_cols) scanned_nn_occs(3) = 0

      end do

      select case (sum(scanned_nn_occs))

        case (0)
          !print*
          !print*,"Case 0"
          largest_label = largest_label + 1
          cluster_label(row,col) = largest_label

        case (1)
          !print*
          !print*,"Case 1"
          do m=1,nnn2
            if (scanned_nn_occs(m) == 1) then
              cluster_label(row,col) = &
                      lfind( cluster_label(row_nn(m),col_nn(m)), labels)
              end if
          end do

        case (2)
          !print*
          !print*,"Case 2"
          do m=1,nnn2-1
            itemp = cluster_label(row_nn(m),col_nn(m))
          do n=m+1,nnn2
            if (scanned_nn_occs(m) == 1 .and. scanned_nn_occs(n) == 1) then
              call lunion(itemp,cluster_label(row_nn(n),col_nn(n)),labels)
              cluster_label(row,col) = lfind(itemp, labels)
            end if
          end do
          end do

        case (3)
          !print*
          !print*,"Case 3"
          itemp = cluster_label(row_nn(1),col_nn(1))
          call lunion(itemp, cluster_label(row_nn(2),col_nn(2)),labels)
          call lunion(itemp, cluster_label(row_nn(3),col_nn(3)),labels)
          cluster_label(row,col) = lfind(itemp, labels)

        case default
          stop 'Hoshen-Kopelman: too many neighbors'

      end select

    end if

  end do
  end do

!  print*
!  call this%print_ocs
!  print*
!  write(*,'(8i4)') transpose(cluster_label)
!  print*
!  do i=1,largest_label
!      print('(i3,a5,i3)'), i, ' --> ' , labels(i)
!  end do

  ! Apply PBC
  do row=1,this%n_rows
    ip = cluster_label(row,1)
    if (ip > 0) then
      do m=1,nnn2
        call this%hop(this%occupations(row,1),nnn2+m,row_new,col_new)
        jp = cluster_label(row_new,col_new)
        if (jp > 0) call lunion(ip,jp,labels)
      end do
    end if
  end do

  do col=2,this%n_cols-1
    ip = cluster_label(1,col)
    if (ip > 0) then
      do m=1,nnn2
        call this%hop(this%occupations(1,col),nnn2+m,row_new,col_new)
        jp = cluster_label(row_new,col_new)
        if (jp > 0) call lunion(ip,jp,labels)
      end do
    end if
  end do

!  print*
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
!  write(*,'(8i4)') transpose(cluster_label)

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

end module mc_lat_class




