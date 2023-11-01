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
    integer, dimension(n_max_lat_site_types,n_shells) :: n_nn
    integer, dimension(:,:,:,:), allocatable  :: shell_list
    ! List of additional nn directions to scan after reaction
    integer, dimension(:,:), allocatable :: nn_new
    ! number of adsorbates per species
    integer, dimension(:), allocatable :: n_ads
    ! initial and current adsorbates lists
    type(adsorbate), dimension(:), allocatable  :: ads_list
    ! Available adsorbtion sites list (n_species x n_max_lat_site_types x (n_avail_ads_sites))
    type(v_list), dimension(:,:), allocatable :: avail_ads_sites
    integer :: n_ini_confs  ! Number of available confs in the start-conf file

    contains

      procedure :: print_ocs  => mc_lat_print_ocs
      procedure :: print_st   => mc_lat_print_st
      procedure :: print_ads  => mc_lat_print_ads
      procedure :: hop        => mc_lat_hop_with_pbc
      procedure :: neighbor   => mc_lat_neighbor_with_pbc
      procedure :: distance   => mc_lat_distance_with_pbc
      procedure :: n_ads_tot  => mc_lat_n_ads_total
!      procedure :: hoshen_kopelman
      procedure :: cluster_size => count_cluster_sizes
      procedure :: rdf_hist
      procedure :: conf_init  => mc_lat_conf_init
      procedure, nopass, public :: mc_lat_init

  end type mc_lat


contains

  function mc_lat_init(c_pars, e_pars) result(lat)

    type(control_parameters), intent(inout) :: c_pars
    type(energy_parameters),     intent(in) :: e_pars
    type(mc_lat)                            :: lat

    integer :: i, j, k, m
    integer :: n_site_types, counter

    lat%n_rows = c_pars%n_rows
    lat%n_cols = c_pars%n_cols

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
    lat%lat_vec_1 = [ cos(0.0_dp),    -sin(0.0_dp)   ]
    lat%lat_vec_2 = [ cos(pi/3.0_dp), -sin(pi/3.0_dp)]

    ! Shellwise number of neighbors for hex lattice
    lat%n_nn(  terrace_site,:) = [6,6,6]
    lat%n_nn(     step_site,:) = [5,5,5]
    lat%n_nn(   corner_site,:) = [4,4,6]
    lat%n_nn(terrace_s_site,:) = [6,7,4]
    lat%n_nn(terrace_c_site,:) = [5,5,5]

    allocate(lat%shell_list(n_max_lat_site_types, n_shells, maxval(lat%n_nn), 2))

    ! terrace

    ! NN list for the hexagonal 111-structure (row,col)
    !  11    12    13***  14**  15***
    !
    !     21    22**   23*   24*   25**
    !
    !        31*** 32*   33    34*   35***
    !
    !           41**  42*   43*   44**  45
    !
    !              51*** 52**  53*** 54    55

    ! Nearest-neigbour (1st) shell (d = 1))
    lat%shell_list(terrace_site,1,1,:) = (/ 0, 1/)
    lat%shell_list(terrace_site,1,2,:) = (/ 1, 0/)
    lat%shell_list(terrace_site,1,3,:) = (/ 1,-1/)
    lat%shell_list(terrace_site,1,4,:) = (/ 0,-1/)
    lat%shell_list(terrace_site,1,5,:) = (/-1, 0/)
    lat%shell_list(terrace_site,1,6,:) = (/-1, 1/)
    ! Next-Nearest-neigbour (2nd) shell  (d = sqrt(3))
    lat%shell_list(terrace_site,2,1,:) = (/ 1, 1/)
    lat%shell_list(terrace_site,2,2,:) = (/ 2,-1/)
    lat%shell_list(terrace_site,2,3,:) = (/ 1,-2/)
    lat%shell_list(terrace_site,2,4,:) = (/-1,-1/)
    lat%shell_list(terrace_site,2,5,:) = (/-2, 1/)
    lat%shell_list(terrace_site,2,6,:) = (/-1, 2/)
    ! Next-Next-Nearest-neigbour (3rd) shell  (d = 2)
    lat%shell_list(terrace_site,3,1,:) = (/ 0, 2/)
    lat%shell_list(terrace_site,3,2,:) = (/ 2, 0/)
    lat%shell_list(terrace_site,3,3,:) = (/ 2,-2/)
    lat%shell_list(terrace_site,3,4,:) = (/ 0,-2/)
    lat%shell_list(terrace_site,3,5,:) = (/-2, 0/)
    lat%shell_list(terrace_site,3,6,:) = (/-2, 2/)

    ! step

    ! NN list for the hexagonal structure (row,col)
    ! col 3 is a step,      col 4 is a corner
    ! col 2 is a terrace_s, col 5 is a terrace_c

    !  t     ts    s      c     tc
    !  11    12    13***  14**  15
    !
    !     21    22**   23*   24*   25***
    !
    !        31*** 32*   33    34**   35
    !
    !           41**  42*   43*   44   45
    !
    !              51*** 52**  53*** 54    55

    ! Nearest-neigbour (1st) shell)
    lat%shell_list(step_site,1,1,:) = (/-1, 1/) ! step -> corner    (d=2/3)
    lat%shell_list(step_site,1,2,:) = (/ 1, 0/) ! step -> step      (d=1)
    lat%shell_list(step_site,1,3,:) = (/-1, 0/) ! step -> step      (d=1)
    lat%shell_list(step_site,1,4,:) = (/ 1,-1/) ! step -> terrace_s (d=1)
    lat%shell_list(step_site,1,5,:) = (/ 0,-1/) ! step -> terrace_s (d=1)
    ! Next-Nearest-neigbour (2nd) shell
    lat%shell_list(step_site,2,1,:) = (/ 0, 1/) ! step -> corner (d=sqrt(13)/3)
    lat%shell_list(step_site,2,2,:) = (/-2, 1/) ! step -> corner (d=sqrt(13)/3)
    lat%shell_list(step_site,2,3,:) = (/ 2,-1/) ! step -> terrace_s (d=sqrt(3))
    lat%shell_list(step_site,2,4,:) = (/-1,-1/) ! step -> terrace_s (d=sqrt(3))
    lat%shell_list(step_site,2,5,:) = (/ 1,-2/) ! step -> terrace   (d=sqrt(3))
    ! Next-Next-Nearest-neigbour (3rd) shell
    lat%shell_list(step_site,3,1,:) = (/-1, 2/) ! step -> terrace_c (d=4/3)
    lat%shell_list(step_site,3,2,:) = (/ 2, 0/) ! step -> step      (d=2)
    lat%shell_list(step_site,3,3,:) = (/-2, 0/) ! step -> step      (d=2)
    lat%shell_list(step_site,3,4,:) = (/ 2,-2/) ! step -> terrace   (d=2)
    lat%shell_list(step_site,3,5,:) = (/ 0,-2/) ! step -> terrace   (d=2)

    ! corner

    ! NN list for the hexagonal 111-structure (row,col)
    ! col 3 is a corner, col 4 is a terrace_c
    ! col 2 is a step,   col 1 is a terrace_s

    !  ts    s     c      tc    t
    !  11    12    13***  14    15
    !
    !     21    22    23*   24**  25***
    !
    !        31    32**  33    34*   35***
    !
    !           41*** 42*   43*   44**  45
    !
    !              51*** 52**  53*** 54    55

    ! Nearest-neigbour (1st) shell
    lat%shell_list(corner_site,1,1,:) = (/ 1,-1/) ! corner -> step      (d=2/3)
    lat%shell_list(corner_site,1,2,:) = (/ 0, 1/) ! corner -> terrace_c (d=2/3)
    lat%shell_list(corner_site,1,3,:) = (/ 1, 0/) ! corner -> corner    (d=1)
    lat%shell_list(corner_site,1,4,:) = (/-1, 0/) ! corner -> corner    (d=1)

    ! Next-Nearest-neigbour (2nd) shell
    lat%shell_list(corner_site,2,1,:) = (/-1, 1/) ! corner -> terrace_c (d=sqrt(3)-2/3)
    lat%shell_list(corner_site,2,2,:) = (/ 1, 1/) ! corner -> terrace_c (d=sqrt(3)-2/3)
    lat%shell_list(corner_site,2,3,:) = (/ 2,-1/) ! corner -> step      (d=sqrt(13)/3)
    lat%shell_list(corner_site,2,4,:) = (/ 0,-1/) ! corner -> step      (d=sqrt(13)/3)

    ! Next-Next-Nearest-neigbour (3rd) shell
    lat%shell_list(corner_site,3,1,:) = (/ 0, 2/) ! corner -> terrace   (d=sqrt( (sqrt(3)-1/3)^2 + 1/4 ) )
    lat%shell_list(corner_site,3,5,:) = (/-1, 2/) ! corner -> terrace   (d=sqrt( (sqrt(3)-1/3)^2 + 1/4 ) )
    lat%shell_list(corner_site,3,2,:) = (/ 2, 0/) ! corner -> corner    (d=2)
    lat%shell_list(corner_site,3,5,:) = (/-2, 0/) ! corner -> corner    (d=2)
    lat%shell_list(corner_site,3,3,:) = (/ 2,-2/) ! corner -> terrace_s (d=sqrt( (5/3)^2 + 1/4 ) )
    lat%shell_list(corner_site,3,4,:) = (/ 1,-2/) ! corner -> terrace_s (d=sqrt( (5/3)^2 + 1/4 ) )

    ! terrace_s

    ! NN list for the hexagonal structure (row,col)
    ! col 3 is a terrace_s, col 4 is a step
    ! col 2 is a terrace,   col 5 is a corner

    !  t     t     ts     s     c
    !  11    12    13***  14**  15**
    !
    !     21    22**   23*   24*   25**
    !
    !        31*** 32*   33    34*   35***
    !
    !           41**  42*   43*   44**  45
    !
    !              51*** 52**  53*** 54    55

    ! Nearest-neigbour (1st) shell
    lat%shell_list(terrace_s_site,1,1,:) = (/ 0, 1/) ! ts -> step (d=1)
    lat%shell_list(terrace_s_site,1,2,:) = (/ 1, 0/) ! ts -> ts   (d=1)
    lat%shell_list(terrace_s_site,1,3,:) = (/ 1,-1/) ! ts -> t    (d=1)
    lat%shell_list(terrace_s_site,1,4,:) = (/ 0,-1/) ! ts -> t    (d=1)
    lat%shell_list(terrace_s_site,1,5,:) = (/-1, 0/) ! ts -> ts   (d=1)
    lat%shell_list(terrace_s_site,1,6,:) = (/-1, 1/) ! ts -> step (d=1)
    ! Next-Nearest-neigbour (2nd) shell
    lat%shell_list(terrace_s_site,2,1,:) = (/ 2,-1/) ! ts -> t (d=sqrt(3))
    lat%shell_list(terrace_s_site,2,2,:) = (/ 1,-2/) ! ts -> t (d=sqrt(3))
    lat%shell_list(terrace_s_site,2,3,:) = (/-1,-1/) ! ts -> t (d=sqrt(3))
    lat%shell_list(terrace_s_site,2,4,:) = (/ 1, 1/) ! ts -> s (d=sqrt(3))
    lat%shell_list(terrace_s_site,2,5,:) = (/-2, 1/) ! ts -> s (d=sqrt(3))
    lat%shell_list(terrace_s_site,2,6,:) = (/-2, 2/) ! ts -> c (d=sqrt( (5/3)^2 + 1/4 ) )
    lat%shell_list(terrace_s_site,2,7,:) = (/-1, 2/) ! ts -> c (d=sqrt( (5/3)^2 + 1/4 ) )
    ! Next-Next-Nearest-neigbour (3rd) shell
    lat%shell_list(terrace_s_site,3,2,:) = (/ 2, 0/) ! ts -> ts (d=2)
    lat%shell_list(terrace_s_site,3,3,:) = (/ 2,-2/) ! ts -> t  (d=2)
    lat%shell_list(terrace_s_site,3,4,:) = (/ 0,-2/) ! ts -> t  (d=2)
    lat%shell_list(terrace_s_site,3,5,:) = (/-2, 0/) ! ts -> ts (d=2)

    ! terrace_c

    ! NN list for the hexagonal structure (row,col)
    ! col 3 is a terrace_c, col 4 is a terrace
    ! col 2 is a corner,    col 5 is a terrace

    !  s     c     tc     t     t
    !  11    12    13***  14**  15***
    !
    !     21   22**   23*   24*   25**
    !
    !        31    32*   33    34*   35***
    !
    !           41*** 42**  43*   44**  45
    !
    !              51    52    53*** 54    55

    ! Nearest-neigbour (1st) shell
    lat%shell_list(terrace_c_site,1,1,:) = (/ 0,-1/) ! tc -> c    (d=2/3)
    lat%shell_list(terrace_c_site,1,2,:) = (/ 1, 0/) ! tc -> tc   (d=1)
    lat%shell_list(terrace_c_site,1,3,:) = (/-1, 0/) ! tc -> tc   (d=1)
    lat%shell_list(terrace_c_site,1,4,:) = (/ 0, 1/) ! tc -> t    (d=1)
    lat%shell_list(terrace_c_site,1,5,:) = (/-1, 1/) ! tc -> t    (d=1)
    ! Next-Nearest-neigbour (2nd) shell
    lat%shell_list(terrace_c_site,2,1,:) = (/-1,-1/) ! tc -> c    (d = sqrt(3)-2/3)
    lat%shell_list(terrace_c_site,2,2,:) = (/ 1,-1/) ! tc -> c    (d = sqrt(3)-2/3)
    lat%shell_list(terrace_c_site,2,3,:) = (/ 1, 1/) ! tc -> t    (d = sqrt(3))
    lat%shell_list(terrace_c_site,2,4,:) = (/-2, 1/) ! tc -> t    (d = sqrt(3))
    lat%shell_list(terrace_c_site,2,5,:) = (/-1, 2/) ! tc -> t    (d = sqrt(3))
    ! Next-Next-Nearest-neigbour (3rd) shell  (d = 2)
    lat%shell_list(terrace_c_site,3,1,:) = (/ 0, 2/) ! tc -> t    (d=2)
    lat%shell_list(terrace_c_site,3,2,:) = (/ 2, 0/) ! tc -> t    (d=2)
    lat%shell_list(terrace_c_site,3,3,:) = (/-2, 0/) ! tc -> t    (d=2)
    lat%shell_list(terrace_c_site,3,4,:) = (/-2, 2/) ! tc -> t    (d=2)
    lat%shell_list(terrace_c_site,3,5,:) = (/ 1,-2/) ! s  -> c    (d=sqrt(13)/3)

    ! Check the distances to the neighbours
!    print*, sqrt( &
!                 ( lat%shell_list(2,:,1)*cos(pi/3.0_dp) &
!                  +lat%shell_list(2,:,2)               )**2    &
!                +                                                    &
!                 ( lat%shell_list(2,:,1)*sin(pi/3.0_dp))**2    &
!        )
!    stop

    ! List of additional nn directions for hex lattice to scan after a reaction
    allocate(lat%nn_new( lat%n_nn(terrace_site,1), lat%n_nn(terrace_site,1)/2) )
    do m=1,lat%n_nn(terrace_site,1)
    do i=1,lat%n_nn(terrace_site,1)/2
      lat%nn_new(m,i) = modulo( m + i - lat%n_nn(terrace_site,1)/2, lat%n_nn(terrace_site,1) ) + 1
    end do
    end do

    lat%n_max_ads_sites = n_max_ads_sites

    ! Allocate arrays
    allocate(lat%occupations(c_pars%n_rows,c_pars%n_cols))
    allocate(lat%lst(        c_pars%n_rows,c_pars%n_cols))
    allocate(lat%n_ads(c_pars%n_species))
    ! Warning: the allocation assumes 1 ads. per unit cell
    allocate(lat%ads_list(    c_pars%n_rows*c_pars%n_cols))
    ! Initialize arrays
    lat%ads_list     = adsorbate(0,0,0,0)

    lat%lst = terrace_site
    ! Define where steps, corners, terrace_c's and terrace_s's are
    if ( c_pars%step_period > 0) then
      do j=1,lat%n_cols,c_pars%step_period
          lat%lst(:,j)   = step_site
          lat%lst(:,j+1) = corner_site
          lat%lst(:,j+2) = terrace_c_site
          lat%lst(:,j+c_pars%step_period-1) = terrace_s_site
      end do
      n_site_types = n_max_lat_site_types
    else
      n_site_types = 1
    end if

    allocate(lat%avail_ads_sites(c_pars%n_species,n_site_types))

    do i=1,c_pars%n_species
    do j=1,n_site_types
      counter = 0
      do k=1,n_max_ads_sites
        if (e_pars%ads_energy(i,j,k) < e_pars%undefined_energy) then
          counter = counter + 1
        end if
      end do
      allocate(lat%avail_ads_sites(i,j)%list(counter))
      counter = 0
      do k=1,n_max_ads_sites
        if (e_pars%ads_energy(i,j,k) < e_pars%undefined_energy) then
          counter = counter + 1
          lat%avail_ads_sites(i,j)%list(counter) = k
        end if
      end do
    end do
    end do

!!----------------------------------------------------------------------------------
!! debug printout for construction of avail_ads_sites
!print *, 'debug printout for construction of avail_ads_sites mc_lat_class line 191'
!do i=1,c_pars%n_species
!do j=1,n_site_types
!  write(*, '(A, A5, A, A10)',advance='no')' species = ',c_pars%ads_names(i),'site type =',lat_site_names(j)
!  print '(A,A5,1x,A,1x,A,1x,A,1x,A,1x)','  list = ',(ads_site_names(lat%avail_ads_sites(i,j)%list(k)), &
!                    k=1,size(lat%avail_ads_sites(i,j)%list))
!  end do
!end do
!----------------------------------------------------------------------------------

    ! initialize lattice
    lat%n_ini_confs = 0
    call lat%conf_init(c_pars)

  end function mc_lat_init

! ---------------------------------------------------------------------
! Subroutine setting mc_lat structure as defined by (itraj mod nconfs)th record in conf-file
!   call without optional argument
!     - selecting the last conf from conf-file and counts the number of confs in it
!     - producing artificial conf if no conf-file defined
!
! ---------------------------------------------------------------------
  subroutine mc_lat_conf_init(this,c_pars, itraj)

    class(mc_lat), intent(inout) :: this
    type(control_parameters), intent(in) :: c_pars
    integer, optional, intent(in) :: itraj

    integer :: n_rows_in, n_cols_in, step_period_in
    character(len=max_string_length) :: line
    integer :: n_species_in
    character(len=max_string_length) :: ads_names_in(100)
    integer, allocatable :: n_ads_in(:)
    integer :: ios, iconf, i, j, idummy, conf
    integer :: counter, s_counter, current_species, ads_site
    integer :: i_rand, species, row, col

    ! Clean up lattice structure
    this%occupations = 0
    ! Set n_ads as in control file
    this%n_ads = c_pars%n_ads

    ! add adsorbates to the lattice
    if (c_pars%cfg_file_name=='none') then

      ! create random population
      ! WARNING: consider adding a special case for higher coverages
      counter = 0
      do species=1,size(this%n_ads)

        i = 0
        do while (i < this%n_ads(species) )

          ! select a random site
          i_rand = irand(this%n_cols*this%n_rows)
          row = (i_rand-1)/this%n_cols + 1
          col = i_rand - (row - 1)*this%n_cols

          if (this%occupations(row,col) == 0) then
            i = i + 1
            counter = counter + 1
            this%occupations(row,col) = counter
            ! Warning: arbitrary choice for the ads. site (the 1st available)!
            ads_site = this%avail_ads_sites(species,this%lst(row,col))%list(1)
            this%ads_list(counter) = adsorbate(row,col,ads_site,species)
          end if

        end do
      end do

      ! create artificial population
!      counter = 0
!      s_counter = 0
!      current_species = 1
!      loop1: do j=1,this%n_cols
!      do i=1,this%n_rows
!        counter   = counter   + 1
!        s_counter = s_counter + 1
!        if (s_counter>this%n_ads(current_species)) then
!          s_counter = 1
!          current_species = current_species + 1
!        end if
!        this%occupations(i,j) = counter
!        ! Warning: arbitrary choice for the ads. site (the 1st available)!
!        ads_site = this%avail_ads_sites(current_species,this%lst(i,j))%list(1)
!        this%ads_list(counter) = adsorbate(i,j,ads_site,current_species)
!        if (counter == this%n_ads_tot()) exit loop1
!      end do
!      end do loop1

    else
      ! Read in a conf from conf-file

      if (present(itraj)) then
        if (this%n_ini_confs == 0) then
          print*, 'Error: mc_lat_conf_init called with two arguments without prior initialization.'
          stop
        end if
        ! number of conf to read from conf-file
        conf = mod(itraj-1,this%n_ini_confs)+1
      else
        ! allows to read the last conf from conf-file
        conf = huge(1)
      end if

      ! Read info from the configuration file
      call open_for_read(inp_unit,trim(c_pars%cfg_file_name))
      read(inp_unit,'(A)') line
      read(inp_unit,*) n_rows_in, n_cols_in, step_period_in
      ! Check consistency with control parameters
      if (n_rows_in      /= c_pars%n_rows .OR.&
          n_cols_in      /= c_pars%n_cols .OR.&
          step_period_in /= c_pars%step_period) then
        print*, 'Error in mc_lat_conf_init:'
        print*, ' inconsistency in control and configuration files:'
        print*, '  lattice definition'
        stop
      end if

      read(inp_unit,'(A)') line
      call split_string (line, ads_names_in, n_species_in)
      if ( n_species_in /= c_pars%n_species ) then
        print*, 'Error in mc_lat_conf_init:'
        print*, '  inconsistent number of species in control and configuration files.'
        stop
      end if
      if ( any(ads_names_in(1:n_species_in) /= c_pars%ads_names) ) then
        print*, 'Error in mc_lat_conf_init:'
        print*, '  inconsistent adsorbate names in control and configuration files.'
        stop
      end if
      allocate(n_ads_in(n_species_in))

      do iconf=1,conf

        read(inp_unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        read(inp_unit,*) n_ads_in
        if ( any(n_ads_in /= c_pars%n_ads) ) then
          print*, 'Error in mc_lat_conf_init:'
          print*, '  inconsistent number of adsorbates in control and configuration files.'
          stop
        end if
        ! Read the configuration
        do i=1, sum(n_ads_in)
          read(inp_unit, *, iostat=ios) idummy,this%ads_list(i)
          if (ios == -1) then
            print*, 'mc_lat_conf_init:'
            stop " Error in control file: premature end of the final configuration"
          end if
          if (ios /= 0) then
            print*, 'mc_lat_conf_init:'
            print*, " Error reading conf-file: ios = ",ios
            stop
          end if
        end do

      end do
      ! save number of confs in the conf-file
      if (.not.present(itraj)) this%n_ini_confs = iconf - 1

      close(inp_unit)

      ! put adsorbates to lattice
      do i=1,this%n_ads_tot()
        this%occupations(this%ads_list(i)%row, this%ads_list(i)%col) = i
      end do

    end if

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
    integer :: row_old, col_old, list_size, id, site_type

    row_old = this%ads_list(i)%row
    col_old = this%ads_list(i)%col

    row = modulo(row_old + this%shell_list(this%lst(row_old,col_old),1,ihop,1) - 1, this%n_rows) + 1
    col = modulo(col_old + this%shell_list(this%lst(row_old,col_old),1,ihop,2) - 1, this%n_cols) + 1

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
    integer :: row_old, col_old

    row_old = this%ads_list(i)%row
    col_old = this%ads_list(i)%col

    ishell = 1
    if (present(shell)) ishell = shell
    row = modulo(row_old + this%shell_list(this%lst(row_old,col_old),ishell,ihop,1) - 1, this%n_rows) + 1
    col = modulo(col_old + this%shell_list(this%lst(row_old,col_old),ishell,ihop,2) - 1, this%n_cols) + 1

  end subroutine

! ---------------------------------------------------------------------
! Subroutine applying PBCs to calculate distance
! ---------------------------------------------------------------------
  subroutine mc_lat_distance_with_pbc(this, ads1, ads2, r)

    class(mc_lat), intent(in) :: this
    integer, intent(in)   :: ads1, ads2
    real(dp), intent(out) :: r

    integer :: d_row, d_col, n_rows_2, n_cols_2

    n_rows_2 = this%n_rows/2
    n_cols_2 = this%n_cols/2
    ! calculate delta_row and delta_col
    d_row = this%ads_list(ads1)%row - this%ads_list(ads2)%row
    d_col = this%ads_list(ads1)%col - this%ads_list(ads2)%col
    ! apply periodic boundary conditions
    if (d_row > n_rows_2 ) then
      d_row = d_row - this%n_rows
    else if (d_row < -n_rows_2 ) then
      d_row = d_row + this%n_rows
    end if
    if (d_col > n_cols_2 ) then
      d_col = d_col - this%n_cols
    else if (d_col < -n_cols_2 ) then
      d_col = d_col + this%n_cols
    end if
    ! calculate the distance
    r = sqrt(  ( d_col*this%lat_vec_1(1) + d_row*this%lat_vec_2(1) )**2 &
             + ( d_col*this%lat_vec_1(2) + d_row*this%lat_vec_2(2) )**2  )

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
!  subroutine hoshen_kopelman(this, species, cluster_label, largest_label)
!
!  implicit none
!
!  class (mc_lat), intent(in) :: this
!  integer, intent(in) :: species
!  integer, dimension(:,:), intent(out) :: cluster_label
!  integer, intent(out) :: largest_label
!
!  type nn_position
!    integer :: shell,m
!  end type nn_position
!
!  integer :: i, m, n, row, col, row_new,col_new, shell
!  integer :: itemp, ip, jp, i_ads, i_ads_nn, n_neighbors
!
!  integer, dimension(n_shells) :: nnn2
!  integer, dimension(this%n_ads(species)) :: labels
!  integer, dimension(n_shells, maxval(this%n_nn)/2) :: scanned_nn_occs, row_nn, col_nn
!  type(nn_position), dimension(sum(this%n_nn)/2) :: scanned_nn_list
!
!
!  nnn2 = this%n_nn(terrace_site,:)/2
!  scanned_nn_occs = 0
!  cluster_label   = 0
!  do i=1,this%n_ads(species)
!      labels(i) = i
!  end do
!
!  largest_label = 0
!  do row=1,this%n_rows
!  do col=1,this%n_cols
!
!    i_ads = this%occupations(row,col)
!
!    if ( i_ads > 0 ) then
!      if ( this%ads_list(i_ads)%id == species ) then
!
!        do shell=1,n_shells
!        do m=1,nnn2(shell)
!
!          ! Take previously scanned neighbours
!          call this%neighbor(i_ads,nnn2(shell)+m,row_nn(shell,m),col_nn(shell,m), shell)
!          i_ads_nn = this%occupations(row_nn(shell,m), col_nn(shell,m))
!
!          scanned_nn_occs(shell,m) = 0
!          if ( i_ads_nn > 0 ) then
!            if ( this%ads_list(i_ads_nn)%id == species ) scanned_nn_occs(shell,m) = 1
!          end if
!
!        end do
!        end do
!
!        ! Switch off periodic boundary conditions
!        ! Warning: works only for hexagonal lattice, because of Corona virus
!
!        ! 1st shell
!        if (row==1) scanned_nn_occs(1,2:3) = 0
!        if (col==1) scanned_nn_occs(1,1) = 0
!        if (col==this%n_cols) scanned_nn_occs(1,3) = 0
!
!        ! 2nd shell
!        if (row==1) scanned_nn_occs(2,:) = 0
!        if (row==2) scanned_nn_occs(2,2) = 0
!        if (col==1) scanned_nn_occs(2,1) = 0
!        if (col==this%n_cols) scanned_nn_occs(2,2:3) = 0
!        if (col==this%n_cols-1) scanned_nn_occs(2,3) = 0
!
!        ! 3rd shell
!        if (row==1 .or. row==2) scanned_nn_occs(3,2:3) = 0
!        if (col==1 .or. col==2) scanned_nn_occs(3,1) = 0
!        if (col==this%n_cols .or. col==this%n_cols-1) scanned_nn_occs(3,3) = 0
!
!        ! create list of neighbors
!        n_neighbors = 0
!        scanned_nn_list = nn_position(0,0)
!        do shell=1,n_shells
!        do m=1,nnn2(shell)
!          if ( scanned_nn_occs(shell,m) == 1 ) then
!            n_neighbors = n_neighbors + 1
!            scanned_nn_list(n_neighbors)%shell = shell
!            scanned_nn_list(n_neighbors)%m     = m
!          end if
!        end do
!        end do
!
!       select case (n_neighbors)
!
!          case (0)
!!            print*
!!            print*,"Case 0"
!            largest_label = largest_label + 1
!            cluster_label(row,col) = largest_label
!
!          case (1)
!!            print*
!!            print*,"Case 1"
!            shell = scanned_nn_list(1)%shell
!            m     = scanned_nn_list(1)%m
!            cluster_label(row,col) = lfind( cluster_label(row_nn(shell,m),col_nn(shell,m)), labels)
!
!          case default
!            shell = scanned_nn_list(1)%shell
!            m     = scanned_nn_list(1)%m
!            itemp = cluster_label(row_nn(shell,m),col_nn(shell,m))
!            do i=2,n_neighbors
!              shell = scanned_nn_list(i)%shell
!              n     = scanned_nn_list(i)%m
!              call lunion(itemp,cluster_label(row_nn(shell,n),col_nn(shell,n)),labels)
!            end do
!            cluster_label(row,col) = lfind(itemp, labels)
!
!        end select
!
!      end if
!    end if
!
!  end do
!  end do
!
!!  print*
!!  call this%print_ocs
!!  print*
!!  print*,'cluster labels before applying PBC:'
!!  write(*,'(20i4)') transpose(cluster_label)
!!  print*,'cluster labels list:'
!!  do i=1,largest_label
!!      print('(i3,a5,i3)'), i, ' --> ' , labels(i)
!!  end do
!
!  ! Apply PBC
!  do row=1,this%n_rows
!  do col=1,2
!    ip = cluster_label(row,col)
!    if (ip > 0) then
!      do shell=1,n_shells
!      do m=1,this%n_nn(terrace_site,shell)
!        call this%neighbor(this%occupations(row,col),m,row_new,col_new,shell)
!        jp = cluster_label(row_new,col_new)
!        if (jp > 0) call lunion(ip,jp,labels)
!      end do
!      end do
!    end if
!  end do
!  end do
!
!  do col=3,this%n_cols-2
!  do row=1,2
!    ip = cluster_label(row,col)
!    if (ip > 0) then
!      do shell=1, n_shells
!      do m=1,this%n_nn(terrace_site,shell)
!        call this%neighbor(this%occupations(row,col),m,row_new,col_new,shell)
!        jp = cluster_label(row_new,col_new)
!        if (jp > 0) call lunion(ip,jp,labels)
!      end do
!      end do
!    end if
!  end do
!  end do
!
!!  print*
!!  print*,'cluster labels list after PBC:'
!!  do i=1,largest_label
!!      print('(i3,a5,i3)'), i, ' --> ' , labels(i)
!!  end do
!
!  ! Going down to roots
!  do row=1,this%n_rows
!  do col=1,this%n_cols
!    if (cluster_label(row,col) > 0) &
!       cluster_label(row,col) = lfind(cluster_label(row,col), labels)
!  end do
!  end do
!
!!  print*
!!  print*,'cluster labels after PBC and rooting:'
!!  write(*,'(20i4)') transpose(cluster_label)
!!  pause
!
!  contains
!
!    integer function lfind(x, labels)
!
!      integer, intent(inout) :: x
!      integer, dimension(:), intent(inout) :: labels
!
!      integer :: y, z
!
!      y = x
!        do while (labels(y) /= y)
!            y = labels(y)
!        end do
!
!    !    do while (labels(x) /= x)
!    !        z = labels(x)
!    !        labels(x) = y
!    !        x = z
!    !    end do
!    !
!        lfind = y
!
!    end function
!
!    subroutine lunion(x, y, labels)
!
!      integer, intent(inout) :: x, y
!      integer, dimension(:), intent(inout) :: labels
!
!        if (x > y) then
!            labels( lfind(x, labels) ) = lfind(y, labels)
!        else if ( x < y) then
!            labels( lfind(y, labels) ) = lfind(x, labels)
!        endif
!
!    end subroutine
!
!  end subroutine
!
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
! Subroutine calculating a histogram for the radial distribution function
! ---------------------------------------------------------------------
  subroutine rdf_hist(this, cpars, hist)

  implicit none

  class (mc_lat), intent(in) :: this
  class(control_parameters), intent(in) :: cpars
  integer, dimension(:,:,:),  intent(inout) :: hist

  integer :: ads1, spec1
  integer :: ads2, spec2
  integer :: n_ads_total, bin, n_bins
  real(dp):: r, inv_bin_size

  n_ads_total = this%n_ads_tot()
  inv_bin_size = 1.0_dp/cpars%rdf_bin_size
  n_bins = cpars%rdf_n_bins

  do ads1 = 1, n_ads_total-1
  do ads2 = ads1+1, n_ads_total

    spec1 = this%ads_list(ads1)%id
    spec2 = this%ads_list(ads2)%id
    call this%distance(ads1, ads2, r)
    bin = floor(r*inv_bin_size + 0.5_dp) + 1
    if (bin <= n_bins) then
      hist(spec1,spec2,bin) = hist(spec1,spec2,bin) + 1
      hist(spec2,spec1,bin) = hist(spec2,spec1,bin) + 1
    end if

  end do
  end do

  end subroutine

end module mc_lat_class




