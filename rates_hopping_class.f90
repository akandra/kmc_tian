module rates_hopping_class

  use constants
  use control_parameters_class
  use mc_lat_class
  use energy_parameters_class
  use energy_mod
  use open_file
  use utilities
  use temperature_laws

  implicit none

  private
  public    :: hopping_init, hopping_type

  type :: v_list_dp
    real(dp), dimension(:), allocatable :: list ! n_avail_ads_sites
  end type

  type :: int_law_pars
    integer                              :: id   ! law id
    real(dp), dimension(n_max_rcic_pars) :: pars ! parameter list
  end type

  type :: hopping_type

    logical :: is_defined           = .false.

    ! Hopping Rates
    !                          n_adsorbate              -> which particle
    !                          .  n_neighbor            -> where to
    !                          .  .
    type(v_list_dp), dimension(:, :), allocatable :: rates


          !                    (   n_species                       -> which species
          !                        .  n_site_type                  -> where from
          !                        .  .  n_adsorption_sites
          !                        .  .  .  n_site_type            -> where to
          !                        .  .  .  .  n_adsorption_sites )
          !                        .  .  .  .  .
    real(dp),            dimension(:, :, :, :, :), allocatable :: process
    type(int_law_pars),  dimension(:, :, :, :, :), allocatable :: rate_corr_pars


          !                    (   n_species                       -> which species
          !                        .  n_site_type                  -> where from
          !                        .  .  n_adsorption_sites
          !                        .  .  .  n_adsorption_sites )   -> where to
          !                        .  .  .  .
    real(dp),            dimension(:, :, :, :), allocatable :: process_intra
    type(int_law_pars),  dimension(:, :, :, :), allocatable :: rate_corr_pars_intra

  contains
    procedure :: construct
    procedure :: print

  end type

contains
!------------------------------------------------------------------------------
  function hopping_init(c_pars, lat, e_pars)
!------------------------------------------------------------------------------
    type(hopping_type) hopping_init

    type(control_parameters), intent(in)    :: c_pars
    type(mc_lat)            , intent(in)    :: lat
    type(energy_parameters) , intent(in)    :: e_pars

    integer :: i, ios, nwords, line_number, i1, i2, i3, i4, m

    integer :: species, st1, st2, ast1, ast2, col_st1(2)
    logical :: e_defined1, e_defined2, r_defined, undefined_rate, undefined_energy

    character(len=max_string_length)                :: buffer
    character(len=max_string_length)                :: words(100)
    character(len=len(trim(c_pars%rate_file_name))) :: file_name

    character(len=10)     :: current_species_name
    integer               :: current_species_id
    integer               :: rct_law_id, rct_law_id_glob
    integer               :: rcic_law_id, rcic_law_id_glob

    integer               :: parse_state
    integer, parameter    :: parse_state_ignore   = -1
    integer, parameter    :: parse_state_default  =  0
    integer, parameter    :: parse_state_hopping  =  1

    logical :: rct_law_defined  = .false.
    logical :: rcic_law_defined = .false.

    real(dp), dimension(n_max_rct_pars )::  rct_pars,  rct_pars_glob
    real(dp), dimension(n_max_rcic_pars):: rcic_pars, rcic_pars_glob

    integer,  parameter   :: default_int = 0
    real(dp), parameter   :: default_rate  = -1.0_dp

    integer :: max_avail_ads_sites

    character(len=2) :: check_str
    logical :: is_same_lst
    character(len=100) :: status_str


    ! determine maximum number of available ads. sites
    max_avail_ads_sites = 1

    do i=1,c_pars%n_species
    do m=1,size(lat%avail_ads_sites(i,:))

      i1 = size(lat%avail_ads_sites(i,m)%list)
      if (max_avail_ads_sites < i1) max_avail_ads_sites = i1

    end do
    end do

    ! Allocate and initialize rates array
    allocate( hopping_init%rates(lat%n_rows*lat%n_cols,lat%n_max_nn+1) ) ! '+1' accounts for intra-site hops
    do i=1,lat%n_rows*lat%n_cols
    do m=1,lat%n_max_nn+1
      allocate( hopping_init%rates(i,m)%list(max_avail_ads_sites) )
      hopping_init%rates(i,m)%list = 0.0_dp
    end do
    end do

    allocate(hopping_init%process( c_pars%n_species,&
                                 n_max_lat_site_types, n_max_ads_sites,&
                                 n_max_lat_site_types, n_max_ads_sites) )
    allocate(hopping_init%process_intra( c_pars%n_species,&
                                 n_max_lat_site_types, n_max_ads_sites, n_max_ads_sites) )

    hopping_init%process       = default_rate
    hopping_init%process_intra = default_rate

    allocate(hopping_init%rate_corr_pars( c_pars%n_species,&
                                 n_max_lat_site_types, n_max_ads_sites,&
                                 n_max_lat_site_types, n_max_ads_sites) )
    allocate(hopping_init%rate_corr_pars_intra( c_pars%n_species,&
                                 n_max_lat_site_types, n_max_ads_sites, n_max_ads_sites) )

    hopping_init%rate_corr_pars       = int_law_pars(default_int,default_rate)
    hopping_init%rate_corr_pars_intra = int_law_pars(default_int,default_rate)

    !  read rate definitions from the input file
    file_name = c_pars%rate_file_name
    call open_for_read(inp_unit, file_name )

    ios = 0
    parse_state = parse_state_default
    line_number = 0
    undefined_energy = .false.


    !---------------------------------------------------------------------------
    ! loop over all lines of rate input file
    !---------------------------------------------------------------------------

    do while (ios == 0)

      read(inp_unit, '(A)', iostat=ios) buffer

      line_number = line_number + 1
      ! ios = 0: valid record read
      ! ios < 0: end of record condition encountered or end of file condition detected
      ! ios > 0: an error is detected
      ! print *, 'line: ', line_number, ' ios=', ios, ' buffer=', trim(buffer), ' length=', len(trim(buffer))
      if (ios < 0) buffer=''  ! treat end of file like blank line

      ! Split an input string
      words = ''
      call split_string(buffer, words, nwords)

      ! skip comments but preserve blank lines (signify end)
      if (nwords == 0 .and. buffer /= '') cycle

      select case (parse_state)

        case(parse_state_default)
          ! in parse state default:
          !    word 'hopping' to mark beginning of a hopping section
          !    ignore anything else until hopping section begins

          if (words(1) == 'hopping') then
            hopping_init%is_defined = .true.
            parse_state = parse_state_hopping
            rct_law_defined  = .false.  ! reset necessary to allow multiple hopping sections
            rcic_law_defined = .false.  ! reset necessary to allow multiple hopping sections

            if (nwords == 2) then
              read(words(2),'(A)') current_species_name
              current_species_id = get_index(current_species_name, c_pars%ads_names)
              if (current_species_id == 0) call error_message(file_name, line_number, buffer, &
                                                              "unknown species in hopping section definition")
            else
              call error_message(file_name, line_number, buffer, &
                         "hopping key must have 1 parameter -- species")
            end if !nwords == 2

          end if !words(1)=='hopping'


        case(parse_state_hopping)
          ! process:
          !    temperature law records,
          !    interaction law records,
          !    to-from records,
          !    section end


          if (words(1) == section_end) then
            parse_state = parse_state_default

          elseif (words(1) =='temperature_law') then
            rct_law_id_glob = get_index(words(2), rct_law_names)
            if (rct_law_id_glob == 0) then
              call error_message(file_name, line_number, buffer,&
                                 "invalid temperature law statement")
            else
              rct_law_defined = .true.
              select case (rct_law_id_glob)
                case (Arrhenius_id)
                  if (nwords/=4) call error_message(file_name, line_number, buffer,&
                                                    "Arrhenius must have 2 parameters")
                  read(words(3),*) rct_pars_glob(1)
                  read(words(4),*) rct_pars_glob(2)
                case (extArrhenius_id)
                  if (nwords/=5) call error_message(file_name, line_number, buffer,&
                                                    "extArrhenius must have 3 parameters")
                  read(words(3),*) rct_pars_glob(1)
                  read(words(4),*) rct_pars_glob(2)
                  read(words(5),*) rct_pars_glob(3)
              end select
            endif

          elseif (words(1) =='interaction_law') then
            rcic_law_id_glob = get_index(words(2), rcic_law_names)
            if (rcic_law_id_glob == 0) then
              call error_message(file_name, line_number, buffer,&
                                 "invalid interaction law statement")
            else
              rcic_law_defined = .true.
              select case (rcic_law_id_glob)
                case (rcic_linear_id)
                  if (nwords/=4) call error_message(file_name, line_number, buffer,&
                                                    "linear interaction law must have 2 parameters")
                  read(words(3),*) rcic_pars_glob(1)
                  read(words(4),*) rcic_pars_glob(2)
              end select
            endif

          else
            ! ---------------------------------------------------------
            ! check if we have a valid from-to rate record
            ! ---------------------------------------------------------
            if (nwords < 4) &
              call error_message(file_name, line_number, buffer, &
                                 "invalid to-from rate record in the hopping section")
            i1 = get_index(words(1),lat_site_names)
            i2 = get_index(words(2),ads_site_names)
            i3 = get_index(words(3),lat_site_names)
            i4 = get_index(words(4),ads_site_names)
            is_same_lst = (words(3) == same_lst_mark)
            if (is_same_lst) i3 = i1

            ! check for invalid site name (invalid lst or ast in either from or to site)
            if ( i1==0 .or. i2==0 .or. i3==0 .or. i4==0) then
              print *, 'parse state: ', parse_state
              call error_message(file_name, line_number, buffer, &
                                 "wrong site name in the hopping section")
            end if

            ! ---------------------------------------------------------
            ! check for duplicate entry
            ! ---------------------------------------------------------
            if (is_same_lst) then
              if (hopping_init%process_intra(current_species_id,i1,i2,i4 ) /= default_rate)&
                call error_message(file_name, line_number, buffer, "duplicated entry (check symmetry duplicates)")
            else
              if (hopping_init%process(current_species_id,i1,i2,i3,i4 ) /= default_rate)&
              call error_message(file_name, line_number, buffer, "duplicated entry (check symmetry duplicates)")
            end if

            ! check energy is defined for initial and final site_type and ads_site
            if( e_pars%ads_energy(current_species_id, i1, i2) == e_pars%undefined_energy .or. &
                e_pars%ads_energy(current_species_id, i3, i4) == e_pars%undefined_energy ) then
              call error_message(file_name, line_number, buffer, &
                               "rate defined for site with undefined adsorption energy", &
                               stop=.false., warning=.false.)

              undefined_energy = .true.
            end if

            ! we have a valid rate record. Process it
            if (nwords == 4) then
              if (rct_law_defined .and. rcic_law_defined) then
                rct_law_id  = rct_law_id_glob
                rct_pars    = rct_pars_glob
                rcic_law_id = rcic_law_id_glob
                rcic_pars   = rcic_pars_glob
              else
                call error_message(file_name, line_number, buffer, &
                                 "temperature and/or interaction laws are not defined")
              end if
            else
              do i=5,nwords
                ! find rct law
                rct_law_id = get_index(words(i), rct_law_names)
                if (rct_law_id /= 0) then
                  select case (rct_law_id)
                    case (Arrhenius_id)
                      if (nwords<i+2) call error_message(file_name, line_number, buffer,&
                                                        "Arrhenius must have 2 parameters")
                      read(words(i+1),*) rct_pars(1)
                      read(words(i+2),*) rct_pars(2)
                    case (extArrhenius_id)
                      if (nwords<i+3) call error_message(file_name, line_number, buffer,&
                                                  "extArrhenius must have 3 parameters")
                      read(words(i+1),*) rct_pars(1)
                      read(words(i+2),*) rct_pars(2)
                      read(words(i+3),*) rct_pars(3)
                    case default
                      call error_message(file_name, line_number, buffer, "This should not happen! Check the code!")
                  end select
                end if
                ! find rcic law
                rcic_law_id = get_index(words(i), rcic_law_names)
                if (rcic_law_id /= 0) then
                  select case (rcic_law_id)
                    case (rcic_linear_id)
                      if (nwords<i+2) call error_message(file_name, line_number, buffer,&
                                                        "linear interaction law must have 2 parameters")
                      read(words(i+1),*) rcic_pars(1)
                      read(words(i+2),*) rcic_pars(2)
                    case default
                      call error_message(file_name, line_number, buffer, "This should not happen! Check the code!")
                  end select
                end if
              end do
            end if

            ! Check if rct and rcic laws are set
            if (rct_law_id == 0) &
              call error_message(file_name, line_number, buffer, "invalid temperature law parameters")
            if (rcic_law_id == 0) &
              call error_message(file_name, line_number, buffer, "invalid interaction law parameters")
            ! Set rate constants for hopping channels
            if (is_same_lst) then
              hopping_init%process_intra(current_species_id,i1,i2,i4 ) = &
                          rct_law(rct_law_id, c_pars%temperature, rct_pars)
              ! symmetrize
              hopping_init%process_intra(current_species_id,i1,i4,i2 ) = &
              hopping_init%process_intra(current_species_id,i1,i2,i4 )
            else
              hopping_init%process(current_species_id,i1,i2,i3,i4 ) = &
                          rct_law(rct_law_id, c_pars%temperature, rct_pars)
              ! symmetrize
              hopping_init%process(current_species_id,i3,i4,i1,i2 ) = &
              hopping_init%process(current_species_id,i1,i2,i3,i4 )
            end if

          endif ! words(1) == section_end

      end select
    end do ! while ios=0

    close(inp_unit)

    if (parse_state /= parse_state_default) then
      !write(*, '(A)') ' hopping: error, incomplete hopping section'
      print *, 'parse state: ', parse_state
      !stop 9961
      call error_message(file_name, line_number, buffer, &
           'hopping: error, incomplete hopping section')
    endif

    if (undefined_energy) then
      write(*, '(A)') ' hopping: error, rates defined for sites with undefined energies'
      stop 9962

    else
      write(*, '(A)') ' hopping: passed check that energies are defined for all rates'

    end if

  ! ----------------------------------------------------------------------------
  ! Hopping rates report
  ! ----------------------------------------------------------------------------

  write(*,'(A)') '  Hopping Rates Report.'
  write(*,'(A)') '  --------------------------'
  do species=1,c_pars%n_species
  do st1       = 1, n_max_lat_site_types

    col_st1 = get_indices(st1,lat%lst(1,:))
    if (col_st1(1) > 0) then

      do i1=1,size(lat%avail_ads_sites(species,st1)%list)
        ast1 = lat%avail_ads_sites(species,st1)%list(i1)

          write(*,'(2X,A11,A11,A4)') c_pars%ads_names(species),lat_site_names(st1),ads_site_names(ast1)

          ! Intra-site hops
          do i2=1,size(lat%avail_ads_sites(species,st1)%list)
            ast2 = lat%avail_ads_sites(species,st1)%list(i2)
            if (ast2 /= ast1) then
              check_str = '  '
              if (hopping_init%process_intra(species,st1,ast1,ast2) /= default_rate) check_str = check_mark
              write(*,'(5X,A2,X,A11,A11,A4)') check_str, c_pars%ads_names(species), same_lst_mark, ads_site_names(ast2)
            end if

          end do

          ! Inter-site hops
          do m = -1,1

            st2 = lat%lst(1, modulo(col_st1(1) + m - 1, lat%n_cols) + 1)

            do i2=1,size(lat%avail_ads_sites(species,st2)%list)
              ast2 = lat%avail_ads_sites(species,st2)%list(i2)
              check_str = '  '
              if (hopping_init%process(species,st1,ast1,st2,ast2) /= default_rate) check_str = check_mark
              write(*,'(5X,A2,X,A11,A11,A4)') check_str, c_pars%ads_names(species), lat_site_names(st2), ads_site_names(ast2)
            end do

          end do
          write(*,*) ''

      end do

    end if

  end do
  end do



  ! Replace default rates with zeros to escape negative rates

    do species=1,c_pars%n_species
    do st1=1, n_max_lat_site_types
    do st2=1, n_max_lat_site_types
    do i1=1,size(lat%avail_ads_sites(species,st1)%list)
    do i2=1,size(lat%avail_ads_sites(species,st2)%list)

      ast1 = lat%avail_ads_sites(species,st1)%list(i1)
      ast2 = lat%avail_ads_sites(species,st2)%list(i2)

      if (hopping_init%process(species,st1,ast1,st2,ast2) == default_rate)&
        hopping_init%process(species,st1,ast1,st2,ast2) = 0.0_dp

      if (hopping_init%process_intra(species,st1,ast1,ast2) == default_rate)&
        hopping_init%process_intra(species,st1,ast1,ast2) = 0.0_dp

    end do
    end do
    end do
    end do
    end do

  end function hopping_init

!-----------------------------------------------------------------------------
  subroutine construct(this, ads, lat, e_pars, beta)
!-----------------------------------------------------------------------------
    class(hopping_type), intent(inout) :: this
    integer, intent(in) :: ads
    class(mc_lat), intent(inout) :: lat
    class(energy_parameters), intent(in) :: e_pars
    real(dp), intent(in) :: beta

    integer :: id, m, iads
    integer :: row_old, col_old, lst_old, ast_old
    integer :: row_new, col_new, lst_new, ast_new
    real(dp) :: energy_old, energy_new

    ! energy for particle ads in its old position
    energy_old = energy(ads, lat, e_pars)
    ! Save the old configuration
    row_old = lat%ads_list(ads)%row
    col_old = lat%ads_list(ads)%col
    lst_old = lat%lst(row_old,col_old)
    ast_old = lat%ads_list(ads)%ast
    id      = lat%ads_list(ads)%id

    ! Construct rates for hops to neighbors

    ! Delete particle ads from the old position
    ! we do it here since we never work with occupations inside the following loop
    lat%occupations(row_old,col_old) = 0

    ! Loop over possible new positions of particle ads
    do m=1,lat%n_nn(lst_old, 1)

      ! Get position and site type of neighbour m
      call lat%neighbor(ads, m, row_new, col_new)
      lst_new  = lat%lst(row_new, col_new)

      ! Check if the cell is free
      if (lat%occupations(row_new, col_new) > 0) then

        this%rates(ads,m)%list = 0.0d0

      else

        ! Put particle ads to site m
        lat%ads_list(ads)%row = row_new
        lat%ads_list(ads)%col = col_new
        lat%occupations(row_new,col_new) = ads

        ! Loop over adsorption sites
        do iads = 1, size(lat%avail_ads_sites(id,lst_new)%list)

          ! Move particle ads to adsorption site list(iads)
          ast_new = lat%avail_ads_sites(id,lst_new)%list(iads)
          lat%ads_list(ads)%ast = ast_new

          ! Calculate energy of ads in new position
          energy_new = energy(ads, lat, e_pars)

          if (debug(2)) then
            print*, '-----------construct debug(2)'
            print*, 'pos new 1 ', lat%ads_list(1)%row, lat%ads_list(1)%col
            print*, 'pos new 2 ', lat%ads_list(2)%row, lat%ads_list(2)%col
            print *, 'ads ', ads, 'm ', m,  'E_old = ', energy_old, 'E_new = ', energy_new
          end if

          ! Apply detailed balance when
          ! energy in the old position < energy in the new position
          if (energy_old < energy_new) then
              this%rates(ads,m)%list(iads) = this%process(id, lst_old, ast_old, lst_new, ast_new)&
                  *exp( -beta*(energy_new - energy_old) )
!                  if (row_old == 1 .and. lat%ads_list(ads)%row == 6.and. lat%ads_list(2)%row == 5  ) then
!                    print*, '-----------------hopping'
!                    print*, 'id ',id, ' old site ',lst_old,' old ads. site ', ast_old
!                    print*, ' new site ',lst_new,' new ads. site ', ast_new
!                    print*, 'pos_old ', row_old, col_old
!                    print*, 'pos 1 ', lat%ads_list(1)%row, lat%ads_list(1)%col
!                    print*, 'pos 2 ', lat%ads_list(2)%row, lat%ads_list(2)%col
!                    print*, 'rate ',this%process(id, lst_old, ast_old, lst_new, ast_new)
!                    print*, 'E_old ', energy_old, 'E_new ', energy_new, 'db factor ', exp( -beta*(energy_new - energy_old) )
!                    print*, 'ads ', ads, 'm ', m, 'iads', iads, 'rate constant ', this%rates(ads,m)%list(iads)
!                  end if
          else
              this%rates(ads,m)%list(iads) = this%process(id, lst_old, ast_old, lst_new, ast_new)
!                print*, id,lst_old, ast_old, lst_new, ast_new,this%process(id, lst_old, ast_old, lst_new, ast_new)
          end if

!            print '(A,i4,A,i4,A,i4)', 'ads ',ads,' neighbor ',m,' ads. site ',lat%ads_list(ads)%ast
!            print '(A,e18.4,A,e18.4)',' E_old = ', energy_old, ' E_new = ',energy_new
!            print *,' r_hop = ',this%process(id, lst_old, ast_old, lst_new, ast_new),&
!                                      ' rate = ', this%rates(ads,m)%list(iads)
!            write(*,*) 'pause'
!            read(*,*)

        end do ! iads

        ! Return particle ads to the old position
        lat%ads_list(ads)%row = row_old
        lat%ads_list(ads)%col = col_old
        lat%ads_list(ads)%ast = ast_old
        lat%occupations(row_new,col_new) = 0

      end if ! occupations

    end do ! m

    lat%occupations(row_old,col_old) = ads


    ! Construct rates for intra-site hops

    m = lat%n_nn(lst_old, 1) + 1

    ! Loop over adsorption site types
    do iads = 1, size(lat%avail_ads_sites(id,lst_old)%list)

      ast_new = lat%avail_ads_sites(id,lst_old)%list(iads)

      if ( ast_new == ast_old) then ! exclude a hop into the same ast

        this%rates(ads,m)%list(iads) = 0

      else

        ! Move particle ads to adsorption site list(iads)
        lat%ads_list(ads)%ast = ast_new

        ! Calculate energy of ads in new position
        energy_new = energy(ads, lat, e_pars)

        if (debug(2)) then
          print*, '-----------construct debug(2)'
          print*, 'pos new 1 ', lat%ads_list(1)%row, lat%ads_list(1)%col
          print*, 'pos new 2 ', lat%ads_list(2)%row, lat%ads_list(2)%col
          print *, 'ads ', ads, 'm ', m,  'E_old = ', energy_old, 'E_new = ', energy_new
        end if

        ! Apply detailed balance when
        ! energy in the old position < energy in the new position
        if (energy_old < energy_new) then
            this%rates(ads,m)%list(iads) = this%process_intra(id, lst_old, ast_old, ast_new)&
                *exp( -beta*(energy_new - energy_old) )
        else
            this%rates(ads,m)%list(iads) = this%process_intra(id, lst_old, ast_old, ast_new)
        end if

      end if

    end do ! iads

    ! Return particle ads to the old position
    lat%ads_list(ads)%ast = ast_old


  end subroutine construct

!------------------------------------------------------------------------------
  subroutine print(this, c_pars)
!------------------------------------------------------------------------------
    class(hopping_type), intent(in) :: this

    class(control_parameters), intent(in) :: c_pars

    integer :: i, i1, i2, i3, i4

    print*, 'Hopping Rates:'
    do i=1,size(this%process,1)
      print '(A)',' ---------------------------'
      print '( A,A)',' species: ', c_pars%ads_names(i)
      print '(A)', ' ---------------------------'
      do i1=1,n_max_lat_site_types
      do i2=1,n_max_ads_sites
      do i3=1,n_max_lat_site_types
      do i4=1,n_max_ads_sites
        if (this%process(i,i1,i2,i3,i4)< 0.0_dp) then
          cycle
        else
          write(*,'(1x,A,A,2X,A,A,6(1pe11.2))') &
              lat_site_names(i1), ads_site_names(i2), &
              lat_site_names(i3), ads_site_names(i4), &
              this%process(i,i1,i2,i3,i4)
        end if
      end do
      end do
      end do
      end do
    end do
    print*

  end subroutine print

end module rates_hopping_class
