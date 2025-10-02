module rates_hopping_class

  !---------------------------------------------------------------------------
  !  Module for hopping  r --> p
  !---------------------------------------------------------------------------
  !  * controlled by the hopping keyword section of the .rates file
  !  * r and p are the same species, since it is a hop
  !  * direction of hop is described by a vector (drow, dcol) and 
  !    in general is not restricted
  !---------------------------------------------------------------------------

  use constants
  use control_parameters_class
  use mc_lat_class
  use energy_parameters_class
  use energy_mod
  use open_file
  use utilities
  use rate_constant_laws

  implicit none

  private
  public    :: hopping_init, hopping_type

  type rate_info_hopping
    integer(1) :: proc
    integer(1) :: m
    real(dp)   :: rate
  end type

  type :: v_list_rate_info
    type(rate_info_hopping), dimension(:), allocatable :: list
    integer(1) :: n_channels
  end type

  type :: hopping_def
    integer  :: r         ! reactant's id
    integer  :: r_lst     ! reactant's lst (lattice site type)
    integer  :: r_ast     ! reactant's ast (adsorption site type)
    integer  :: p_lst     ! lst after hop
    integer  :: p_ast     ! ast after hop
    integer  :: p_vec(2)  ! hop vector
    real(dp) :: rate      ! hopping rate
    type(int_law_pars) :: rcic
  end type


  type :: hopping_type


    logical :: is_defined = .false.


    !---------------------------------------------------------------------------
    ! Hopping paths structure
    !---------------------------------------------------------------------------
    ! Gives the rate of hopping for each set of where-from (lst,ast)
    ! and where-to (vector, ast) given in the hopping section of the
    ! .rates input file
    !---------------------------------------------------------------------------
    type(hopping_def), dimension(:), allocatable :: paths

    !---------------------------------------------------------------------------
    ! Hopping channel structure
    !---------------------------------------------------------------------------
    ! Gives a list of hopping information for a given adsorbate
    !   * list elements give the hopping channel which specifies the rate,
    !     where-from lst and ast and the where-to vector and ast
    !---------------------------------------------------------------------------
    !
    !                                 adsorbate      -> which particle
    !                                 .
    type(v_list_rate_info), dimension(:), allocatable :: rate_info

    ! Number of hopping paths
    integer :: n_paths

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

    integer :: i, ios, ntokens, line_number, i1, i2, i3, i4, m, j, n1, n2
    integer :: d1, d2, n_d1, n_d2

    integer :: species, st1, st2, ast1, ast2, col_st1(2)
    logical :: undefined_energy

    character(len=max_string_length)                :: buffer
    character(len=max_string_length)                :: tokens(100)
    character(len=len(trim(c_pars%rate_file_name))) :: file_name

    character(len=10)     :: current_species_name
    integer               :: current_species_id
    integer               :: rct_law_id, rct_law_id_glob
    integer               :: rcic_law_id, rcic_law_id_glob

    integer               :: parse_state
    integer, parameter    :: parse_state_ignore   = -1
    integer, parameter    :: parse_state_default  =  0
    integer, parameter    :: parse_state_hopping  =  hopping_id

    logical :: rct_law_defined  = .false.
    logical :: rcic_law_defined = .false.
    logical :: duplicate_error  = .false.

    real(dp), dimension(n_max_rct_pars )::  rct_pars,  rct_pars_glob
    real(dp), dimension(n_max_rcic_pars):: rcic_pars, rcic_pars_glob

    integer,  parameter   :: default_int = 0
    real(dp), parameter   :: default_rate  = -1.0_dp

    integer :: max_avail_ads_sites, pass
    integer :: n_hopping_paths
    real(dp) :: delta_eps

    character(len=2) :: check_str
    character(len=max_string_length) :: stemp

    ! determine maximum number of available ads. sites
    max_avail_ads_sites = 1

    do i=1,c_pars%n_species
    do m=1,size(lat%avail_ads_sites(i,:))

      i1 = size(lat%avail_ads_sites(i,m)%list)
      if (max_avail_ads_sites < i1) max_avail_ads_sites = i1

    end do
    end do

    !---------------------------------------------------------------------------
    !  Allocate and Intialize rate_info structure
    !---------------------------------------------------------------------------
    allocate( hopping_init%rate_info(lat%n_rows*lat%n_cols) )

INTENDED FOR PASS 2 ----------------
    do i=1,lat%n_rows*lat%n_cols
      allocate( hopping_init%rate_info(i)%list( max_avail_ads_sites * &
                                                     max_avail_ads_sites * &
                                                     lat%n_max_nn) )
      hopping_init%rate_info(i)%list = rate_info_hopping( default_int, default_int, default_rate )
    end do
-------------------------------------------------


    ! read rate definitions from the input file
    file_name = c_pars%rate_file_name

    do pass = 1,2

      ! Allocate and initialize rates array before the second pass
      if (pass == 2) then

        allocate( hopping_init%paths(n_hopping_paths) )
        hopping_init%n_paths = n_hopping_paths
        hopping_init%paths = hopping_def( default_int, default_int, default_int, &
                                          default_int, default_int, &
                                          [default_int, default_int], &
                                          default_rate, &
                                          int_law_pars( default_int, [0.0_dp, 0.0_dp] ) )
      endif

      ! reset counter of hopping paths
      n_hopping_paths = 0

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
        !if (ios < 0) buffer = section_end  ! treat end of file as the section end

        ! Split an input string
        tokens = ''
        call split_string(buffer, tokens, ntokens)

        ! skip comments
        if (ntokens == 0) cycle

        select case (parse_state)

          case(parse_state_default)
            ! in parse state default:
            !    word 'hopping' to mark beginning of a hopping section
            !    ignore anything else until hopping section begins

            if (tokens(1) == reaction_names(hopping_id)) then

              hopping_init%is_defined = .true.
              parse_state = parse_state_hopping
              ! reset necessary variables to allow multiple hopping sections
              rct_law_defined  = .false.
              rcic_law_defined = .false.
              rct_law_id_glob  = 0
              rcic_law_id_glob = 0

              if (ntokens == 2) then
                read(tokens(2),'(A)') current_species_name
                current_species_id = get_index(current_species_name, c_pars%ads_names)
                if (current_species_id == 0) call error_message(file_name, line_number, buffer, &
                                                                "unknown species in hopping section definition")
              else
                call error_message(file_name, line_number, buffer, &
                          "hopping key must have 1 parameter -- species")
              end if ! ntokens == 2

            end if ! tokens(1)

          case(parse_state_hopping)
            ! process:
            !    temperature law records,
            !    interaction law records,
            !    from-to records,
            !    hopping vector records:
            !      - 0 1
            !      - 0 1 Arrhenius 1.0E13 0.5
            !      - 0 1 linear 0.5 0.5
            !      - 0 1 Arrhenius 1.0E13 0.5 linear 0.5 0.5
            !    section end

            if (tokens(1) == section_end) then
              parse_state = parse_state_default

            elseif (tokens(1) =='temperature_law') then
              rct_law_id_glob = get_index(tokens(2), rct_law_names)
              if (rct_law_id_glob == 0) then
                call error_message(file_name, line_number, buffer,&
                                  "invalid temperature law statement")
              else
                rct_law_defined = .true.
                select case (rct_law_id_glob)
                  case (Arrhenius_id)
                    if (ntokens/=4) call error_message(file_name, line_number, buffer,&
                                                      "Arrhenius must have 2 parameters")
                    read(tokens(3),*) rct_pars_glob(1)
                    read(tokens(4),*) rct_pars_glob(2)
                  case (extArrhenius_id)
                    if (ntokens/=5) call error_message(file_name, line_number, buffer,&
                                                      "extArrhenius must have 3 parameters")
                    read(tokens(3),*) rct_pars_glob(1)
                    read(tokens(4),*) rct_pars_glob(2)
                    read(tokens(5),*) rct_pars_glob(3)
                end select
              endif

            elseif (tokens(1) =='interaction_law') then
              rcic_law_id_glob = get_index(tokens(2), rcic_law_names)
              if (rcic_law_id_glob == 0) then
                call error_message(file_name, line_number, buffer,&
                                  "invalid interaction law statement")
              else
                rcic_law_defined = .true.
                select case (rcic_law_id_glob)
                  case (rcic_linear_id)
                    if (ntokens/=4) call error_message(file_name, line_number, buffer,&
                                                      "linear interaction law must have 2 parameters")
                    read(tokens(3),*) rcic_pars_glob(1)
                    read(tokens(4),*) rcic_pars_glob(2)
                end select
              endif

            else

              ! check if we have a valid from-to record
              if (ntokens == 4) then
                i1 = get_index(tokens(1),lat_site_names)
                i2 = get_index(tokens(2),ads_site_names)
                i3 = get_index(tokens(3),lat_site_names)
                i4 = get_index(tokens(4),ads_site_names)

                ! check for invalid lst or ast
                if ( i1==0 .or. i2==0 .or. i3==0 .or. i4==0) then
                  print *, 'parse state: ', parse_state
                  call error_message(file_name, line_number, buffer, &
                                    "wrong site name in the hopping section")
                end if

              ! check if we have a valid hopping vector and to-ast record
              elseif (ntokens > 1             .and. &
                      read_int(tokens(1), n1) .and. &
                      read_int(tokens(2), n2)         )  then

                ! we have a valid hopping vector record. Process it

                ! increment counter of hopping paths with account for reversibility
                n_hopping_paths = n_hopping_paths + 2

                if (pass == 2) then

                  ! check for duplicate entry
                  do i=1,n_hopping_paths - 2
                    if ( hopping_init%paths(i)%r     == current_species_id .and. &
                         hopping_init%paths(i)%r_lst == i1                 .and. &
                         hopping_init%paths(i)%r_ast == i2                 .and. &
                         hopping_init%paths(i)%p_lst == i3                 .and. &
                         hopping_init%paths(i)%p_ast == i4                 .and. &
                         hopping_init%paths(i)%p_vec == [n1, n2]           .and. &
                         hopping_init%paths(i)%r_lst == i3                 .and. &
                         hopping_init%paths(i)%r_ast == i4                 .and. &
                         hopping_init%paths(i)%p_lst == i1                 .and. &
                         hopping_init%paths(i)%p_ast == i2                 .and. &
                         hopping_init%paths(i)%p_vec == [-n1, -n2] ) then
                      call error_message(file_name, line_number, buffer, &
                                          "duplicated entry", stop = .false.)
                      duplicate_error = .true.
                    end if
                  end do

                  ! check if energy is defined for all sites involved in the hop
                  if( e_pars%ads_energy(current_species_id, i1, i2) == e_pars%undefined_energy .or. &
                      e_pars%ads_energy(current_species_id, i3, i4) == e_pars%undefined_energy ) then

                      call error_message(file_name, line_number, buffer, &
                                        "rate defined for site with undefined adsorption energy", &
                                        stop=.false., warning=.false.)

                      undefined_energy = .true.
                  end if

                end if

                !   if record has hopping vector information only, 
                !      set the law ids and pars to global default values
                if (ntokens == 2) then
                  if (rct_law_defined .and. rcic_law_defined) then
                    rct_law_id  = rct_law_id_glob
                    rct_pars    = rct_pars_glob
                    rcic_law_id = rcic_law_id_glob
                    rcic_pars   = rcic_pars_glob
                  elseif ( .not. rct_law_defined .and. .not. rcic_law_defined ) then
                      call error_message(file_name, line_number, buffer, &
                                    "temperature and interaction laws are not defined")
                  elseif (.not. rct_law_defined) then
                      call error_message(file_name, line_number, buffer, &
                                    "temperature law is not defined")
                  else
                      call error_message(file_name, line_number, buffer, &
                                    "interaction law is not defined")
                  end if

                else
                  ! set rct_law_id and rcic_law_id to 0  
                  ! to indicate no law on this line found yet
                  rct_law_id  = 0
                  rcic_law_id = 0

                  do i=3,ntokens
                    ! check  for rct law on this line
                    if (get_index(tokens(i), rct_law_names) /= 0) then
                      rct_law_id = get_index(tokens(i), rct_law_names)
                      select case (rct_law_id)
                        case (Arrhenius_id)
                          do j=1,2
                            if ( .not. read_num(tokens(i+j),rct_pars(j)) )&
                              call error_message(file_name, line_number, buffer,&
                                                      "Arrhenius must have 2 numerical parameters")
                          end do
                        case (extArrhenius_id)
                          do j=1,3
                            if ( .not. read_num(tokens(i+j),rct_pars(j)) )&
                              call error_message(file_name, line_number, buffer,&
                                                      "extArrhenius must have 3 numerical parameters")
                          end do
                        case default
                          call error_message(file_name, line_number, buffer, "This should not happen! Check the code!")
                      end select

                    end if

                    ! check if we have an rcic law on this line
                    if (get_index(tokens(i), rcic_law_names) /= 0) then
                      rcic_law_id = get_index(tokens(i), rcic_law_names)
                      select case (rcic_law_id)
                        case (rcic_linear_id)
                          do j=1,2
                            if ( .not. read_num(tokens(i+j),rcic_pars(j)) )&
                              call error_message(file_name, line_number, buffer,&
                                                      "linear interaction must have 2 numerical parameters")
                          end do
                        case default
                          call error_message(file_name, line_number, buffer, "This should not happen! Check the code!")
                      end select
                    end if

                  end do ! i=3,ntokens

                  ! if rct_law  or rcic_law are not defined on this line, 
                  ! set to global defaults
                  if (rct_law_id == 0) then
                    rct_law_id = rct_law_id_glob
                    rct_pars   = rct_pars_glob
                  end if
                  if (rcic_law_id == 0) then
                    rcic_law_id = rcic_law_id_glob
                    rcic_pars   = rcic_pars_glob
                  end if

                  ! check for valid temperature and interaction laws
                  if (rct_law_id == 0) &
                    call error_message(file_name, line_number, buffer, &
                                      "no temperature law is specified")
                  if (rcic_law_id == 0) &
                    call error_message(file_name, line_number, buffer, &
                                      "no interaction law is specified")

                end if ! (ntokens==2)

                ! Check if rct and rcic laws are properly set
                if (rct_law_id == 0) &
                  call error_message(file_name, line_number, buffer, "invalid temperature law")
                if (rcic_law_id == 0) &
                  call error_message(file_name, line_number, buffer, "invalid interaction law")

                if (pass == 2) then

                  ! Set rate constants and rcic for hopping paths
                  hopping_init%paths(n_hopping_paths-1)%r      = current_species_id
                  hopping_init%paths(n_hopping_paths-1)%r_lst  = i1
                  hopping_init%paths(n_hopping_paths-1)%r_ast  = i2
                  hopping_init%paths(n_hopping_paths-1)%p_lst  = i3
                  hopping_init%paths(n_hopping_paths-1)%p_ast  = i4
                  hopping_init%paths(n_hopping_paths-1)%p_vec  = [n1, n2]

                  hopping_init%paths(n_hopping_paths-1)%rcic%id   = rcic_law_id
                  hopping_init%paths(n_hopping_paths-1)%rcic%pars = rcic_pars

                  hopping_init%paths(n_hopping_paths-1)%rate  = 
                                rct_law(rct_law_id, c_pars%temperature, rct_pars)

                  ! Reverse process

                  hopping_init%paths(n_hopping_paths)%r      = current_species_id
                  hopping_init%paths(n_hopping_paths)%r_lst  = i3
                  hopping_init%paths(n_hopping_paths)%r_ast  = i4
                  hopping_init%paths(n_hopping_paths)%p_lst  = i1
                  hopping_init%paths(n_hopping_paths)%p_ast  = i2
                  hopping_init%paths(n_hopping_paths)%p_vec  = [-n1, -n2]

                  hopping_init%paths(n_hopping_paths)%rcic%id   = rcic_law_id
                  ! Reverse rcic parameters array
                  do i=1,size(rcic_pars)
                    hopping_init%paths(n_hopping_paths)%rcic%pars(i) = &
                                                            rcic_pars(size(rcic_pars)-i+1)
                  end do

                  ! State "to" energy minus state "from" energy
                  delta_eps = e_pars%ads_energy(current_species_id, i3, i4) - &
                              e_pars%ads_energy(current_species_id, i1, i2)
                  ! detailed balance                              
                  hopping_init%paths(n_hopping_paths)%rate  = 
                                hopping_init%paths(n_hopping_paths-1)%rate&
                                *exp(c_pars%beta*delta_eps)

                end if ! (pass==2)

              ! we have an invalid record
              else
                call error_message(file_name, line_number, buffer, &
                                  "invalid from-to record in the hopping section")


              end if ! valid from-to record

            endif ! tokens(1) == section_end

          end select ! parse_state
      end do ! while ios=0

      close(inp_unit)

      if (parse_state /= parse_state_default) then
        write(stemp,'(A,I2)') 'parse state: ', parse_state
        call error_message(file_name, line_number, buffer, &
            'hopping: incomplete hopping section. ' // trim(stemp))
      endif

    end do ! pass = 1,2

    if (duplicate_error) &
      call error_message(file_name, 0, '', &
            '  hopping: duplicated hopping entries found')

    if (undefined_energy) then
      call error_message(file_name, 0, '', &
            '  hopping: rates defined for sites with undefined energies')
    else
      if (hopping_init%is_defined) then
        write(*, '(A)') ' hopping: passed check that energies are defined for all rates'
      else
        write(*, '(A)') ' no hopping'
      end if

    end if


    ! ----------------------------------------------------------------------------
  ! Hopping rates report
  ! ----------------------------------------------------------------------------

  if  (hopping_init%is_defined .and. debug(6)) then

    write(*,'(A)') '  Hopping Rates Report.'
    write(*,'(A)') '  --------------------------'
    write(*,*)

    do i=1, hopping_init%n_paths
      write(*,'(5X, A4, A, A11, A11, A, A11, A11, A, F10.3, A, 10F10.3)' &
        c_pars%ads_names(hopping_init%paths(i)%r), ' hops from ',&
        lat_site_names(hopping_init%paths(i)%r_lst), &
        ads_site_names(hopping_init%paths(i)%r_ast), ' to ', &
        lat_site_names(hopping_init%paths(i)%p_lst), &
        ads_site_names(hopping_init%paths(i)%p_ast), ' with rate ', &
        hopping_init%paths(i)%rate , ' and rcic ' &
        trim(rcic_law_names(hopping_init%paths(i)%rcic%id)), &
         hopping_init%paths(i)%rcic%pars )
    end do

  end if

  end function hopping_init

  
!-----------------------------------------------------------------------------
  subroutine construct(this, ads, lat, c_pars, e_pars, beta)
!-----------------------------------------------------------------------------
    class(hopping_type), intent(inout) :: this
    integer, intent(in) :: ads
    class(mc_lat), intent(inout) :: lat
    class(control_parameters), intent(in) :: c_pars
    class(energy_parameters), intent(in) :: e_pars
    real(dp), intent(in) :: beta

    integer :: species, m, iads
    integer :: row_old, col_old, lst_old, ast_old
    integer :: row_new, col_new, lst_new, ast_new
    real(dp) :: energy_old, energy_new, int_energy_old, int_energy_new, int_energy_ts, delta_eps

    ! energy for particle ads in its old position
    energy_old = energy(ads, lat, c_pars, e_pars)
    ! Save the old configuration
    row_old = lat%ads_list(ads)%row
    col_old = lat%ads_list(ads)%col
    lst_old = lat%lst(row_old,col_old)
    ast_old = lat%ads_list(ads)%ast
    species = lat%ads_list(ads)%id

    ! Construct rates for hops to neighbors

    ! Delete particle ads from the old position
    ! we do it here since we never work with occupations inside the following loop
    lat%occupations(row_old,col_old) = 0

    ! Loop over possible new positions of particle ads
    print*, size(this%direction_list(species,lst_old,ast_old,:)), ' should be equal to ', &
            size(this%direction_list(species,lst_old,ast_old), 1), ' and to ', &
            size(this%direction_list, 4)

    do m=1,size(this%direction_list(species,lst_old,ast_old,:))

      ! Get position and site type of neighbour m
      call lat%neighbor2(ads, this%direction_list(species,lst_old,ast_old,m,1:2), row_new, col_new)
      lst_new  = lat%lst(row_new, col_new)

      ! Check if the cell is free
      if (lat%occupations(row_new, col_new) > 0) then

        this%rates(ads,m)%list = 0.0d0

      else

        ! Put particle ads to site m
        lat%ads_list(ads)%row = row_new
        lat%ads_list(ads)%col = col_new
        lat%occupations(row_new,col_new) = ads
        ast_new = this%direction_list(species,lst_old,ast_old,m,3)
        lat%ads_list(ads)%ast = ast_new

        ! Calculate energy of ads in new position
        energy_new = energy(ads, lat, c_pars, e_pars)

        if (debug(10)) then
          write(*,*) ''
          write(*,'(A,I5,A,I5)') "ads_id", id, " ads_no", iads
          write(*,'(A,A4,A4,A,A4,A4)') "from ", lat_site_names(lst_old), ads_site_names(ast_old), &
                                        " to ", lat_site_names(lst_new), ads_site_names(ast_new)
          write(*,'(A,F8.3,A,F8.3)') "E0_i=", e_pars%ads_energy(id, lst_old, ast_old),&
                                      " E0_f=", e_pars%ads_energy(id, lst_new, ast_new)
          write(*,'(A,F8.3,A,F8.3)') "E_i =", energy_old, " E_f =", energy_new
        end if

        ! Calculate interaction correction
        int_energy_old = energy_old - e_pars%ads_energy(id, lst_old, ast_old)
        int_energy_new = energy_new - e_pars%ads_energy(id, lst_new, ast_new)
        int_energy_ts  = rcic_law(this%rate_corr_pars(id, lst_old, ast_old, lst_new, ast_new), &
                                  int_energy_old, int_energy_new)

        ! Barrier correction due to the perturbation
        delta_eps = int_energy_ts - int_energy_old

        ! Add barrier correction if the unperturbed process is uphill
        if ( e_pars%ads_energy(id, lst_new, ast_new) > e_pars%ads_energy(id, lst_old, ast_old) )&
          delta_eps = delta_eps + &
              e_pars%ads_energy(id, lst_new, ast_new) - e_pars%ads_energy(id, lst_old, ast_old)

        this%rates(ads,m)%list(iads) = &
          this%process(id, lst_old, ast_old, lst_new, ast_new)*exp( -beta*delta_eps )

        if (debug(10) .and. ads == 1) then
          write(*,'(A,F8.3)') "delta_eps =", delta_eps
          write(*,'(A,ES10.3)') "uncorrected rate =", this%process(id, lst_old, ast_old, lst_new, ast_new)
          write(*,'(A,F8.3)') "correction energy =",&
          -log(this%rates(ads,m)%list(iads)/this%process(id, lst_old, ast_old, lst_new, ast_new))/beta
          write(*,'(A,F6.3,A,F6.3,A,F6.3,A,F6.3)') "rcic law:", &
                    this%rate_corr_pars(id, lst_old, ast_old, lst_new, ast_new)%pars(1),&
                    ' *', int_energy_old, ' + ', this%rate_corr_pars(id, lst_old, ast_old, lst_new, ast_new)%pars(2),&
                    ' *', int_energy_new
          write(*,'(A,F8.3)') "V_TS=", int_energy_ts
        end if

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
        energy_new = energy(ads, lat, c_pars, e_pars)

!        if (debug(2)) then
!          print*, '-----------construct debug(2)'
!          print*, 'pos new 1 ', lat%ads_list(1)%row, lat%ads_list(1)%col
!          print*, 'pos new 2 ', lat%ads_list(2)%row, lat%ads_list(2)%col
!          print *, 'ads ', ads, 'm ', m,  'E_old = ', energy_old, 'E_new = ', energy_new
!        end if

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
