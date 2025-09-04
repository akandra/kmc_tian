module rates_bimolecular_class

  !---------------------------------------------------------------------------
  !  Module for binary bimolecular reaction  r1 + r2 --> p1 + p2
  !---------------------------------------------------------------------------
  !  * controlled by the bimolecular keyword section of the .reaction file
  !  * p1 is assumed to be at the same lattice site as r1 with its adsorption
  !    site type (ast) choosen from those specified in the .reaction file
  !  * if there exists information that p1 is exclusively formed at one of the
  !    reactant sites, this reactant should be listed as r1; if there are two different
  !    rates, they should be listed separately in the input file
  !  * the symmetric case (like O+O -> O2) is tricky.  THINK MORE ABOUT THIS
  !  *
  !  * THINK AND THEN ACT, NOT OTHER WAY AROUND
  !  *
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
  public    :: bimolecular_init, bimolecular_type, rate_info_bimolecular


  type rate_info_bimolecular
    integer(1) :: proc
    integer(1) :: m
    real(dp)   :: rate
  end type


  type :: v_list_rate_info
    type(rate_info_bimolecular), dimension(:), allocatable :: list
    integer(1) :: n_channels
  end type


  type :: bimolecular_def
    integer  :: r1        ! reactant 1 id
    integer  :: r1_lst    ! reactant 1 lattice site type ( terrace, step, corner)
    integer  :: r1_ast    ! reactant 1 adsorption site type (top, fcc, hcp, etc.)
    integer  :: r2        ! reactant 2 id
    integer  :: r2_lst    ! reactant 2 lattice site type ( terrace, step, corner)
    integer  :: r2_ast    ! reactant 2 adsorption site type (top, fcc, hcp, etc.)
    integer  :: p1        ! product 1 id
    integer  :: p1_lst    ! product 1 lattice site type
    integer  :: p1_ast    ! product 1 adsorbate site type
    integer  :: p2        ! product 2 id
    integer  :: p2_lst    ! product 2 lattice site type
    integer  :: p2_ast    ! product 2 adsorbate site type
    real(dp) :: rate      ! bimolecular rate
    type(int_law_pars) :: rcic
  end type


  type :: bimolecular_type


    logical :: is_defined = .false.


    !---------------------------------------------------------------------------
    ! Reaction channels structure
    !---------------------------------------------------------------------------
    ! Gives the rate of reaction for each set of reactant and product lattice
    ! and adsorption site types given in the bimolecular section of the
    ! .reaction input file
    !---------------------------------------------------------------------------
    type(bimolecular_def), dimension(:), allocatable :: channels

    !---------------------------------------------------------------------------
    ! Reaction rate information structure
    !---------------------------------------------------------------------------
    ! Gives a list of reaction information for a given adsorbate
    !   * list elements give the reaction channel which specifies the rate,
    !     latice site type, and adsorption site type (ast) of reactants and
    !     products
    !---------------------------------------------------------------------------
    !
    !                                 adsorbate      -> which particle
    !                                 .
    type(v_list_rate_info), dimension(:), allocatable :: rate_info

    ! Number of bimolecular channels
    integer :: n_processes

  contains
    procedure :: construct
    procedure :: print

  end type

contains

!------------------------------------------------------------------------------
  function bimolecular_init(c_pars, lat, e_pars)
!------------------------------------------------------------------------------
    type(bimolecular_type) bimolecular_init

    type(control_parameters), intent(in)    :: c_pars
    type(mc_lat)            , intent(in)    :: lat
    type(energy_parameters) , intent(in)    :: e_pars

    integer :: i, j, ios, nwords, line_number, m

    logical :: undefined_energy

    character(len=max_string_length)                :: buffer
    character(len=max_string_length)                :: words(100)
    character(len=len(trim(c_pars%rate_file_name))) :: file_name

    character(len=10)     :: r1_name, r2_name, p1_name, p2_name
    integer               :: r1_id,   r2_id,   p1_id,   p2_id
    integer               :: rct_law_id, rct_law_id_glob
    integer               :: rcic_law_id, rcic_law_id_glob

    integer               :: parse_state
    integer, parameter    :: parse_state_default =  0
    integer, parameter    :: parse_state_bimolecular =  bimolecular_id
    integer, parameter    :: n_bimolecular_species = 4

    logical :: rct_law_defined  = .false.
    logical :: rcic_law_defined = .false.

    real(dp), dimension(n_max_rct_pars )::  rct_pars,  rct_pars_glob
    real(dp), dimension(n_max_rcic_pars):: rcic_pars, rcic_pars_glob

    integer,  parameter   :: default_int = 0
    real(dp), parameter   :: default_rate  = -1.0_dp

    integer, dimension(n_bimolecular_species) :: lst, ast

    integer :: max_avail_ads_sites
    integer :: n_bimolecular_channels = 0
    integer :: bimolecular_counter = 0

    logical :: duplicate_error = .false.
    logical ::        rp_error = .false.

    ! maximal number of available ads. sites
    max_avail_ads_sites = 1
    do i=1,c_pars%n_species
    do m=1,size(lat%avail_ads_sites(i,:))
      j = size(lat%avail_ads_sites(i,m)%list)
      if (max_avail_ads_sites < j) max_avail_ads_sites = j
    end do
    end do

    !---------------------------------------------------------------------------
    !  Allocate and Intialize rate_info structure
    !---------------------------------------------------------------------------
    allocate( bimolecular_init%rate_info(lat%n_rows*lat%n_cols) )

    do i=1,lat%n_rows*lat%n_cols
      allocate( bimolecular_init%rate_info(i)%list( lat%n_max_nn * &
                                                     max_avail_ads_sites  ) )
      bimolecular_init%rate_info(i)%list = rate_info_bimolecular( default_int, default_int, default_rate )
    end do

    !  read rate definitions from the input file
    file_name = c_pars%rate_file_name
    call open_for_read(inp_unit, file_name )


!-------------------------------------------------------------------------------
!   First pass of processing the .reaction file:
!     * define the number of bimolecular reactions
!     * check for parsing errors
!-------------------------------------------------------------------------------
    ios = 0
    parse_state = parse_state_default
    line_number = 0

    do while (ios == 0)

      read(inp_unit, '(A)', iostat=ios) buffer
      line_number = line_number + 1
        ! ios < 0: end of record or endfile condition detected
        ! ios > 0: an error is detected
        ! ios = 0  otherwise

        ! Split an input string
        words = ''
        call split_string(buffer, words, nwords)

      ! skip comments
      if (nwords == 0) cycle

      select case (parse_state)

        case(parse_state_default)
          ! in parse state default:
          !    word 'bimolecular' to mark beginning of a bimolecular section
          !    ignore anything else until bimolecular section begins

          if (words(1) == reaction_names(bimolecular_id)) then

            parse_state = parse_state_bimolecular

          end if  ! words(1)

        case(parse_state_bimolecular)

          if (words(1) == section_end) then
            parse_state = parse_state_default
          ! don't count lines defining rct and rcic
          elseif (words(1) =='temperature_law' .or. &
                  words(1) =='interaction_law') then
            cycle

          else
            n_bimolecular_channels = n_bimolecular_channels + 1

          endif ! words(1) == section_end

      end select                                      ! select case(words(1))

    end do ! while ios=0

    close(inp_unit)

    allocate(bimolecular_init%channels(n_bimolecular_channels))
    bimolecular_init%n_processes = n_bimolecular_channels
    bimolecular_init%channels = bimolecular_def(0,0,0,0,0,0,0,0,0,0,0,0,0.0_dp,int_law_pars(0,0.0_dp))

!-------------------------------------------------------------------------------
!   Second pass of processing the .reaction file:
!     * build the channel structure
!-------------------------------------------------------------------------------
    call open_for_read(inp_unit, file_name )

    ios = 0
    parse_state = parse_state_default
    line_number = 0
    undefined_energy = .false.

    do while (ios == 0)

      read(inp_unit, '(A)', iostat=ios) buffer
      line_number = line_number + 1
        ! ios < 0: end of record or endfile
        ! ios > 0: an error is detected
        ! ios = 0  otherwise

        ! Split an input string
        words = ''
        call split_string(buffer, words, nwords)

      ! skip comments
      if (nwords == 0) cycle

      select case (parse_state)

        case(parse_state_default)
          ! in parse state default:
          !    word 'bimolecular' to mark beginning of a bimolecular section
          !    ignore anything else until bimolecular section begins

          if (words(1) == reaction_names(bimolecular_id)) then

            parse_state = parse_state_bimolecular
            bimolecular_init%is_defined = .true.
            ! reset necessary to allow multiple bimolecular sections
            rct_law_defined  = .false.
            rcic_law_defined = .false.
            rct_law_id_glob  = 0
            rcic_law_id_glob = 0

            if (nwords == n_bimolecular_species + 1) then

              read(words(2),'(A)') r1_name
              r1_id = get_index(r1_name, c_pars%ads_names )

              read(words(3),'(A)') r2_name
              r2_id = get_index(r2_name, c_pars%ads_names )

              read(words(4),'(A)') p1_name
              p1_id = get_index(p1_name, c_pars%ads_names )

              read(words(5),'(A)') p2_name
              p2_id = get_index(p2_name, c_pars%ads_names )


              if (r1_id == 0 .or. r2_id == 0 .or. p1_id == 0 .or. p2_id == 0) &
                  call error_message(file_name, line_number, buffer, &
                            "unknown species in bimolecular section definition")
            else
              call error_message(file_name, line_number, buffer, &
                         "bimolecular key must have 4 parameter -- reactant_1 reactant_2 product_1 product_2")
            end if !nwords == n_bimolecular_species + 1

          end if  ! words(1)

        case(parse_state_bimolecular)
          ! process:
          !    temperature law records,
          !    interaction law records,
          !    'from-from-to-to' records,
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
            ! check if we have a valid from-from-to-to rate record
            ! ---------------------------------------------------------
            if (nwords < 2*n_bimolecular_species) &
              call error_message(file_name, line_number, buffer, &
                                 "invalid 'from' rate record in the bimolecular section")

            bimolecular_counter = bimolecular_counter + 1

            do i=1,n_bimolecular_species
              lst(i) = get_index(words(2*i-1),lat_site_names)
              ast(i) = get_index(words(2*i  ),ads_site_names)
            end do

            if ( any(lst(1:n_bimolecular_species)==0) ) &
              call error_message(file_name, line_number, buffer, &
                         "wrong lattice site name in the bimolecular section")
            if ( any(ast(1:n_bimolecular_species)==0) ) &
              call error_message(file_name, line_number, buffer, &
                         "wrong adsorption site name in the bimolecular section")

            do i=1,bimolecular_counter-1
              if ( bimolecular_init%channels(i)%r1     == r1_id  .and. &
                   bimolecular_init%channels(i)%r1_lst == lst(1) .and. &
                   bimolecular_init%channels(i)%r1_ast == ast(1) .and. &
                   bimolecular_init%channels(i)%r2     == r2_id  .and. &
                   bimolecular_init%channels(i)%r2_lst == lst(2) .and. &
                   bimolecular_init%channels(i)%r2_ast == ast(2) .and. &
                   bimolecular_init%channels(i)%p1     == p1_id  .and. &
                   bimolecular_init%channels(i)%p1_lst == lst(3) .and. &
                   bimolecular_init%channels(i)%p1_ast == ast(3) .and. &
                   bimolecular_init%channels(i)%p2     == p2_id  .and. &
                   bimolecular_init%channels(i)%p2_lst == lst(4) .and. &
                   bimolecular_init%channels(i)%p2_ast == ast(4)       ) then
                call error_message(file_name, line_number, buffer, &
                                   "duplicated entry", stop = .false.)
                  duplicate_error = .true.
              end if
            end do

            if ( lst(1) /= lst(3) ) then
              call error_message(file_name, line_number, buffer, &
                  "r1 and p1 lattice site types have to be the same", stop = .false.)
              rp_error = .true.
            end if
            if ( lst(2) /= lst(4) ) then
              call error_message(file_name, line_number, buffer, &
                  "r2 and p2 lattice site types have to be the same", stop = .false.)
              rp_error = .true.
            end if

            ! check energy is defined for reactant's and for products' lat_site_type and ads_site
            if( e_pars%ads_energy(r1_id, lst(1), ast(1)) == e_pars%undefined_energy .or. &
                e_pars%ads_energy(r2_id, lst(2), ast(2)) == e_pars%undefined_energy .or. &
                e_pars%ads_energy(p1_id, lst(3), ast(3)) == e_pars%undefined_energy .or. &
                e_pars%ads_energy(p2_id, lst(4), ast(4)) == e_pars%undefined_energy ) then

              call error_message(file_name, line_number, buffer, &
                                   "rate defined for site with undefined adsorption energy", &
                                   stop=.false., warning=.false.)
                undefined_energy = .true.
            end if

            ! ---------------------------------------------------------
            ! we have a valid rate record. Process it
            ! ---------------------------------------------------------
            ! if record only has 'from-from-to-to' site information then set the law ids and pars to global default values
            if (nwords == 2*n_bimolecular_species) then
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
              ! set rct_law_id and rcic_law_id to 0 to indicate no law on this line found yet
              rct_law_id  = 0
              rcic_law_id = 0

              do i=2*n_bimolecular_species+1,nwords
                ! check  for rct law on this line
                if (get_index(words(i), rct_law_names) /= 0) then
                  rct_law_id = get_index(words(i), rct_law_names)
                  select case (rct_law_id)
                    case (Arrhenius_id)
                      do j=1,2
                        if ( .not. read_num(words(i+j),rct_pars(j)) )&
                          call error_message(file_name, line_number, buffer,&
                                                  "Arrhenius must have 2 numerical parameters")
                      end do
                    case (extArrhenius_id)
                      do j=1,3
                        if ( .not. read_num(words(i+j),rct_pars(j)) )&
                          call error_message(file_name, line_number, buffer,&
                                                  "extArrhenius must have 3 numerical parameters")
                      end do
                    case default
                      call error_message(file_name, line_number, buffer, "This should not happen! Check the code!")
                  end select

                end if

                ! check if we have an rcic law on this line
                if (get_index(words(i), rcic_law_names) /= 0) then
                  rcic_law_id = get_index(words(i), rcic_law_names)
                  select case (rcic_law_id)
                    case (rcic_linear_id)
                      do j=1,2
                        if ( .not. read_num(words(i+j),rcic_pars(j)) )&
                          call error_message(file_name, line_number, buffer,&
                                                  "linear interaction must have 2 numerical parameters")
                      end do
                    case default
                      call error_message(file_name, line_number, buffer, "This should not happen! Check the code!")
                  end select
                end if

              end do ! i=2*n_bimolecular_species+1,nwords

              ! if rct_law  or rcic_law are not defined on this line, set to global defaults
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

            end if ! (nwords==2*n_bimolecular_species)

            ! Check if rct and rcic laws are properly set
            if (rct_law_id == 0) &
              call error_message(file_name, line_number, buffer, "invalid temperature law")
            if (rcic_law_id == 0) &
              call error_message(file_name, line_number, buffer, "invalid interaction law")

            ! Set rate constants and rcic for bimolecular channels
            bimolecular_init%channels(bimolecular_counter)%r1      = r1_id
            bimolecular_init%channels(bimolecular_counter)%r1_lst  = lst(1)
            bimolecular_init%channels(bimolecular_counter)%r1_ast  = ast(1)
            bimolecular_init%channels(bimolecular_counter)%r2      = r2_id
            bimolecular_init%channels(bimolecular_counter)%r2_lst  = lst(2)
            bimolecular_init%channels(bimolecular_counter)%r2_ast  = ast(2)
            bimolecular_init%channels(bimolecular_counter)%p1      = p1_id
            bimolecular_init%channels(bimolecular_counter)%p1_lst  = lst(3)
            bimolecular_init%channels(bimolecular_counter)%p1_ast  = ast(3)
            bimolecular_init%channels(bimolecular_counter)%p2      = p2_id
            bimolecular_init%channels(bimolecular_counter)%p2_lst  = lst(4)
            bimolecular_init%channels(bimolecular_counter)%p2_ast  = ast(4)


            bimolecular_init%channels(bimolecular_counter)%rate = &
                          rct_law(rct_law_id, c_pars%temperature, rct_pars)

            bimolecular_init%channels(bimolecular_counter)%rcic%id   = rcic_law_id
            bimolecular_init%channels(bimolecular_counter)%rcic%pars = rcic_pars

          endif ! words(1) == section_end

      end select

    end do ! while ios=0

    close(inp_unit)

    if (duplicate_error) stop 990
    if (       rp_error) stop 989

    if (undefined_energy) then
!      write(*, '(A)') ' bimolecular: error, rates defined for sites with undefined energies'
      stop 988

    else

      if (bimolecular_init%is_defined) then
        write(6, '(A/)') ' bimolecular: passed check that energies are defined for all rates'
      else
        write(6, '(A)') ' no bimolecular reactions'
      end if

    end if

  end function bimolecular_init

!-------------------------------------------------------------------------------
!  subroutine construct
!-------------------------------------------------------------------------------
  subroutine construct(this, ads, lat, c_pars, e_pars, beta)

    class(bimolecular_type), intent(inout) :: this
    integer, intent(in) :: ads
    class(mc_lat), intent(inout) :: lat
    class(control_parameters), intent(in) :: c_pars
    class(energy_parameters), intent(in) :: e_pars
    real(dp), intent(in) :: beta

    integer :: m, iprocs
    integer :: id_r1, row_r1, col_r1, lst_r1, ast_r1
    integer :: id_r2, row_r2, col_r2, lst_r2, ast_r2, ads_r2
    integer :: channel
    real(dp):: int_energy, int_energy_ts, delta_eps

    ! Consider ads as reactant 1
    row_r1 = lat%ads_list(ads)%row
    col_r1 = lat%ads_list(ads)%col
    lst_r1 = lat%lst(row_r1,col_r1)
    ast_r1 = lat%ads_list(ads)%ast
    id_r1  = lat%ads_list(ads)%id

    ! Loop over possible positions for reactant 2
    channel = 0
    do m=1,lat%n_nn(lst_r1,1)

      ! Get position and site type of neighbour m
      call lat%neighbor(ads, m, row_r2, col_r2)

      ! Check if the cell is occupied
      ads_r2 = lat%occupations(row_r2, col_r2)

      if (ads_r2/=0) then

        ! Consider neighbour m as a possible reactant 2
        lst_r2 = lat%lst(row_r2,col_r2)
        ast_r2 = lat%ads_list(ads_r2)%ast
        id_r2  = lat%ads_list(ads_r2)%id

        do iprocs=1, this%n_processes

          if ( this%channels(iprocs)%r1     == id_r1  .and. &
               this%channels(iprocs)%r2     == id_r2  .and. &
               this%channels(iprocs)%r1_lst == lst_r1 .and. &
               this%channels(iprocs)%r2_lst == lst_r2 .and. &
               this%channels(iprocs)%r1_ast == ast_r1 .and. &
               this%channels(iprocs)%r2_ast == ast_r2     ) &
          then

            channel = channel + 1

            ! Calculate interaction correction
            ! WARNING! here we assume that the correction does not depend on products
            !          consider to introduce it later
            int_energy    = energy(ads,lat,c_pars,e_pars)    - e_pars%ads_energy(id_r1, lst_r1, ast_r1)&
                          + energy(ads_r2,lat,c_pars,e_pars) - e_pars%ads_energy(id_r2, lst_r2, ast_r2)
            int_energy_ts = rcic_law(this%channels(iprocs)%rcic, int_energy, 0.0_dp)
            ! Barrier correction due to the perturbation
            delta_eps = int_energy_ts - int_energy

!            ! Get the information on products
!            id_p1  = this%channels(iprocs)%p1
!            lst_p1 = this%channels(iprocs)%p1_lst
!            ast_p1 = this%channels(iprocs)%p1_ast
!            id_p2  = this%channels(iprocs)%p2
!            lst_p2 = this%channels(iprocs)%p2_lst
!            ast_p2 = this%channels(iprocs)%p2_ast

            this%rate_info(ads)%list(channel)%proc   = iprocs
            this%rate_info(ads)%list(channel)%m      = m
            ! WARNING: Decide later if we need the rate field
            this%rate_info(ads)%list(channel)%rate  = this%channels(iprocs)%rate&
                                                    *exp( -beta*delta_eps)

          end if

        end do ! iprocs

      end if ! occupations

    end do ! m

    ! save the number of channels in the rate_info structure
    this%rate_info(ads)%n_channels = channel

end subroutine construct

!------------------------------------------------------------------------------
  subroutine print(this, c_pars)
!------------------------------------------------------------------------------
    class(bimolecular_type), intent(in) :: this

    class(control_parameters), intent(in) :: c_pars

    integer :: i
    character(:), allocatable :: r1_name, r2_name, p1_name, p2_name

    print*
    print*, 'Bimolecular Rates:'
    print '(A)', ' ---------------------------'
    do i=1,size(this%channels,1)
      r1_name = c_pars%ads_names(this%channels(i)%r1 )
      r2_name = c_pars%ads_names(this%channels(i)%r2)
      p1_name = c_pars%ads_names(this%channels(i)%p1)
      p2_name = c_pars%ads_names(this%channels(i)%p2)
      write(*,'(1x,A,A,A,2X,A,A,2X,A,A,2X,A,A,1pe11.2)')  &
              trim(r1_name)// ' + '// trim(r2_name)// ' -> '// trim(p1_name)//' + '// trim(p2_name)//': ',&
              lat_site_names(this%channels(i)%r1_lst), ads_site_names(this%channels(i)%r1_ast), &
              lat_site_names(this%channels(i)%r2_lst), ads_site_names(this%channels(i)%r2_ast), &
              lat_site_names(this%channels(i)%p1_lst), ads_site_names(this%channels(i)%p1_ast), &
              lat_site_names(this%channels(i)%p2_lst), ads_site_names(this%channels(i)%p2_ast), &
              this%channels(i)%rate
    end do
    print*

  end subroutine print

end module rates_bimolecular_class
