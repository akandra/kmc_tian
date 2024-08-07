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
    integer               :: law_id

    integer               :: parse_state
    integer, parameter    :: parse_state_ignore  = -1
    integer, parameter    :: parse_state_default =  0
    integer, parameter    :: parse_state_bimolecular =  bimolecular_id
    integer, parameter    :: n_bimolecular_species = 4


    real(dp), dimension(3):: pars = 0.0_dp

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

      if (ios /= 0) exit

        ! Split an input string
        words = ''
        call split_string(buffer, words, nwords)

        select case (words(1)) ! take a keyword
!------------------------------------------------------------------------------
          case('bimolecular')                       ! of select case (words(1)
!------------------------------------------------------------------------------
            bimolecular_init%is_defined = .true.

            if (parse_state /= parse_state_default) &
              call error_message(file_name, line_number, buffer, &
                         "invalid ending of the reaction section")
            parse_state = parse_state_bimolecular
            if ( nwords/=6) call error_message(file_name, line_number, buffer, &
                               "bimolecular key must have 5 parameters")

            read(words(2),'(A)') r1_name
            r1_id = get_index(r1_name, c_pars%ads_names )
            if (r1_id == 0) call error_message(file_name, line_number, buffer, &
                                         "inconsistent bimolecular reactant 1 definition")

            read(words(3),'(A)') r2_name
            r2_id = get_index(r2_name, c_pars%ads_names )
            if (r2_id == 0) call error_message(file_name, line_number, buffer, &
                                         "inconsistent bimolecular reactant 2 definition")

            read(words(4),'(A)') p1_name
            p1_id = get_index(p1_name, c_pars%ads_names )
            if (p1_id == 0) call error_message(file_name, line_number, buffer, &
                                         "inconsistent bimolecular product 1 definition")

            read(words(5),'(A)') p2_name
            p2_id = get_index(p2_name, c_pars%ads_names )
            if (p2_id == 0) call error_message(file_name, line_number, buffer, &
                                         "inconsistent bimolecular product 2 definition")

            law_id = get_index(words(nwords), rct_law_names )
            if (law_id == 0) call error_message(file_name, line_number, buffer, &
                                                  "bimolecular: unknown temperature law")
!-------------------------------------------------------------------------------
          case ('terrace','step','corner')            ! of select case(words(1))
!-------------------------------------------------------------------------------

            select case (parse_state)

              case(parse_state_ignore)
                ! ignore
                ! print *, 'warning ignoring line', line_number, buffer

              case(bimolecular_id)

                n_bimolecular_channels = n_bimolecular_channels + 1

              case default
                call error_message(file_name, line_number, buffer, "bimolecular: invalid site type statement")

            end select

!-------------------------------------------------------------------------------
          case('')                                    ! of select case(words(1))
!-------------------------------------------------------------------------------
            if (buffer == '') then
              parse_state = parse_state_default
!              print*, 'blank line '
!            else
!              print*, 'comment: ', trim(buffer)
            end if

!-------------------------------------------------------------------------------
          case default                                ! of select case(words(1))
!-------------------------------------------------------------------------------
            if ( parse_state == parse_state_default .and. get_index(words(1),reaction_names) /= 0 ) &
              parse_state = parse_state_ignore

            if (parse_state /= parse_state_ignore) &
              call error_message(file_name, line_number, buffer, "bimolecular: unknown key")

        end select                                      ! select case(words(1))

    end do ! while ios=0

    close(inp_unit)

    allocate(bimolecular_init%channels(n_bimolecular_channels))
    bimolecular_init%n_processes = n_bimolecular_channels
    bimolecular_init%channels = bimolecular_def(0,0,0,0,0,0,0,0,0,0,0,0,0.0_dp)

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

      if (ios /= 0) exit

        ! Split an input string
        words = ''
        call split_string(buffer, words, nwords)

        select case (words(1)) ! take a keyword
!-------------------------------------------------------------------------------
          case('bimolecular')                        ! of select case (words(1)
!-------------------------------------------------------------------------------
            bimolecular_init%is_defined = .true.
            parse_state = parse_state_bimolecular

            read(words(2),'(A)') r1_name
            r1_id = get_index(r1_name, c_pars%ads_names )

            read(words(3),'(A)') r2_name
            r2_id = get_index(r2_name, c_pars%ads_names )

            read(words(4),'(A)') p1_name
            p1_id = get_index(p1_name, c_pars%ads_names )

            read(words(5),'(A)') p2_name
            p2_id = get_index(p2_name, c_pars%ads_names )

            law_id = get_index(words(nwords), rct_law_names )

!            ! debugging printout
!            print*, 'reactant1 name and id: ', r1_name, r1_id
!            print*, 'reactant2 name and id: ', r2_name, r2_id
!            print*, 'product1  name and id: ', p1_name, p1_id
!            print*, 'product2  name and id: ', p2_name, p2_id
!            print*, c_pars%ads_names
!            stop 111

!-------------------------------------------------------------------------------
          case ('terrace','step','corner')            ! of select case(words(1))
!-------------------------------------------------------------------------------

            select case (parse_state)

              case(parse_state_ignore)
                ! ignore
                ! print *, 'warning ignoring line', line_number, buffer

              case(bimolecular_id)

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
                    call error_message(file_name, line_number, buffer, "duplicated entry", stop = .false.)
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

                select case (law_id)

                  case (Arrhenius_id)
                    if (nwords/=10) call error_message(file_name, line_number, buffer,&
                                              "Arrhenius must have 2 parameters")
                    read(words( 9),*) pars(1)
                    read(words(10),*) pars(2)
                    bimolecular_init%channels(bimolecular_counter)%rate = arrhenius(c_pars%temperature, pars(1:2))

                  case (extArrhenius_id)
                    if (nwords/=11) call error_message(file_name, line_number, buffer,&
                                              "extArrhenius must have 3 parameters")
                    read(words( 9),*) pars(1)
                    read(words(10),*) pars(2)
                    read(words(11),*) pars(3)
                    bimolecular_init%channels(bimolecular_counter)%rate = extArrhenius(c_pars%temperature, pars(1:3))

                  case default
                    call error_message(file_name, line_number, buffer, "bimolecular: this should not happen! Check the code!")

                end select

                 ! Debugging printout
!                 print*, 'reaction: ', reaction_names(parse_state),&
!                        ' for species:', r1_name, r2_name, p1_name, p2_name
!                 print*, 'law: ',  rct_law_names(law_id),&
!                        ' from: ', lat_site_names(lst(1)),ads_site_names(ast(1)),&
!                        ' and: ',  lat_site_names(lst(2)),ads_site_names(ast(2)),&
!                        ' to: ',   lat_site_names(lst(3)),ads_site_names(ast(3)),&
!                        ' and: ',  lat_site_names(lst(4)),ads_site_names(ast(4))
!                print'(A,3f16.3)', 'with pars: ', pars
!                pause

              case default
                call error_message(file_name, line_number, buffer, "bimolecular: invalid site type statement")

            end select

!-------------------------------------------------------------------------------
          case('')                                    ! of select case(words(1))
!-------------------------------------------------------------------------------
            if (buffer == '') then
              parse_state = parse_state_default
!              print*, 'blank line '
!            else
!              print*, 'comment: ', trim(buffer)
            end if

!-------------------------------------------------------------------------------
          case default                                ! of select case(words(1))
!-------------------------------------------------------------------------------
            if ( parse_state == parse_state_default .and. get_index(words(1),reaction_names) /= 0 ) &
              parse_state = parse_state_ignore

            if (parse_state /= parse_state_ignore) &
              call error_message(file_name, line_number, buffer, "bimolecular: unknown key")

        end select                                       ! select case(words(1))

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
!dja  subroutine construct(this, ads, lat, e_pars, beta)
  subroutine construct(this, ads, lat)

    class(bimolecular_type), intent(inout) :: this
    integer, intent(in) :: ads
    class(mc_lat), intent(inout) :: lat
!dja    class(energy_parameters), intent(in) :: e_pars
!dja    real(dp), intent(in) :: beta

    integer :: m, iprocs
    integer :: id_r1, row_r1, col_r1, lst_r1, ast_r1
    integer :: id_r2, row_r2, col_r2, lst_r2, ast_r2, ads_r2
    integer :: channel

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
            this%rate_info(ads)%list(channel)%rate  = this%channels(iprocs)%rate

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
