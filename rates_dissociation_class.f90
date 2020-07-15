module rates_dissociation_class

  !---------------------------------------------------------------------------
  !  Module for binary dissocation reaction  r --> p1 + p2
  !---------------------------------------------------------------------------
  !  * controlled by the dissociation keyword section of the .reaction file
  !  * p1 is assumed to be at the same lattice site as r with its adsorption
  !    site type (ast) choosen from those specified in the .reaction file
  !  * p2 lattice site is taken as a nearest neighbor of r with its ast
  !    taken from those specified in the .reaction file
  !---------------------------------------------------------------------------

  use constants
  use control_parameters_class
  use mc_lat_class
  use energy_parameters_class
  use energy_mod
  use open_file
  use utilities

  implicit none

  private
  public    :: dissociation_init, dissociation_type, rate_info


  type rate_info
    integer(1) :: proc
    integer(1) :: m
    real(dp)   :: rate
  end type


  type :: v_list_rate_info
    type(rate_info), dimension(:), allocatable :: list
    integer(1) :: n_channels
  end type


  type :: dissociation_def
    integer  :: r         ! reactant id
    integer  :: r_lst     ! reactant lattice site type ( terrace, step, corner)
    integer  :: r_ast     ! reactant adsorption site type (top, fcc, hcp, etc.)
    integer  :: p1        ! product 1 id
    integer  :: p1_lst    ! product 1 lattice site type
    integer  :: p1_ast    ! product 1 adsorbate site type
    integer  :: p2        ! product 2 id
    integer  :: p2_lst    ! product 2 lattice site type
    integer  :: p2_ast    ! product 2 adsorbate site type
    real(dp) :: rate      ! dissociation rate
  end type


  type :: dissociation_type


    logical :: is_defined = .false.


    !---------------------------------------------------------------------------
    ! Reaction channels structure
    !---------------------------------------------------------------------------
    ! Gives the rate of reaction for each set of reactant and product lattice
    ! and adsorption site types given in the dissociation section of the
    ! .reaction input file
    !---------------------------------------------------------------------------
    type(dissociation_def), dimension(:), allocatable :: channels

    !---------------------------------------------------------------------------
    ! Reaction rate information structure
    !---------------------------------------------------------------------------
    ! Gives a list of reaction information for a given adsorbate
    !   * list elements give the reaction channel which specifies the rate,
    !     latice site type, and adsorption site type (ast) of reactants and
    !     products
    !   * the structure also gives the direction of the site for product 2
    !     relative to that of the reactant
    !---------------------------------------------------------------------------
    !
    !                                 adsorbate      -> which particle
    !                                 .
    type(v_list_rate_info), dimension(:), allocatable :: rate_info



    ! List of additional nn directions to scan after hop
    integer, dimension(:,:), allocatable :: nn_new

    ! Number of dissociation channels
    integer :: n_processes

  contains
    procedure :: construct
    procedure :: print

  end type

contains

!------------------------------------------------------------------------------
  function dissociation_init(c_pars, lat, e_pars)
!------------------------------------------------------------------------------
    type(dissociation_type) dissociation_init

    type(control_parameters), intent(in)    :: c_pars
    type(mc_lat)            , intent(in)    :: lat
    type(energy_parameters) , intent(in)    :: e_pars

    integer :: i, ios, nwords, line_number, i1, i2, i3, i4, i5, i6, m

    integer :: species, st1, st2, ast1, ast2
    logical :: e_defined1, e_defined2, r_defined, undefined_rate, undefined_energy

    character(len=max_string_length)                :: buffer
    character(len=max_string_length)                :: words(100)
    character(len=len(trim(c_pars%rate_file_name))) :: file_name
    character(1)                                    :: answer

    character(len=10)     :: current_reactant_name, current_product1_name, current_product2_name
    integer               :: current_reactant_id, current_product1_id, current_product2_id
    character(len=20)     :: current_law_name
    integer               :: current_law_id

    integer               :: parse_state
    integer, parameter    :: parse_state_ignore  = -1
    integer, parameter    :: parse_state_default =  0
    integer, parameter    :: parse_state_dissociation =  dissociation_id


    real(dp), dimension(3):: pars = 0.0_dp

    integer,  parameter   :: default_int = 0
    real(dp), parameter   :: default_rate  = -1.0_dp

    integer :: row, col, site, id
    integer :: n_nn, n_nn2, max_avail_ads_sites
    integer :: n_dissociation_channels = 0
    integer :: dissociation_counter = 0

    logical :: duplicate_error = .false.

    n_nn  = lat%n_nn(1)
    n_nn2 = n_nn/2
    ! List of additional nn directions to scan after dissociation
    allocate(dissociation_init%nn_new(n_nn,n_nn2))
    do m=1,n_nn
    do i=1,n_nn2
      dissociation_init%nn_new(m,i) = modulo( m+i-n_nn2, n_nn ) + 1
    end do
    end do

    ! maximal number of available ads. sites
    max_avail_ads_sites = 1
    do i=1,c_pars%n_species
    do m=1,n_max_site_types
      i1 = size(lat%avail_ads_sites(i,m)%list)
      if (max_avail_ads_sites < i1) max_avail_ads_sites = i1
    end do
    end do


    !---------------------------------------------------------------------------
    !  Allocate and Intialize rate_info structure
    !---------------------------------------------------------------------------
    allocate( dissociation_init%rate_info(lat%n_rows*lat%n_cols) )

    do i=1,lat%n_rows*lat%n_cols
      allocate( dissociation_init%rate_info(i)%list( max_avail_ads_sites * &
                                                     max_avail_ads_sites * &
                                                     lat%n_nn(1)) )
      dissociation_init%rate_info(i)%list = rate_info( 0, 0, 0.0_dp )
    end do

    !  read rate definitions from the input file
    file_name = c_pars%rate_file_name
    call open_for_read(inp_unit, file_name )


!-------------------------------------------------------------------------------
!   First pass of processing the .reaction file:
!     * define the number of dissociation reactions
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
          case('dissociation')                       ! of select case (words(1)
!------------------------------------------------------------------------------
            dissociation_init%is_defined = .true.

            if (parse_state /= parse_state_default) &
              call error_message(file_name, line_number, buffer, &
                         "invalid ending of the reaction section")
            parse_state = parse_state_dissociation
            if (nwords/=5) call error_message(file_name, line_number, buffer, &
                               "dissociation key must have 4 parameters")

            read(words(2),'(A)') current_reactant_name
            current_reactant_id = get_index(current_reactant_name, c_pars%ads_names )
            if (current_reactant_id == 0) call error_message(file_name, line_number, buffer, &
                                                  "inconsistent dissociation reactant definition")

            read(words(3),'(A)') current_product1_name
            current_product1_id = get_index(current_product1_name, c_pars%ads_names )
            if (current_product1_id == 0) call error_message(file_name, line_number, buffer, &
                                                  "inconsistent dissociation product 1 definition")

            read(words(4),'(A)') current_product2_name
            current_product2_id = get_index(current_product2_name, c_pars%ads_names )
            if (current_product2_id == 0) call error_message(file_name, line_number, buffer, &
                                                  "inconsistent dissociation product 2 definition")

            current_law_id = get_index(words(5), law_names )
            if (current_law_id == 0) call error_message(file_name, line_number, buffer, &
                                                  "dissociation: unknown temperature law")
!            print*, 'reactant name and id: ', current_reactant_name, current_reactant_id
!            print*, 'product1 name and id: ', current_product1_name, current_product1_id
!            print*, 'product2 name and id: ', current_product2_name, current_product2_id
!            print*, c_pars%ads_names
!            stop 111

!-------------------------------------------------------------------------------
          case ('terrace','step','corner')            ! of select case(words(1))
!-------------------------------------------------------------------------------


            select case (parse_state)

              case(parse_state_ignore)
                ! ignore
                ! print *, 'warning ignoring line', line_number, buffer

              case(dissociation_id)

                n_dissociation_channels = n_dissociation_channels + 1

              case default
                call error_message(file_name, line_number, buffer, "Dissociation: invalid site type statement")

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
              call error_message(file_name, line_number, buffer, "Dissociation: unknown key")

        end select                                      ! select case(words(1))

    end do ! while ios=0

    close(inp_unit)

    allocate(dissociation_init%channels(n_dissociation_channels))
    dissociation_init%n_processes = n_dissociation_channels


!-------------------------------------------------------------------------------
!   Second pass of processing the .reaction file:
!     * build the process structure
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
          case('dissociation')                        ! of select case (words(1)
!-------------------------------------------------------------------------------
            dissociation_init%is_defined = .true.
            parse_state = parse_state_dissociation

            read(words(2),'(A)') current_reactant_name
            current_reactant_id = get_index(current_reactant_name, c_pars%ads_names )

            read(words(3),'(A)') current_product1_name
            current_product1_id = get_index(current_product1_name, c_pars%ads_names )

            read(words(4),'(A)') current_product2_name
            current_product2_id = get_index(current_product2_name, c_pars%ads_names )

            current_law_id = get_index(words(5), law_names )

!            ! debugging printout
!            print*, 'reactant name and id: ', current_reactant_name, current_reactant_id
!            print*, 'product1 name and id: ', current_product1_name, current_product1_id
!            print*, 'product2 name and id: ', current_product2_name, current_product2_id
!            print*, c_pars%ads_names
!            stop 111

!-------------------------------------------------------------------------------
          case ('terrace','step','corner')            ! of select case(words(1))
!-------------------------------------------------------------------------------


            select case (parse_state)

              case(parse_state_ignore)
                ! ignore
                ! print *, 'warning ignoring line', line_number, buffer

              case(dissociation_id)

                dissociation_counter = dissociation_counter + 1

                i1 = get_index(words(1),    site_names)
                i2 = get_index(words(2),ads_site_names)
                i3 = get_index(words(3),    site_names)
                i4 = get_index(words(4),ads_site_names)
                i5 = get_index(words(5),    site_names)
                i6 = get_index(words(6),ads_site_names)

                if ( i1==0 .or. i2==0 .or. i3==0 .or. i4==0 .or. i5==0 .or. i6==0) &
                  call error_message(file_name, line_number, buffer, &
                             "wrong site name in the dissociation section")


                do i=1,dissociation_counter - 1
                  if ( dissociation_init%channels(i)%r      == current_reactant_id .and. &
                       dissociation_init%channels(i)%r_lst  == i1                  .and. &
                       dissociation_init%channels(i)%r_ast  == i2                  .and. &
                       dissociation_init%channels(i)%p1     == current_product1_id .and. &
                       dissociation_init%channels(i)%p1_lst == i3                  .and. &
                       dissociation_init%channels(i)%p1_ast == i4                  .and. &
                       dissociation_init%channels(i)%p2     == current_product2_id .and. &
                       dissociation_init%channels(i)%p2_lst == i5                  .and. &
                       dissociation_init%channels(i)%p2_ast == i6                      ) then
                    call error_message(file_name, line_number, buffer, "duplicated entry", stop = .false.)
                    duplicate_error = .true.
                  end if
                end do
                if (duplicate_error) stop 993

                ! check energy is defined for reactant's and for products' site_type and ads_site
                if( e_pars%ads_energy(current_reactant_id, i1, i2) == e_pars%undefined_energy .or. &
                    e_pars%ads_energy(current_product1_id, i3, i4) == e_pars%undefined_energy .or. &
                    e_pars%ads_energy(current_product2_id, i5, i6) == e_pars%undefined_energy ) then

                    call error_message(file_name, line_number, buffer, &
                                       " rate defined for site with undefined adsorption energy", &
                                       stop=.false., warning=.false.)

                    undefined_energy = .true.
                end if

                dissociation_init%channels(dissociation_counter)%r      = current_reactant_id
                dissociation_init%channels(dissociation_counter)%r_lst  = i1
                dissociation_init%channels(dissociation_counter)%r_ast  = i2
                dissociation_init%channels(dissociation_counter)%p1     = current_product1_id
                dissociation_init%channels(dissociation_counter)%p1_lst = i3
                dissociation_init%channels(dissociation_counter)%p1_ast = i4
                dissociation_init%channels(dissociation_counter)%p2     = current_product2_id
                dissociation_init%channels(dissociation_counter)%p2_lst = i5
                dissociation_init%channels(dissociation_counter)%p2_ast = i6

                select case (current_law_id)

                  case (Arrhenius_id)
                    if (nwords/=8) call error_message(file_name, line_number, buffer,&
                                              "Arrhenius must have 2 parameters")
                    read(words(7),*) pars(1)
                    read(words(8),*) pars(2)
                    dissociation_init%channels(dissociation_counter)%rate = arrhenius(c_pars%temperature, pars(1:2))

                  case (extArrhenius_id)
                    if (nwords/=9) call error_message(file_name, line_number, buffer,&
                                              "extArrhenius must have 3 parameters")
                    read(words(7),*) pars(1)
                    read(words(8),*) pars(2)
                    read(words(9),*) pars(3)
                    dissociation_init%channels(dissociation_counter)%rate = extArrhenius(c_pars%temperature, pars(1:3))

                  case default
                    call error_message(file_name, line_number, buffer, "Dissociation: this should not happen! Check the code!")

                end select

!                 ! Debugging printout
!                 print*, 'reaction: ', reaction_names(parse_state),&
!                        ' for species:', current_reactant_name, current_product1_name, current_product2_name
!                 print*, 'law: ', law_names(current_law_id),&
!                        ' from:', site_names(i1),ads_site_names(i2),&
!                        ' to:'  , site_names(i3),ads_site_names(i4),&
!                        ' and '  , site_names(i5),ads_site_names(i6)
!                print'(A,3f16.3)', 'with pars: ', pars
!                stop 112

              case default
                call error_message(file_name, line_number, buffer, "Dissociation: invalid site type statement")

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
              call error_message(file_name, line_number, buffer, "Dissociation: unknown key")

        end select                                       ! select case(words(1))

    end do ! while ios=0

    close(inp_unit)


    if (undefined_energy) then
      write(*, '(A)') ' Dissociation: error, rates defined for sites with undefined energies'
      stop 995

    else
      write(*, '(A)') ' Dissociation: passed check that energies are defined for all rates'

    end if

  end function dissociation_init




!-------------------------------------------------------------------------------
!  subroutine construct
!-------------------------------------------------------------------------------
  subroutine construct(this, ads, lat, e_pars, beta)

    class(dissociation_type), intent(inout) :: this
    integer, intent(in) :: ads
    class(mc_lat), intent(inout) :: lat
    class(energy_parameters), intent(in) :: e_pars
    real(dp), intent(in) :: beta

    integer :: id_r, id_p1, id_p2, m, iprocs, i_ast_p1, i_ast_p2
    integer :: row, col, lst, ast
    integer :: row_2, col_2, lst_2, ast_2
    integer :: n_channels
    integer :: i
    character(8) :: lst_name, lst_p2_name, ast_r_name, ast_p1_name, ast_p2_name

!debug(1) = (ads == 10)
!debug(1) = .false.

! debug printout header
!if (debug(1)) then
!  print *
!  print '(a)'    ,' Debug printout from rates_dissociation _class subroutine construct'
!  print '(a,i0)' ,' Print of rates as they are calculated.  ads =', ads
!  print '(a)'    ,' --------------------------------------------------------------------'
!  print '(a)'    ,'  ads# proc#  m  rate       lst      ast_r   ast_p1   lst_p2   ast_p2'
!  print '(a)'    ,' --------------------------------------------------------------------'
!
!end if

    ! Get the reactant information
    row  = lat%ads_list(ads)%row
    col  = lat%ads_list(ads)%col
    lst  = lat%lst(row,col)
    ast  = lat%ads_list(ads)%ast
    id_r = lat%ads_list(ads)%id

    ! Loop over possible positions for product 2
    n_channels=0
    do m=1,lat%n_nn(1)

      ! Get position and site type of neighbour m
      call lat%neighbor(ads, m, row_2, col_2)

      ! Check if the cell is free
      if (lat%occupations(row_2, col_2) == 0) then

        ! Get the lattice site type for product 2
        lst_2 = lat%lst(row_2, col_2)

        do iprocs=1, this%n_processes

          if ( this%channels(iprocs)%r     == id_r .and. &
               this%channels(iprocs)%r_lst == lst  .and. &
               this%channels(iprocs)%r_ast == ast  .and. &
               this%channels(iprocs)%p2_lst== lst_2       ) &
          then
            id_p1   = this%channels(iprocs)%p1
            id_p2   = this%channels(iprocs)%p2

            ! Determine all the values of the adsorption sites for products p1 and p2
            ! that match an entry in the process structure.
            ! If the rate for this entry is >0, add channel to the list
            do i_ast_p1=1, size(lat%avail_ads_sites(id_p1, lst  )%list)
            do i_ast_p2=1, size(lat%avail_ads_sites(id_p2, lst_2)%list)
              if ( this%channels(iprocs)%p1_ast == lat%avail_ads_sites(id_p1, lst  )%list(i_ast_p1) .and. &
                   this%channels(iprocs)%p2_ast == lat%avail_ads_sites(id_p2, lst_2)%list(i_ast_p2) .and. &
                   this%channels(iprocs)%rate > 0 )                                                        &
              then
                n_channels=n_channels + 1
                this%rate_info(ads)%list(n_channels)%proc  = iprocs
                this%rate_info(ads)%list(n_channels)%m     = m
                ! WARNING: Decide later if we need the rate field
                this%rate_info(ads)%list(n_channels)%rate  = this%channels(iprocs)%rate

!if(debug(1))then
!  lst_name    = site_names(lst)
!  lst_p2_name = site_names(lst_2)
!  ast_r_name  = ads_site_names(lat%avail_ads_sites(id_r , lst  )%list(ast))
!  ast_p1_name = ads_site_names(lat%avail_ads_sites(id_p1, lst  )%list(i_ast_p1))
!  ast_p2_name = ads_site_names(lat%avail_ads_sites(id_p2, lst_2)%list(i_ast_p2))
!
!  !  ads# chan#  m  rate       lst      ast_r   ast_p1   lst_p2   ast_p2
!  !   1    1     1__2.00e+00   terrace  top     fcc      terrace  fcc
!  !0        1         2         3         4        5        6         7
!  !1234567890123456789012345678901234567890123456890124567890123456789012345
!  print '(t4,i0, t9,i0, t15,i0, t16,1pe10.2, t29,A7, 2x,A3, t46,A3, t55,A7, 2x,A3)', &
!        ads, iprocs, m, this%channels(iprocs)%rate, lst_name, ast_r_name, ast_p1_name, lst_p2_name, ast_p2_name
!end if

              end if
            end do
            end do

            ! save the number of channels in the rate_info structure
            this%rate_info(ads)%n_channels = n_channels

          end if

        end do

      end if ! occupations

    end do ! m

 end subroutine construct

!------------------------------------------------------------------------------
  subroutine print(this, c_pars)
!------------------------------------------------------------------------------
    class(dissociation_type), intent(in) :: this

    class(control_parameters), intent(in) :: c_pars

    integer :: i
    character(:), allocatable :: reactant_name, product1_name, product2_name

    print*
    print*, 'Dissociation Rates:'
    print '(A)', ' ---------------------------'
    do i=1,size(this%channels,1)
      reactant_name = c_pars%ads_names(this%channels(i)%r )
      product1_name = c_pars%ads_names(this%channels(i)%p1)
      product2_name = c_pars%ads_names(this%channels(i)%p2)
      write(*,'(1x,A,A,A,2X,A,A,2X,A,A,e12.3)')  &
              trim(reactant_name)// ' -> '// trim(product1_name)// ' + '// trim(product2_name)//': ',&
              site_names(this%channels(i)%r_lst),  ads_site_names(this%channels(i)%r_ast), &
              site_names(this%channels(i)%p1_lst), ads_site_names(this%channels(i)%p1_ast), &
              site_names(this%channels(i)%p2_lst), ads_site_names(this%channels(i)%p2_ast), &
              this%channels(i)%rate
    end do
    print*

  end subroutine print

!-----------------------------------------------------------------------------
!             Temperature dependence law subroutines
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  real(dp) function arrhenius(temperature, parameters)
!-----------------------------------------------------------------------------
    real(dp), intent(in) :: temperature
    real(dp), dimension(:), intent(in) :: parameters
    real(dp) :: prefactor, act_energy

    prefactor  = parameters(1)
    act_energy = parameters(2)

    arrhenius = prefactor*exp(-act_energy/(kB*temperature))

  end function arrhenius

!-----------------------------------------------------------------------------
  real(dp) function extArrhenius(temperature, parameters)
!-----------------------------------------------------------------------------
    real(dp), intent(in) :: temperature
    real(dp), dimension(:), intent(in) :: parameters
    real(dp) :: prefactor, act_energy, power, kT

      prefactor  = parameters(1)
      act_energy = parameters(2)
      power = parameters(3)
      kT = kB*temperature

      extArrhenius = prefactor*exp(-act_energy/kT)/(kT**power)

  end function extArrhenius

end module rates_dissociation_class
