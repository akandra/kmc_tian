module rates_desorption_class

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
  public    :: desorption_init, desorption_type

  type :: desorption_type

    logical :: is_defined = .false.
    ! Desorption Rates
    !                          n_adsorbate              -> which particle
    !                          .
    real(dp), dimension(:), allocatable :: rates

    !               (   n_species                       -> which species
    !                   .  n_site_type                  -> where from
    !                   .  .  n_adsorption_sites )
    !                   .  .  .
    real(dp), dimension(:, :, :), allocatable :: process

  contains
    procedure :: construct
    procedure :: print

  end type


contains
!------------------------------------------------------------------------------
  function desorption_init(c_pars, lat, e_pars)
!------------------------------------------------------------------------------
    type(desorption_type) desorption_init

    type(control_parameters), intent(in)    :: c_pars
    type(mc_lat)            , intent(in)    :: lat
    type(energy_parameters) , intent(in)    :: e_pars

    integer :: ios, nwords, line_number, i1, i2

    integer :: species, st1, ast1
    logical :: e_defined1, r_defined, undefined_rate, undefined_energy

    character(len=max_string_length)                :: buffer
    character(len=max_string_length)                :: words(100)
    character(len=len(trim(c_pars%rate_file_name))) :: file_name

    character(len=10)     :: current_species_name
    integer               :: current_species_id
    integer               :: current_law_id

    integer               :: parse_state
    integer, parameter    :: parse_state_ignore     = -1
    integer, parameter    :: parse_state_default    =  0
    integer, parameter    :: parse_state_desorption =  desorption_id


    real(dp), dimension(3):: pars = 0.0_dp

    integer,  parameter   :: default_int = 0
    real(dp), parameter   :: default_rate  = -1.0_dp


    ! Allocate rates array
    allocate( desorption_init%rates(lat%n_rows*lat%n_cols) )
    allocate( desorption_init%process( c_pars%n_species,&
                                 n_max_lat_site_types, n_max_ads_sites) )

    desorption_init%process = default_rate
    desorption_init%rates    = 0.0_dp

    !  read rate definitions from the input file
    file_name = c_pars%rate_file_name
    call open_for_read(inp_unit, file_name )

    ios = 0
    parse_state = parse_state_default
    line_number = 0
    undefined_energy = .false.



    do while (ios == 0)

      read(inp_unit, '(A)', iostat=ios) buffer
      line_number = line_number + 1
        ! ios < 0: end of record condition encountered or endfile condition detected
        ! ios > 0: an error is detected
        ! ios = 0  otherwise

      if (ios /= 0) exit

        ! Split an input string
        words = ''
        call split_string(buffer, words, nwords)

        select case (words(1)) ! take a keyword
!------------------------------------------------------------------------------
          case('desorption')                               ! select case (words(1)
!------------------------------------------------------------------------------
            desorption_init%is_defined = .true.

            if (parse_state /= parse_state_default) &
              call error_message(file_name, line_number, buffer, &
                         "invalid ending of the reaction section")
            parse_state = parse_state_desorption
            if (nwords/=3) call error_message(file_name, line_number, buffer, &
                               "desorption key must have 2 parameters")

            read(words(2),'(A)') current_species_name
            current_species_id = get_index(current_species_name, c_pars%ads_names )
            if (current_species_id == 0) call error_message(file_name, line_number, buffer, &
                                                  "inconsistent desorption definition")

            current_law_id = get_index(words(3), rct_law_names )
            if (current_law_id == 0) call error_message(file_name, line_number, buffer, &
                                                  "unknown temperature law")
!            print*, 'name     =', current_species_name
!            print*, 'id       =', current_species_id
!            print*, c_pars%ads_names

!------------------------------------------------------------------------------
          case ('terrace','step','corner')              ! select case(words(1))
!------------------------------------------------------------------------------


            select case (parse_state)

              case(parse_state_ignore)
                ! ignore
                ! print *, 'warning ignoring line', line_number, buffer

              case(desorption_id)

                i1 = get_index(words(1),lat_site_names)
                i2 = get_index(words(2),ads_site_names)

                if ( i1==0 .or. i2==0) &
                  call error_message(file_name, line_number, buffer, &
                             "wrong site name in the desorption section")

                ! check for duplicate entry
                if (desorption_init%process(current_species_id,i1,i2 ) /= default_rate)&
                  call error_message(file_name, line_number, buffer, "duplicated entry")

                ! check energy is defined for initial site_type and ads_site
                if( e_pars%ads_energy(current_species_id, i1, i2) == e_pars%undefined_energy) then

                    call error_message(file_name, line_number, buffer, &
                                       "rate defined for site with undefined adsorption energy", &
                                       stop=.false., warning=.false.)
                    undefined_energy = .true.
                end if

                select case (current_law_id)

                  case (Arrhenius_id)
                    if (nwords/=4) call error_message(file_name, line_number, buffer,&
                                              "Arrhenius must have 2 parameters")
                    read(words(3),*) pars(1)
                    read(words(4),*) pars(2)
                    desorption_init%process(current_species_id,i1,i2 ) = &
                                arrhenius(c_pars%temperature, pars(1:2))

                  case (extArrhenius_id)
                    if (nwords/=5) call error_message(file_name, line_number, buffer,&
                                              "extArrhenius must have 3 parameters")
                    read(words(3),*) pars(1)
                    read(words(4),*) pars(2)
                    read(words(5),*) pars(3)
                    desorption_init%process(current_species_id,i1,i2 ) = &
                                extArrhenius(c_pars%temperature, pars(1:3))

                  case default
                    call error_message(file_name, line_number, buffer, "This should not happen! Check the code!")

                end select

!                 print*, 'reaction: ', reaction_names(parse_state),&
!                        ' for species:', current_species_name
!                 print*, 'law: ', rct_law_names(current_law_id),&
!                        ' from:', lat_site_names(i1),ads_site_names(i2)
!                print'(A,3f16.3)', 'with pars: ', pars

              case default
                call error_message(file_name, line_number, buffer, "invalid site type statement")

            end select

!------------------------------------------------------------------------------
          case('')                                      ! select case(words(1))
!------------------------------------------------------------------------------
            if (buffer == '') then
              parse_state = parse_state_default
!              print*, 'blank line '
!            else
!              print*, 'comment: ', trim(buffer)
            end if

!------------------------------------------------------------------------------
          case default                                  ! select case(words(1))
!------------------------------------------------------------------------------
            if ( parse_state == parse_state_default .and. get_index(words(1),reaction_names) /= 0 ) &
              parse_state = parse_state_ignore

            if (parse_state /= parse_state_ignore) &
              call error_message(file_name, line_number, buffer, "unknown key")

        end select                                      ! select case(words(1))

    end do ! while ios=0

    close(inp_unit)

    if (.not. undefined_energy) then

      if (desorption_init%is_defined) then
        write(6, '(A/)') ' desorption: passed check that energies are defined for all rates'
      else
        write(6, '(A)') ' no desorption reactions'
      end if

    end if


    ! ---------------------------------------------------------------------------------------------
    ! Check the input consistency
    ! ---------------------------------------------------------------------------------------------
    !
    ! Check if rates are defined for all values of site_types, ads_sites
    !   for which adsorption energies are defined
    ! Adsorpton energies (n_species x n_site_type x n_adsorption_sites)
    ! real(dp), dimension(:,:,:), allocatable :: ads_energy
    !
    ! Note ads_energies could be allocated of n_site_types rather than n_max_lat_site_types
    !

    if (desorption_init%is_defined) then

      undefined_rate = .false.
      do species   = 1, c_pars%n_species
      do st1       = 1, n_max_lat_site_types
      do ast1      = 1, n_max_ads_sites

        e_defined1 = e_pars%ads_energy(species, st1, ast1) /= e_pars%undefined_energy
        r_defined  = desorption_init%process (species, st1, ast1) /= default_rate

        if ( e_defined1 .and. (.not. r_defined)) then
          if (.not. undefined_rate) then
            undefined_rate = .true.
            print '(2A)',  ' Desorption: warning missing rate definitions in the file ', file_name
            print '(A)',  ' Missing definitions:'
            print '(6x, A)', 'ads  lat_site    ads_site'
          end if

          print '(6x, a5, A10, 2x, a3, 6x, a10, 2x, a3, 6x, L1)' ,            &
                  c_pars%ads_names(species),             &
                  lat_site_names(st1), ads_site_names(ast1)
        end if

      end do
      end do
      end do

      if(undefined_rate) then
        print '(A)', ' Add rate definitions if they should be present'
        print *
      else
        print '(A)', ' Desorption: passed rates consistency check'
        print *
      end if

      if(undefined_energy) then
        print '(A)', ' Desorption: error - undefined energy -- aborting execution'
        stop 997
      end if

    end if

  end function desorption_init

!-----------------------------------------------------------------------------
  subroutine construct(this, ads, lat, e_pars, beta)
!-----------------------------------------------------------------------------
    class(desorption_type), intent(inout) :: this
    integer, intent(in) :: ads
    class(mc_lat), intent(inout) :: lat
    class(energy_parameters), intent(in) :: e_pars
    real(dp), intent(in) :: beta

    integer :: id
    integer :: row, col, lst, ast

    row = lat%ads_list(ads)%row
    col = lat%ads_list(ads)%col
    lst = lat%lst(row,col)
    ast = lat%ads_list(ads)%ast
    id  = lat%ads_list(ads)%id

    this%rates(ads) = this%process(id, lst, ast)&
                    *exp( beta*( energy(ads,lat,e_pars) - e_pars%ads_energy(id,lst,ast) ))
!   print*
!   print*, 'id ',id, ' lst ',lst,' ads. site ', ast
!   print*, 'rate ',this%process(id, lst, ast)

!   print '(A,i4)', 'ads ',ads
!   print '(A,e18.4)',' Energy = ', energy(ads,lat,e_pars)
!   print *,' rate = ', this%rates(ads)
!   write(*,*) 'pause'
!   read(*,*)


  end subroutine construct

!------------------------------------------------------------------------------
  subroutine print(this, c_pars)
!------------------------------------------------------------------------------
    class(desorption_type), intent(in) :: this

    class(control_parameters), intent(in) :: c_pars

    integer :: i, i1, i2

    print*, 'Desorption Rates:'
    do i=1,c_pars%n_species
      print '(/A)','---------------------------'
      print '( A,A)','species: ', c_pars%ads_names(i)
      print '(A)', '---------------------------'
      do i1=1,n_max_lat_site_types
      do i2=1,n_max_ads_sites
        if (this%process(i,i1,i2)< 0.0_dp) then
          cycle
        else
          write(*,'(A,A,1pe11.2)') &
              lat_site_names(i1), ads_site_names(i2), this%process(i,i1,i2)
        end if
      end do
      end do
    end do
    print*

  end subroutine print

end module rates_desorption_class
