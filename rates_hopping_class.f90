module rates_hopping_class

  use constants
  use control_parameters_class
  use mc_lat_class
  use energy_parameters_class
  use energy_mod
  use open_file
  use utilities

  implicit none

  private
  public    :: hopping_rates_init, hopping_rates_type

  type :: v_list_dp

    real(dp), dimension(:), allocatable :: list ! n_avail_ads_sites

  end type

  type :: hopping_rates_type

    ! Hopping Rates
    !                          n_adsorbate              -> which particle
    !                          .  n_neighbor            -> where to
    !                          .  .
    type(v_list_dp), dimension(:, :), allocatable :: rates

    !               (   n_species                       -> which species
    !                   .  n_site_type                  -> where from
    !                   .  .  n_adsorption_sites
    !                   .  .  .  n_site_type            -> where to
    !                   .  .  .  .  n_adsorption_sites )
    !                   .  .  .  .  .
    real(dp), dimension(:, :, :, :, :), allocatable :: process

  contains
    procedure :: construct
    procedure ::  print

  end type


  interface hopping_rates_type

    module procedure :: hopping_rates_init

  end interface

contains
!------------------------------------------------------------------------------
  function hopping_rates_init(c_pars, lat, e_pars)
!------------------------------------------------------------------------------
    type(hopping_rates_type) hopping_rates_init

    type(control_parameters), intent(inout) :: c_pars
    type(mc_lat)            , intent(in)    :: lat
    type(energy_parameters) , intent(in)    :: e_pars

    integer :: i, ios, nwords, line_number, i1, i2, i3, i4, m

    integer :: species, st1, st2, ast1, ast2
    logical :: e_defined1, e_defined2, r_defined, undefined_rate, undefined_energy

    character(len=max_string_length)                :: buffer
    character(len=max_string_length)                :: words(100)
    character(len=len(trim(c_pars%rate_file_name))) :: file_name
    character(1)                                    :: answer

    character(len=10)     :: current_species_name
    integer               :: current_species_id
    character(len=20)     :: current_law_name
    integer               :: current_law_id

    integer               :: parse_state
    integer, parameter    :: parse_state_ignore  = -1
    integer, parameter    :: parse_state_default =  0
    integer, parameter    :: parse_state_hopping =  1


    real(dp), dimension(3):: pars = 0.0_dp

    integer,  parameter   :: default_int = 0
    real(dp), parameter   :: default_rate  = -1.0_dp

    integer :: row, col, site, id

    ! Allocate rates array
    allocate( hopping_rates_init%rates(lat%n_rows*lat%n_cols,lat%n_nn(1)) )
    do i=1,lat%n_ads_tot()
      row = lat%ads_list(i)%row
      col = lat%ads_list(i)%col
      site = lat%lst(row,col)
      id  = lat%ads_list(i)%id
      do m=1,lat%n_nn(1)
        allocate( hopping_rates_init%rates(i,m)%list(size(lat%avail_ads_sites(id,site)%list)) )
      end do
    end do

    allocate(hopping_rates_init%process( c_pars%n_species,&
                                 n_max_site_types, n_max_ads_sites,&
                                 n_max_site_types, n_max_ads_sites) )

    hopping_rates_init%process = default_rate

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
          case('hopping')                               ! select case (words(1)
!------------------------------------------------------------------------------
            if (parse_state /= parse_state_default) &
              call error_message(file_name, line_number, buffer, &
                         "invalid ending of the reaction section")
            parse_state = parse_state_hopping
            if (nwords/=3) call error_message(file_name, line_number, buffer, &
                               "hopping key must have 2 parameters")

            read(words(2),'(A)') current_species_name
            current_species_id = get_index(current_species_name, c_pars%ads_names )
            if (current_species_id == 0) call error_message(file_name, line_number, buffer, &
                                                  "inconsistent hopping definition")

            current_law_id = get_index(words(3), law_names )
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
                 print *, 'warning ignoring line', line_number, buffer

              case(hopping_id)

                i1 = get_index(words(1),    site_names)
                i2 = get_index(words(2),ads_site_names)
                i3 = get_index(words(3),    site_names)
                i4 = get_index(words(4),ads_site_names)

                if ( i1==0 .or. i2==0 .or. i3==0 .or. i4==0) &
                  call error_message(file_name, line_number, buffer, &
                             "wrong species name in the hopping section")

                ! check for duplicate entry
                if (hopping_rates_init%process(current_species_id,i1,i2,i3,i4 ) /= default_rate)&
                  call error_message(file_name, line_number, buffer, "duplicated entry (check symetry duplicates)")

                ! check energy is defined for initial and final site_type and ads_site
                if( e_pars%ads_energy(current_species_id, i1, i2) == e_pars%undefined_energy .or. &
                    e_pars%ads_energy(current_species_id, i3, i4) == e_pars%undefined_energy ) then

                    call error_message(file_name, line_number, buffer, &
                                       "rate defined for site with undefined adsorption energy", &
                                       stop=.false., warning=.true.)

                    undefined_energy = .true.
                end if

                select case (current_law_id)

                  case (Arrhenius_id)
                    if (nwords/=6) call error_message(file_name, line_number, buffer,&
                                              "Arrhenius must have 2 parameters")
                    read(words(5),*) pars(1)
                    read(words(6),*) pars(2)
                    hopping_rates_init%process(current_species_id,i1,i2,i3,i4 ) = &
                                arrhenius(c_pars%temperature, pars(1:2))
                    ! symmetrize
                    hopping_rates_init%process(current_species_id,i3,i4,i1,i2 ) = &
                    hopping_rates_init%process(current_species_id,i1,i2,i3,i4 )

                  case (extArrhenius_id)
                    if (nwords/=7) call error_message(file_name, line_number, buffer,&
                                              "extArrhenius must have 3 parameters")
                    read(words(5),*) pars(1)
                    read(words(6),*) pars(2)
                    read(words(7),*) pars(3)
                    hopping_rates_init%process(current_species_id,i1,i2,i3,i4 ) = &
                                extArrhenius(c_pars%temperature, pars(1:3))
                    ! symmetrize
                    hopping_rates_init%process(current_species_id,i3,i4,i1,i2 ) = &
                    hopping_rates_init%process(current_species_id,i1,i2,i3,i4 )

                  case default
                    call error_message(file_name, line_number, buffer, "This should not happen! Check the code!")

                end select

!                 print*, 'reaction: ', reaction_names(parse_state),&
!                        ' for species:', current_species_name
!                 print*, 'law: ', law_names(current_law_id),&
!                        ' from:', site_names(i1),ads_site_names(i2),&
!                        ' to:'  , site_names(i3),ads_site_names(i4)
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

    if (undefined_energy) then
      print *
      write(*, '(A)') 'warnings issued because of extraneous lines in rates file'

    else
      print *
      write(*, '(A)') 'passed check that energies are defined for all rates'
!      pause
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
    ! Note ads_energies could be allocated of n_site_types rather than n_max_site_types
    !

    undefined_rate = .false.
    do species   = 1, c_pars%n_species
    do st1       = 1, n_max_site_types
    do ast1      = 1, n_max_ads_sites
    do st2       = st1, n_max_site_types
    do ast2      = 1, n_max_ads_sites

      e_defined1 = e_pars%ads_energy(species, st1, ast1) /= e_pars%undefined_energy
      e_defined2 = e_pars%ads_energy(species, st2, ast2) /= e_pars%undefined_energy
      r_defined  = hopping_rates_init%process (species, st1, ast1, st2, ast2) /= default_rate

      if ( (e_defined1 .and. e_defined2) .and. (.not. r_defined)) then
        if (.not. undefined_rate) then
          undefined_rate = .true.
          print*
          print '(A)',  '--- Dear Sir, Madam:'
          print '(A)',  '      It is my duty to inform you that there are missing rate definitions in'
          print '(2A)', '      the file ', file_name
          print *
          print '(A)',  '      Missing definitions:'
          print*
          !             123451234567890xx123xxxxxx1234567890xx123xxxxxx1xxxxxx1
          print '(6x, A)', 'ads  lat_site    ads_site lat_site    ads_site'

        end if

        print '(6x, a5, A10, 2x, a3, 6x, a10, 2x, a3, 6x, L1, 7x, L1, 7x, L1)' ,            &
                c_pars%ads_names(species),             &
                site_names(st1), ads_site_names(ast1), &
                site_names(st2), ads_site_names(ast2)
                !e_defined1, e_defined2, r_defined
      end if

    end do
    end do
    end do
    end do
    end do

!    ! check for extraneous rates (rate defined but energy not defined)
!    undefined_energy = .false.
!    do species   = 1, c_pars%n_species
!    do st1       = 1, n_max_site_types
!    do ast1      = 1, n_max_ads_sites
!    do st2       = st1, n_max_site_types  ! check: do we have to explicitly check lower part of matrix?
!    do ast2      = 1, n_max_ads_sites
!
!      e_defined1 = e_pars%ads_energy(species, st1, ast1) /= e_pars%undefined_energy
!      e_defined2 = e_pars%ads_energy(species, st2, ast2) /= e_pars%undefined_energy
!      r_defined  = hopping_rates_init%process (species, st1, ast1, st2, ast2) /= default_rate
!
!      if ( r_defined .and. .not. (e_defined1 .and. e_defined2)) then
!
!        if (.not. undefined_energy) then
!          undefined_energy = .true.
!          print*
!          print '(A)', 'Dear Sir, Madam extraneous rate definition error message'
!          print*
!          !             123451234567890xx123xxxxxx1234567890xx123xxxxxx1xxxxxx1
!          print '(A)', 'ads  lat_site    ads_site lat_site    ads_site'
!        end if
!
!        print '(a5, A10, 2x, a3, 6x, a10, 2x, a3, 6x, L1, 7x, L1, 7x, L1)' ,            &
!                c_pars%ads_names(species),             &
!                site_names(st1), ads_site_names(ast1), &
!                site_names(st2), ads_site_names(ast2)
!                !e_defined1, e_defined2, r_defined
!      end if
!
!
!    end do
!    end do
!    end do
!    end do
!    end do

    if(undefined_rate) then
      print '(/6x, A)', 'Please supply the required rates'
      print '(/6x, A)', 'As always, I remain your humble servant, kMC Code'
      print *
      stop 997

    else
      print '(/A)', 'passed requited rates consistency check'
      print*
!      stop 'debugging stop'

    end if

  end function hopping_rates_init

!-----------------------------------------------------------------------------
  subroutine construct(this, ads, lat, e_pars, beta)
!-----------------------------------------------------------------------------
    class(hopping_rates_type), intent(inout) :: this
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

    ! Delete particle ads from the old position
    ! we do it here since we never work with occupations inside the following loop
    lat%occupations(row_old,col_old) = 0

    ! Loop over possible new positions of particle ads
    do m=1,lat%n_nn(1)

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

          ! Apply detailed balance when
          ! energy in the old position < energy in the new position
          if (energy_old < energy_new) then
              this%rates(ads,m)%list(iads) = this%process(id, lst_old, ast_old, lst_new, ast_new)&
                  *exp( -beta*(energy_new - energy_old) )
!                print*
!                print*, 'id ',id, ' old site ',lst_old,' old ads. site ', ast_old
!                print*, ' new site ',lst_new,' new ads. site ', ast_new
!                print*, 'rate ',this%process(id, lst_old, ast_old, lst_new, ast_new)
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

  end subroutine construct

!------------------------------------------------------------------------------
  subroutine print(this, c_pars)
!------------------------------------------------------------------------------
    class(hopping_rates_type), intent(in) :: this

    class(control_parameters), intent(in) :: c_pars

    integer :: i, i1, i2, i3, i4

    print*, 'Hopping Rates:'
    do i=1,size(this%process,1)
      print '(/A)','---------------------------'
      print '( A,A)','species: ', c_pars%ads_names(i)
      print '(A)', '---------------------------'
      do i1=1,n_max_site_types
      do i2=1,n_max_ads_sites
      do i3=1,n_max_site_types
      do i4=1,n_max_ads_sites
        if (this%process(i,i1,i2,i3,i4)< 0.0_dp) then
          cycle
        else
          write(*,'(A,A,2X,A,A,6e12.3)') &
              site_names(i1), ads_site_names(i2), &
              site_names(i3), ads_site_names(i4), &
              this%process(i,i1,i2,i3,i4)
        end if
      end do
      end do
      end do
      end do
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

end module rates_hopping_class
