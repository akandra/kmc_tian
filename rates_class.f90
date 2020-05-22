module rates_class

  use constants
  use control_parameters_class
  use open_file
  use utilities

  implicit none

  type, public :: rates_type

    ! Hopping rates (n_species x n_site_type x n_adsorption_sites
    !                               x n_site_type x n_adsorption_sites)
    real(dp), dimension(:,:,:,:,:), allocatable :: r_hop
    !
    ! To Be Added
    ! Reaction rates

  contains

!    procedure :: read  => control_parameters_read

  end type


  interface rates_type

    module procedure :: rates_init

  end interface

contains

  function rates_init(control_pars)

    type(rates_type) rates_init

    type(control_parameters), intent(inout) :: control_pars

     integer :: i, ios, nwords, line_number, i1, i2, i3, i4
    character(len=max_string_length) :: buffer
    character(len=max_string_length) :: words(100)
    character(len=len(trim(control_pars%rate_file_name))) :: file_name

    character(len=10) :: current_species_name
    integer           :: current_species_id
    character(len=20) :: current_law_name
    integer           :: current_law_id

    integer           :: parse_state
    integer           :: parse_state_default      =0

    real(dp), dimension(3) :: pars = 0.0_dp

    integer,  parameter:: default_int = 0
    real(dp), parameter:: default_dp  = -1.0_dp

    allocate(rates_init%r_hop( control_pars%n_species,&
                               n_site_types, n_ads_sites,&
                               n_site_types, n_ads_sites) )

    rates_init%r_hop = default_dp

    !  read rate definitions from the input file
    file_name = control_pars%rate_file_name
    call open_for_read(inp_unit, file_name )

    ios = 0
    parse_state = parse_state_default
    line_number = 0

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

          case('hopping')
            if (parse_state /= parse_state_default) &
              call error(file_name, line_number, buffer, &
                         "invalid ending of the reaction section")
            parse_state = get_index('hopping',reaction_names)
            if (nwords/=3) call error(file_name, line_number, buffer, &
                               "hopping key must have 2 parameters")

            read(words(2),'(A)') current_species_name
            current_species_id = get_index(current_species_name, control_pars%ads_names )
            if (current_species_id == 0) call error(file_name, line_number, buffer, &
                                                  "inconsistent hopping definition")

            current_law_id = get_index(words(3), law_names )
            if (current_law_id == 0) call error(file_name, line_number, buffer, &
                                                  "unknown temperature law")
            print*, 'name     =', current_species_name
            print*, 'id       =', current_species_id
            print*, control_pars%ads_names

          case ('terrace','step','corner')

            select case (parse_state)

              case(hopping_id)

                i1 = get_index(words(1),    site_names)
                i2 = get_index(words(2),ads_site_names)
                i3 = get_index(words(3),    site_names)
                i4 = get_index(words(4),ads_site_names)

                if ( i1==0 .or. i2==0 .or. i3==0 .or. i4==0) &
                  call error(file_name, line_number, buffer, &
                             "wrong species name in the hopping section")
                if (rates_init%r_hop(current_species_id,i1,i2,i3,i4 ) /= default_dp)&
                  call error(file_name, line_number, buffer, "duplicated entry")

                select case (current_law_id)

                  case (Arrhenius_id)
                    if (nwords/=6) call error(file_name, line_number, buffer,&
                                              "Arrhenius must have 2 parameters")
                    read(words(5),*) pars(1)
                    read(words(6),*) pars(2)
                    rates_init%r_hop(current_species_id,i1,i2,i3,i4 ) = &
                                0.0_dp ! arrhenius(control_pars%temperature, pars(1:2))
                    rates_init%r_hop(current_species_id,i3,i4,i1,i2 ) = &
                    rates_init%r_hop(current_species_id,i1,i2,i3,i4 )

                  case (extArrhenius_id)
                    if (nwords/=7) call error(file_name, line_number, buffer,&
                                              "extArrhenius must have 3 parameters")
                    read(words(5),*) pars(1)
                    read(words(6),*) pars(2)
                    read(words(7),*) pars(3)
                    rates_init%r_hop(current_species_id,i1,i2,i3,i4 ) = &
                                0.0_dp ! extarrhenius(control_pars%temperature, pars(1:3))
                    rates_init%r_hop(current_species_id,i3,i4,i1,i2 ) = &
                    rates_init%r_hop(current_species_id,i1,i2,i3,i4 )

                  case default
                    call error(file_name, line_number, buffer, "This cannot happen! Check the code!")

                end select

                 print*, 'reaction: ', reaction_names(parse_state),&
                        ' for species:', current_species_name
                 print*, 'law: ', law_names(current_law_id),&
                        ' from:', site_names(i1),ads_site_names(i2),&
                        ' to:'  , site_names(i3),ads_site_names(i4)
                print'(A,3f16.3)', 'with pars: ', pars

              case default
                call error(file_name, line_number, buffer, "invalid site type statement")

            end select


          case('')
            if (buffer == '') then
              parse_state = parse_state_default
              print*, 'blank line '
            else
              print*, 'comment: ', trim(buffer)
            end if

          case default
            print*, 'unprocessed line: ', trim(buffer)
            call error(file_name, line_number, buffer, "unknown key")

        end select

    end do ! while ios=0

    ! Check the input consistency

  end function

!real(8) function arrhenius(temperature, prefactor, act_energy)
!
!real(8), intent(in) :: temperature, prefactor, act_energy
!
!    arrhenius = prefactor*exp(-act_energy/temperature)
!
!
!end function


end module rates_class
