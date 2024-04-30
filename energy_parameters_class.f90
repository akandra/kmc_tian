module energy_parameters_class

  use utilities
  use control_parameters_class
  use open_file

  implicit none

  private
  public :: energy_parameters_init

  type, public :: energy_parameters

    ! Adsorpton energies (n_species x n_site_type x n_adsorption_sites)
    real(dp), dimension(:,:,:), allocatable :: ads_energy
    ! Interaction energy law id (n_species x n_site_type x n_adsorption_sites x n_species x n_site_type x n_adsorption_sites)
    integer, dimension(:,:,:,:,:,:), allocatable :: int_energy_law_id
    ! Interaction energy (n_species x n_site_type x n_adsorption_sites x n_species x n_site_type x n_adsorption_sites x n_shells )
    real(dp), dimension(:,:,:,:,:,:,:),allocatable :: int_energy_pars
    ! Interaction energy mask (n_species x n_species x n_shells)
    ! .true. means to skip interaction
    logical, dimension(:,:,:,:,:,:,:),allocatable :: int_energy_skip
    !
    logical :: is_interaction
    ! Default value for undefined energy
    real(dp) :: undefined_energy


  contains

!    procedure :: read  => control_parameters_read

  end type


contains

  function energy_parameters_init(control_pars)

    type(energy_parameters) energy_parameters_init

    type(control_parameters), intent(inout) :: control_pars

    integer :: i, ios, nwords, line_number, i1, i1s, i1a, i2, i2s, i2a, i_law

    character(len=max_string_length) :: buffer
    character(len=max_string_length) :: words(100)
    character(len=len(trim(control_pars%energy_file_name))) :: file_name

    character(len=10) :: current_species_name
    integer           :: current_species_id

    integer           :: parse_state
    integer           :: parse_state_default      =0
    integer           :: parse_state_adsorption   =1
    integer           :: parse_state_interaction  =2


    integer,  parameter:: default_int = 0
    real(dp), parameter:: default_dp  = huge(0.0_dp)

    i = control_pars%n_species
    allocate(energy_parameters_init%ads_energy(i,n_max_lat_site_types,n_max_ads_sites))
    allocate(energy_parameters_init%int_energy_law_id(i,n_max_lat_site_types,n_max_ads_sites,&
                                                      i,n_max_lat_site_types,n_max_ads_sites))
    allocate(energy_parameters_init%int_energy_pars(i,n_max_lat_site_types,n_max_ads_sites,&
                                                    i,n_max_lat_site_types,n_max_ads_sites,n_shells))
    allocate(energy_parameters_init%int_energy_skip(i,n_max_lat_site_types,n_max_ads_sites,&
                                                    i,n_max_lat_site_types,n_max_ads_sites,n_shells))

    energy_parameters_init%ads_energy = default_dp
    energy_parameters_init%int_energy_law_id = default_int
    energy_parameters_init%int_energy_pars = default_dp
    energy_parameters_init%int_energy_skip = .true.
    energy_parameters_init%is_interaction  = .false.
    energy_parameters_init%undefined_energy = default_dp

    !  read energy definitions from the input file
    file_name = control_pars%energy_file_name
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

        if (words(1) =='adsorption') then

          if (parse_state /= parse_state_default) &
            call error_message(file_name, line_number, buffer, &
                       "invalid ending of the adsorption/interaction section")
          parse_state = parse_state_adsorption
          if (nwords/=2) call error_message(file_name, line_number, buffer, &
                             "adsorption key must have 1 parameter")

          read(words(2),'(A)') current_species_name
          current_species_id = get_index(current_species_name, control_pars%ads_names )
          if (current_species_id == 0) call error_message(file_name, line_number, buffer, &
                                                "inconsistent adsorbate definition")
!            print*, 'name     =', current_species_name
!            print*, 'id       =', current_species_id
!            print*, control_pars%ads_names

        elseif (get_index(words(1),lat_site_names) /= 0) then

          if (parse_state /= parse_state_adsorption) &
            call error_message(file_name, line_number, buffer, "invalid site type statement")

          i1 = get_index(words(1),lat_site_names)
          i2 = get_index(words(2),ads_site_names)

          if ( i1==0 ) &
              call error_message(file_name, line_number, buffer, &
                           "unknown lattice site type in the adsorption section")
          if ( i2==0 ) &
              call error_message(file_name, line_number, buffer, &
                           "unknown adsorption site type in the adsorption section")
          if (energy_parameters_init%ads_energy(current_species_id,i1,i2 ) /= default_dp)&
              call error_message(file_name, line_number, buffer, "duplicated entry")

          read(words(3),*) energy_parameters_init%ads_energy(current_species_id,i1,i2 )
!              print*, 'species  ' ,current_species_id, &
!                      'site     ' ,get_index(words(1),lat_site_names),&
!                      'ads_site ' ,get_index(words(2),ads_site_names)
!              print*, 'energy:  ' ,energy_parameters_init%ads_energy&
!                                    (current_species_id,&
!                                     get_index(words(1),lat_site_names),&
!                                     get_index(words(2),ads_site_names) )

        elseif (words(1) == 'interaction') then

          if (parse_state /= parse_state_default) &
            call error_message(file_name, line_number, buffer, &
                       "invalid ending of the adsorption/interaction section")
          parse_state = parse_state_interaction
          if (nwords/=1) call error_message(file_name, line_number, buffer, &
                             "interaction key must have no parameters")

        elseif (get_index(words(1),int_law_names) /= 0) then

          if (parse_state /= parse_state_interaction) &
            call error_message(file_name, line_number, buffer, &
                      "invalid interaction law statement")
          if (nwords/=7+n_shells) call error_message(file_name, line_number, buffer, &
                            "interaction law key must have (7 + n_shells) parameters")

          i_law  = get_index(words(1),int_law_names)
          i1  = get_index(words(2),control_pars%ads_names)
          i1s = get_index(words(3),lat_site_names)
          i1a = get_index(words(4),ads_site_names)
          i2  = get_index(words(5),control_pars%ads_names)
          i2s = get_index(words(6),lat_site_names)
          i2a = get_index(words(7),ads_site_names)

          if ( i1==0 .or. i2==0 ) call error_message(file_name, line_number, buffer, &
                            "unknown species name in the interaction law")
          if ( i1s==0 .or. i2s==0 ) &
              call error_message(file_name, line_number, buffer, &
                           "unknown lattice site type in the interaction section")
          if ( i1a==0 .or. i2a==0 ) &
              call error_message(file_name, line_number, buffer, &
                           "unknown adsorption site type in the interaction section")
          if (energy_parameters_init%int_energy_law_id(i1,i1s,i1a,i2,i2s,i2a) /= default_int)&
              call error_message(file_name, line_number, buffer, "duplicated entry in the interaction section")
          energy_parameters_init%int_energy_law_id(i1,i1s,i1a,i2,i2s,i2a) = i_law
          energy_parameters_init%int_energy_law_id(i2,i2s,i2a,i1,i1s,i1a) = i_law
          do i=1, n_shells
            if (read_num(words(7+i),energy_parameters_init%int_energy_pars(i1,i1s,i1a,i2,i2s,i2a,i)))&
              energy_parameters_init%int_energy_skip(i1,i1s,i1a,i2,i2s,i2a,i) = .false.
            energy_parameters_init%int_energy_pars(i2,i2s,i2a,i1,i1s,i1a,i) = &
                                energy_parameters_init%int_energy_pars(i1,i1s,i1a,i2,i2s,i2a,i)
            energy_parameters_init%int_energy_skip(i2,i2s,i2a,i1,i1s,i1a,i) = &
                                energy_parameters_init%int_energy_skip(i1,i1s,i1a,i2,i2s,i2a,i)
          end do
!            print*, 'int. law: ', energy_parameters_init%int_energy_law_id(i1,i2),&
!                    ' for species 1:', i1,&
!                    ' and species 2:', i2
!            print'(A,3e16.8,3L)', 'int. pars: ', energy_parameters_init%int_energy_pars(i1,i2,:)&
!                                              , energy_parameters_init%int_energy_skip(i1,i2,:)

        elseif (words(1) == '') then

          if (buffer == '') then
            parse_state = parse_state_default
!              print*, 'blank line '
!            else
!              print*, 'comment: ', trim(buffer)
          end if

        else

!            print*, 'unprocessed line: ', trim(buffer)
          call error_message(file_name, line_number, buffer, "unknown key")

        end if

    end do ! while ios=0

    energy_parameters_init%is_interaction = &
              any(.not.energy_parameters_init%int_energy_skip)
    ! Check the input consistency

  end function

end module energy_parameters_class
