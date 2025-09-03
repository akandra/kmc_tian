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

    ! Interaction energy law id (n_species x n_site_type x n_adsorption_sites x 
    !                            n_species x n_site_type x n_adsorption_sites)
    integer, dimension(:,:,:,:,:,:), allocatable :: int_energy_law_id
    integer, dimension(:,:,:,:,:,:), allocatable :: n_interactions

    ! Interaction energy (n_species x n_site_type x n_adsorption_sites x 
    !                     n_species x n_site_type x n_adsorption_sites x max_n_neighbors )
    real(dp), dimension(:,:,:,:,:,:,:),  allocatable :: int_energy_pars

    ! List of directions to neighbors 
    !                    (n_species x n_site_type x n_adsorption_sites x 
    !                     n_species x n_site_type x n_adsorption_sites x max_n_neighbors x 2)
    integer,  dimension(:,:,:,:,:,:,:,:),allocatable :: neighbor

    logical :: is_interaction
    ! Species importance key (n_species)
    logical, dimension(:), allocatable :: is_essential

    ! Default value for undefined energy
    real(dp) :: undefined_energy


  contains

!    procedure :: read  => control_parameters_read

  end type


contains

  function energy_parameters_init(control_pars)

    type(energy_parameters) energy_parameters_init

    type(control_parameters), intent(inout) :: control_pars

    integer :: i, ios, line_number, i1, i1s, i1a, i2, i2s, i2a, i_law, n

    character(len=max_string_length) :: buffer
    character(len=max_string_length), allocatable :: tokens(:), toks(:)
    integer :: ntokens, ntoks, nn
    character(len=len(trim(control_pars%energy_file_name))) :: file_name

    character(len=10) :: current_species_name
    integer           :: current_species_id

    integer           :: parse_state
    integer           :: parse_state_default      =0
    integer           :: parse_state_adsorption   =1
    integer           :: parse_state_interaction  =2

    logical :: first_warning = .true.
    integer,  parameter:: default_int = 0
    real(dp) :: temp_dp
    integer  :: temp_int

    i = control_pars%n_species
    allocate(energy_parameters_init%ads_energy(i,n_max_lat_site_types,n_max_ads_sites))
    allocate(energy_parameters_init%int_energy_law_id(i,n_max_lat_site_types,n_max_ads_sites,&
                                                      i,n_max_lat_site_types,n_max_ads_sites))
    allocate(energy_parameters_init%n_interactions(i,n_max_lat_site_types,n_max_ads_sites,&
                                                   i,n_max_lat_site_types,n_max_ads_sites))
    allocate(energy_parameters_init%is_essential(i))

    energy_parameters_init%undefined_energy = huge(0.0_dp)
    energy_parameters_init%ads_energy = energy_parameters_init%undefined_energy
    energy_parameters_init%int_energy_law_id = default_int
    energy_parameters_init%n_interactions    = 0
    !energy_parameters_init%int_energy_pars = energy_parameters_init%undefined_energy
    !energy_parameters_init%int_energy_skip = .true.
    energy_parameters_init%is_interaction  = .false.
    energy_parameters_init%is_essential    = .false.

    !  read energy definitions from the input file
    file_name = control_pars%energy_file_name

!-------------------------------------------------------------------------------
!   First pass of processing an .energy file:
!     * define the maximal number of interactions and do allocation for int_energy_pars
!     * check for parsing errors
!-------------------------------------------------------------------------------
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
      if (.not. get_tokens(buffer, ntokens, tokens)) then
        call error_message(file_name, line_number, buffer, "unpaired curly braces")
      else
        if (ntokens == 0) cycle ! skip empty lines

        if (tokens(1) =='adsorption') then

          if (parse_state /= parse_state_default) &
            call error_message(file_name, line_number, buffer, &
                      "invalid ending of the adsorption/interaction section")

          parse_state = parse_state_adsorption

          if (ntokens/=2 .and. ntokens/=3) &
            call error_message(file_name, line_number, buffer, &
                      "adsorption key must have 1 or 2 parameters")

        elseif (get_index(tokens(1),lat_site_names) /= 0) then

          if (parse_state /= parse_state_adsorption) then
            ! print *, 'parse_state: ', parse_state
            call error_message(file_name, line_number, buffer, &
                       "adsorption site definition outside of the adsorption section")
          end if

        elseif (tokens(1) == 'interaction') then

          if (parse_state /= parse_state_default) &
            call error_message(file_name, line_number, buffer, &
                      "invalid ending of the adsorption/interaction section")

          parse_state = parse_state_interaction

          if (ntokens/=1) call error_message(file_name, line_number, buffer, &
                            "interaction key must have no parameters")

        elseif (get_index(tokens(1),int_law_names) /= 0) then

          if (parse_state /= parse_state_interaction) &
            call error_message(file_name, line_number, buffer, &
                      "interaction law statement is outside of the interaction section")

        elseif (read_int(tokens(1), temp_int)) then

          if (parse_state /= parse_state_interaction) &
            call error_message(file_name, line_number, buffer, &
                      "line starts with an integer outside of the interaction section")

          energy_parameters_init%n_interactions(i1,i1s,i1a,i2,i2s,i2a) = &
               energy_parameters_init%n_interactions(i1,i1s,i1a,i2,i2s,i2a) + 1

        elseif (tokens(1) == section_end) then

          parse_state = parse_state_default

        else

!            print*, 'unprocessed line: ', trim(buffer)
          call error_message(file_name, line_number, buffer, "unknown key")

        end if
      end if ! get_tokens
    end do ! while ios=0

    close(inp_unit)

    allocate(energy_parameters_init%int_energy_pars(i,n_max_lat_site_types,n_max_ads_sites,&
                                                    i,n_max_lat_site_types,n_max_ads_sites,&
                                                  maxval(energy_parameters_init%n_interactions)))
    allocate(energy_parameters_init%neighbor(       i,n_max_lat_site_types,n_max_ads_sites,&
                                                    i,n_max_lat_site_types,n_max_ads_sites,&
                                                  maxval(energy_parameters_init%n_interactions),2))
    
    energy_parameters_init%int_energy_pars = 0.0_dp
    energy_parameters_init%neighbor = 0
    
!-------------------------------------------------------------------------------
!   Second pass of processing an .energy file:
!     * build the energy parameters structure
!-------------------------------------------------------------------------------

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
      if (.not. get_tokens(buffer, ntokens, tokens)) then
        call error_message(file_name, line_number, buffer, "unpaired curly braces")
      else
        if (ntokens == 0) cycle ! skip empty lines

        if (tokens(1) =='adsorption') then

          parse_state = parse_state_adsorption

          read(tokens(2),'(A)') current_species_name
          current_species_id = get_index(current_species_name, control_pars%ads_names )
          if (current_species_id == 0) &
            call error_message(file_name, line_number, buffer, &
                      "inconsistent adsorbate definition")
          if (ntokens > 2) then
            if (tokens(3) == essential_name) then
              energy_parameters_init%is_essential(current_species_id) = .true.
            else
              call error_message(file_name, line_number, buffer, &
                      "adsorption: unknown stopping trigger name")
            end if
          end if

!            print*, 'name     =', current_species_name
!            print*, 'id       =', current_species_id
!            print*, control_pars%ads_names

        elseif (get_index(tokens(1),lat_site_names) /= 0) then

          i1 = get_index(tokens(1),lat_site_names)
          i2 = get_index(tokens(2),ads_site_names)

          if ( i1==0 ) &
              call error_message(file_name, line_number, buffer, &
                          "unknown lattice site type in the adsorption section")
          if ( i2==0 ) &
              call error_message(file_name, line_number, buffer, &
                          "unknown adsorption site type in the adsorption section")
          if (energy_parameters_init%ads_energy(current_species_id,i1,i2 ) /= energy_parameters_init%undefined_energy)&
              call error_message(file_name, line_number, buffer, "duplicated entry")

          if (read_num(tokens(3), temp_dp)) then
            energy_parameters_init%ads_energy(current_species_id,i1,i2 ) = temp_dp
          else
            call error_message(file_name, line_number, buffer, &
                        "adsorption energy must be a number")
          end if

        elseif (tokens(1) == 'interaction') then

          parse_state = parse_state_interaction

        elseif (get_index(tokens(1),int_law_names) /= 0) then

          if (ntokens /= 7) then
            call error_message(file_name, line_number, buffer, &
                            "interaction law key must have 7 parameters")
          else

            i_law = get_index(tokens(1),int_law_names)
            i1    = get_index(tokens(2),control_pars%ads_names)
            i1s   = get_index(tokens(3),lat_site_names)
            i1a   = get_index(tokens(4),ads_site_names)
            i2    = get_index(tokens(5),control_pars%ads_names)
            i2s   = get_index(tokens(6),lat_site_names)
            i2a   = get_index(tokens(7),ads_site_names)

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

            if (energy_parameters_init%n_interactions(i1,i1s,i1a,i2,i2s,i2a) > 0) then
              if (energy_parameters_init%n_interactions(i2,i2s,i2a,i1,i1s,i1a) == 0) then
                energy_parameters_init%n_interactions(i2,i2s,i2a,i1,i1s,i1a) = &
                     energy_parameters_init%n_interactions(i1,i1s,i1a,i2,i2s,i2a)
              else
                call error_message(file_name, line_number, buffer, &
                            "duplicated entry in the interaction section",)
              end if

            end if

            do n=1,energy_parameters_init%n_interactions(i1,i1s,i1a,i2,i2s,i2a)

              read(inp_unit, '(A)', iostat=ios) buffer
              ! Split an input string
              if (get_tokens(buffer, ntokens, tokens)) then
                if (ntokens /= 3) &
                  call error_message(file_name, line_number, buffer, &
                              "neighbor interaction line must have 3 parameters")

                if (read_int(tokens(1),temp_int)) then
                  energy_parameters_init%neighbor(i1,i1s,i1a,i2,i2s,i2a,n,1) =  temp_int
                  ! Reverse direction for the opposite interaction
                  energy_parameters_init%neighbor(i2,i2s,i2a,i1,i1s,i1a,n,1) = -temp_int
                else
                  call error_message(file_name, line_number, buffer, &
                              "1st neighbor direction must be an integer")
                end if

                if (read_int(tokens(2),temp_int)) then
                  energy_parameters_init%neighbor(i1,i1s,i1a,i2,i2s,i2a,n,2) =  temp_int
                  ! Reverse direction for the opposite interaction
                  energy_parameters_init%neighbor(i2,i2s,i2a,i1,i1s,i1a,n,2) = -temp_int
                else
                  call error_message(file_name, line_number, buffer, &
                              "2nd neighbor direction must be an integer")
                end if

                if (read_num(tokens(3),temp_dp)) then
                  energy_parameters_init%int_energy_pars(i1,i1s,i1a,i2,i2s,i2a,n) = temp_dp
                  ! Symmetric interaction energy for the opposite interaction
                  energy_parameters_init%int_energy_pars(i2,i2s,i2a,i1,i1s,i1a,n) = temp_dp
                else
                  call error_message(file_name, line_number, buffer, &
                              "interaction energy must be a number")
                end if

              else
                call error_message(file_name, line_number, buffer, "unpaired curly braces")

              end if

              if ( all(energy_parameters_init%neighbor(i1,i1s,i1a,i2,i2s,i2a,n,:)==0) ) &
                error_message(file_name, line_number, buffer, &
                              "invalid neighbor direction (0,0)" )
            end do

          end if

        
        elseif (tokens(1) == section_end) then

          parse_state = parse_state_default

        else

!            print*, 'unprocessed line: ', trim(buffer)
          call error_message(file_name, line_number, buffer, "unknown key")

        end if
      end if ! get_tokens
    end do ! while ios=0


    close(inp_unit)

    energy_parameters_init%is_interaction = &
              any(.not.energy_parameters_init%int_energy_skip)

    ! Check the input consistency
    do i1 =1,control_pars%n_species
    do i1s=1,n_max_lat_site_types
    do i1a=1,n_max_ads_sites
      if (energy_parameters_init%ads_energy(i1,i1s,i1a) == energy_parameters_init%undefined_energy) then
        do i2 =1,control_pars%n_species
        do i2s=1,n_max_lat_site_types
        do i2a=1,n_max_ads_sites
          if (any(.not. energy_parameters_init%int_energy_skip(i1,i1s,i1a, i2,i2s,i2a, :))) then
            if (first_warning) then
              write(*,'(A)') ' Warning:'
              write(*,'(A)') '    file: '//trim(file_name)
              first_warning = .false.
            end if
            write(*,'(A)')   '    interaction section: adsorption energy for '//&
                                  trim(control_pars%ads_names(i1))//' '//&
                                  trim(lat_site_names(i1s))//' '//trim(ads_site_names(i1a))//' '//&
                                  'not defined'
          end if
        end do
        end do
        end do
      end if
    end do
    end do
    end do

  end function

end module energy_parameters_class
