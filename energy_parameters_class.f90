module energy_parameters_class

  use constants
  use control_parameters_class
  use open_file
  use utilities

  implicit none

  type, public :: energy_parameters

    ! Adsorpton energies (n_species x n_site_type x n_adsorption_sites)
    real(dp), dimension(:,:,:), allocatable :: ads_energy
    ! Interaction energy law id (n_species x n_species)
    integer, dimension(:,:), allocatable :: int_energy_law_id
    ! Interaction energy (n_species x n_species x 3)
    real(dp), dimension(:,:,:),allocatable :: int_energy_pars

  contains

!    procedure :: read  => control_parameters_read

  end type


  interface energy_parameters

    module procedure :: energy_parameters_init

  end interface

contains

  function energy_parameters_init(control_pars)

    type(energy_parameters) energy_parameters_init

    type(control_parameters), intent(inout) :: control_pars

    character(len=*), parameter :: err = "Error in the energy file: "
    character(len=*), parameter :: warning = "energy file: "

    integer :: i, ios, nwords, nwords1
    character(len=max_string_length) :: buffer
    character(len=max_string_length) :: words(100), words1(100)

    character(len=10) :: current_species_name
    integer           :: current_species_id

    integer           :: parse_state
    integer           :: parse_state_default      =0
    integer           :: parse_state_adsorption   =1
    integer           :: parse_state_interaction  =2

    logical :: found_it


    i = control_pars%n_species
    allocate(energy_parameters_init%ads_energy(i,n_site_types,n_ads_sites))
    allocate(energy_parameters_init%int_energy_law_id(i,i))
    allocate(energy_parameters_init%int_energy_pars(i,i,n_shells))

    energy_parameters_init%ads_energy = huge(0.0_dp)
!    energy_parameters_init%int_energy_law_id = ''
    energy_parameters_init%int_energy_law_id = 0
    energy_parameters_init%int_energy_pars = 0.0_dp

    !  read energy definitions from the input file

    call open_for_read(inp_unit, trim(control_pars%energy_file_name) )

    ios = 0
    parse_state = parse_state_default

    do while (ios == 0)

      read(inp_unit, '(A)', iostat=ios) buffer
        ! ios < 0: end of record condition encountered or endfile condition detected
        ! ios > 0: an error is detected
        ! ios = 0  otherwise

      if (ios /= 0) exit

        ! Split an input string
        words = ''
        call split_string(buffer, words, nwords)

        select case (words(1)) ! take a keyword

          case('adsorption')

            parse_state = parse_state_adsorption
            if (nwords/=2) stop err // "adsorption must have 1 parameter."

            read(words(2),'(A)') current_species_name
            current_species_id = get_index(current_species_name, control_pars%ads_names )
            if (current_species_id == 0) stop err // "Inconsistent adsorbate definition."

            print*, 'name     =', current_species_name
            print*, 'id       =', current_species_id
            print*, control_pars%ads_names

          case default
            print*, 'unprocessed line: ', trim(buffer)

        end select

    end do ! while ios=0



stop 123


    ! Check the input consistency

  end function

end module energy_parameters_class
