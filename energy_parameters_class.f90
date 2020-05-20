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
    integer, dimension(:,:) :: int_energy_law_id
    ! Interaction energy (n_species x n_species x 3)
    real(dp), dimension(:,:,:) :: int_energy_pars

  contains

!    procedure :: read  => control_parameters_read

  end type


  interface energy_parameters

    module procedure :: energy_parameters_init

  end interface

contains

  function energy_parameters_init(control_pars)

    type(energy_parameters) energy_parameters_init

    type(control_parameters), intent(in) :: control_pars

    character(len=*), parameter :: err = "Error in the energy file: "
    character(len=*), parameter :: warning = "energy file: "

    integer :: i, ios, nwords, nwords1
    character(len=max_string_length) :: buffer
    character(len=max_string_length) :: words(100), words1(100)
    character(len=10) :: token
    logical :: found_it


    i = control_pars%n_species
    allocate(energy_parameters_init%ads_energy(i,n_site_types,n_ads_sites))
    allocate(energy_parameters_init%int_energy_law_id(i,i))
    allocate(energy_parameters_init%int_energy_pars(i,i,n_shells))

    energy_parameters_init%ads_energy = huge(0.0_dp)
    energy_parameters_init%int_energy_law_id = ''
    energy_parameters_init%int_energy_pars = 0.0_dp

    !  read energy definitions from the input file

    call open_for_read(inp_unit, trim(control_pars%energy_file_name) )

    ios = 0
    do while (ios == 0)

      read(inp_unit, '(A)', iostat=ios) buffer
        ! ios < 0: end of record condition encountered or endfile condition detected
        ! ios > 0: an error is detected
        ! ios = 0  otherwise

      if (ios == 0) then

        ! Split an input string
        words = ''
        call split_string(buffer, words, nwords)

        select case (words(1)) ! take a keyword

          case('adsorption')

            if (nwords/=2) stop err // "adsorption must have 1 parameter."
            read(words(2),'(A)') token
            found_it = .false.
            do i=1,control_pars%n_species
              if ( token == control_pars%ads_names(i) ) then
                found_it = .true.
                exit
              end if
            end do
            if (.not. found_it) stop err // "Inconsistent adsorbate definition."
            do
              read(inp_unit, '(A)') buffer
              if ( buffer == '' ) exit
              call split_string(buffer, words1, nwords1)

              select case (words1(1))

                case('terrace')

                  re



              end select

            end do

          case('nlat')

            select case (nwords)

              case (2)
                read(words(2),*) control_parameters_init%n_rows
                control_parameters_init%n_cols = control_parameters_init%n_rows

              case (3)
                read(words(2),*) control_parameters_init%n_rows
                read(words(3),*) control_parameters_init%n_cols

              case default
                stop err // "nlat must have 1 or 2 parameters."

            end select

          case('step_period')

            if (nwords/=2) stop err // "step_period must have 1 parameter."
            read(words(2),*) control_parameters_init%step_period

          case('adsorbates')

            if (nwords==1) stop err // "adsorbates must have at least 1 parameter."
            control_parameters_init%n_species = nwords - 1
            allocate(control_parameters_init%ads_names(nwords - 1))
            do i=1,nwords-1
              read(words(i+1),*) control_parameters_init%ads_names(i)
            end do

          case('coverages')

            if (nwords==1) stop err // "coverages must have at least 1 parameter."
            control_parameters_init%n_species = nwords - 1
            allocate(coverage(nwords - 1))
            do i=1,nwords-1
              read(words(i+1),*) coverage(i)
            end do

          case('temperature')

            if (nwords/=2) stop err // "temperature must have 1 parameter."
            read(words(2),*) control_parameters_init%temperature

          case('energy')
            if (nwords/=2) stop err // "energy must have 1 parameter."
            read(words(2),*) control_parameters_init%energy_file_name

          case('start_conf')
            if (nwords/=2) stop err // "start_conf must have 1 parameter."
            read(words(2),*) control_parameters_init%cfg_file_name

          case('save_period')
            if (nwords/=2) stop err // "save_period must have 1 parameter."
            read(words(2),*) control_parameters_init%save_period

          case('mmc_nsteps')
            if (nwords/=2) stop err // "mmc_nsteps must have 1 parameter."
            read(words(2),*) control_parameters_init%n_mmc_steps

          case('mmc_hist_period')
            if (nwords/=2) stop err // "mmc_hist_period must have 1 parameter."
            read(words(2),*) control_parameters_init%hist_period

          case('kmc_ntrajs')
            if (nwords/=2) stop err // "kmc_ntrajs_period must have 1 parameter."
            read(words(2),*) control_parameters_init%n_trajs

          case('kmc_time')
            if (nwords/=2) stop err // "kmc_time must have 1 parameter."
            read(words(2),*) control_parameters_init%t_end

          case('kmc_nbins')
            if (nwords/=2) stop err // "kmc_nbins must have 1 parameter."
            read(words(2),*) control_parameters_init%n_bins

          case('kmc_rates')
            if (nwords/=2) stop err // "kmc_rates must have 1 parameter."
            read(words(2),*) control_parameters_init%rate_file_name

          case('')

          case default
            print *, warning // 'Skipping invalid key ', trim(words(1))
        end select

      end if

    end do ! ios

    close(inp_unit)

    ! Check the input consistency
    if (size(coverage) /= size(control_parameters_init%ads_names))&
      stop err // "adsorbates and coverages are inconsistent"

    allocate(control_parameters_init%n_ads(control_parameters_init%n_species))
    ! Calculate number of adsorbate particles
    control_parameters_init%n_ads = nint(coverage*control_parameters_init%n_rows&
                                                 *control_parameters_init%n_cols)

  end function

end module energy_parameters_class
