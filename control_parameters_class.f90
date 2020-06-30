module control_parameters_class

  use constants
  use open_file
  use utilities

  implicit none

  private
  public :: control_parameters_init

  type, public :: control_parameters

    character(len=3) :: algorithm   ! MC algorithm to use (bkl or mmc)
    integer :: n_rows               ! number of rows
    integer :: n_cols               ! number of columns
    integer :: step_period          ! = step_density^-1, 0 means no steps
    integer :: n_species            ! number of the adsorbate types
    integer :: save_period          ! period for the output
    character(len=10),&
        allocatable :: ads_names(:) ! adsorbate names
    integer, allocatable :: n_ads(:)! initial number of adsorbates
    real(dp):: temperature          ! in K
    character(len=max_string_length) ::&
                   energy_file_name,& ! name of the file with adsorption and interaction energies
                   cfg_file_name,   & ! name of the file with initial configuration
                   file_name_base     ! filename base for output files

    ! MMC-specific parameters

    integer :: n_mmc_steps          ! number of mmc steps
    integer :: hist_period          ! period for histogram  calculation

    ! kMC-specific parameters

    integer  :: n_trajs             ! number of kMC trajectories
    integer  :: n_bins              ! number of time intervals in kmc histogram
    real(dp) :: t_end               ! kmc simulation time
    character(len=max_string_length) ::&
                   rate_file_name   ! name of the file with rate parameters
  contains

!    procedure :: read  => control_parameters_read

  end type


  interface control_parameters

    ! dja change to make f2003 complient
    !module procedure :: control_parameters_init
    module procedure control_parameters_init

  end interface

contains

  function control_parameters_init(file_name_base)

    type(control_parameters) control_parameters_init

    character(len=max_string_length), intent(in) :: file_name_base

    character(len=*), parameter :: err = "Error in the control file: "
    character(len=*), parameter :: warning = "Control file: "

    integer :: i, ios, nwords
    character(len=max_string_length) :: buffer
    character(len=max_string_length) :: words(100)

    real(dp),allocatable :: coverage(:)
    integer :: default_int = huge(0)




    control_parameters_init%algorithm        = ''
    control_parameters_init%n_rows           = -1
    control_parameters_init%n_cols           = -1
    control_parameters_init%step_period      = -1
    control_parameters_init%n_species        = -1
    control_parameters_init%save_period      = default_int
    control_parameters_init%temperature      = -1.0_dp
    control_parameters_init%energy_file_name = 'none'
    control_parameters_init%cfg_file_name    = 'none'
    control_parameters_init%file_name_base   = 'none'
    ! MMC-specific parameters
    control_parameters_init%n_mmc_steps      = -1
    control_parameters_init%hist_period      = -1
    ! kMC-specific parameters
    control_parameters_init%n_trajs          = -1
    control_parameters_init%n_bins           = default_int
    control_parameters_init%t_end            = -1.0_dp
    control_parameters_init%rate_file_name   = 'none'

    control_parameters_init%file_name_base = file_name_base
    !  read control parameters from the input file
    call open_for_read(inp_unit, trim(file_name_base)//'.control' )

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

          case('algorithm')

            if (nwords/=2) stop err // "algorithm must have 1 parameter."
            read(words(2),'(A)') control_parameters_init%algorithm
            !print*, 'algorithm is .',control_parameters_init%algorithm,'.'

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

          case('mmc_save_period')
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

    if (control_parameters_init%step_period > 0 &
        .and. mod(control_parameters_init%n_cols,&
                  control_parameters_init%step_period) /= 0)&
      stop err // "inconsistent number of columns (nlat's second value) and step_period."

    allocate(control_parameters_init%n_ads(control_parameters_init%n_species))
    ! Calculate number of adsorbate particles
    control_parameters_init%n_ads = nint(coverage*control_parameters_init%n_rows&
                                                 *control_parameters_init%n_cols)

    select case (control_parameters_init%algorithm)

    case('mmc')
      ! Set mmc_save_period to the proper value
      if (control_parameters_init%save_period == 0 )&
        control_parameters_init%save_period = control_parameters_init%n_mmc_steps + 1
      if (control_parameters_init%save_period == default_int )&
        stop err // "mmc_save_period is undefined."
      if (control_parameters_init%save_period < 0 )&
        stop err // "mmc_save_period is negative."

    case('bkl')
      ! Check kmc_nbins
      if (control_parameters_init%n_bins == 0 ) control_parameters_init%n_bins = 1
      if (control_parameters_init%n_bins == default_int )&
        stop err // "kmc_nbins is undefined."
      if (control_parameters_init%n_bins < 0 )&
        stop err // "kmc_nbins is negative."


    end select

  end function

end module control_parameters_class
