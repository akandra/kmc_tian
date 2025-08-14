module control_parameters_class

  use constants
  use open_file
  use utilities

  implicit none

  private
  public :: control_parameters_init

  type, public :: control_parameters

    character(len=10) :: algorithm  ! MC algorithm to use (bkl, mmc, mmc-gc)
    integer :: n_rows               ! number of rows
    integer :: n_cols               ! number of columns
    integer :: step_period          ! = 1/step_density (0 means no steps)
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
    integer  :: rdf_period            ! period for rdf_hist calculations
    real(dp) :: rdf_bin_size          ! size (in init cell) of the bin for rdf calculation
    integer  :: rdf_n_bins            ! number of bins for rdf_hist calculations
    logical  :: show_progress         ! if to show progress bar

    ! MMC-specific parameters

    integer :: n_mmc_steps          ! number of mmc steps
    integer :: hist_period          ! how often to calculate a histogram
    logical :: conf_save            ! key to save confs
    logical :: running_avgs_save    ! key to save running averages

    ! MMC-GC specific parameters

    real(dp), allocatable :: chem_pots(:) ! chemical potentials for the ideal gas (eV)
    integer :: gc_period                  ! period for gc-mmc (0 means no gc steps = canonical mmc)

    ! kMC-specific parameters

    integer  :: start_traj             ! starting kMC trajectory number
    integer  :: n_trajs                ! number of kMC trajectories
    integer, allocatable  :: n_bins(:) ! number of time intervals in kmc simulations
    real(dp) :: log_scale_t1           ! initial time for log binning
    real(dp) :: power_scale_t1         ! initial time for power binning
    real(dp), allocatable :: t_end(:)  ! kmc simulation times
    character(len=max_string_length) ::&
                   rate_file_name   ! name of the file with rate parameters
  contains

!    procedure :: read  => control_parameters_read

  end type

contains

  function control_parameters_init(file_name_base)

    type(control_parameters) control_parameters_init

    character(len=max_string_length), intent(in) :: file_name_base

    character(len=*), parameter :: err = "Error in the control file: "
    character(len=*), parameter :: warning = "Control file: "

    integer :: i, ios, nwords
    character(len=max_string_length) :: buffer
    character(len=max_string_length) :: words(100)

    real(dp),allocatable :: coverages(:)
    integer :: default_int = huge(0)

! TODO: produce proper definitions of defaults and strategy for testing the required keys
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
    control_parameters_init%rdf_bin_size     = -1.0_dp
    control_parameters_init%rdf_n_bins       = -1
    control_parameters_init%rdf_period       = -1
    control_parameters_init%show_progress    = .false.
    ! MMC-specific parameters
    control_parameters_init%gc_period        =  0
    control_parameters_init%n_mmc_steps      = -1
    control_parameters_init%hist_period      = -1
    control_parameters_init%conf_save        = .false.
    control_parameters_init%running_avgs_save= .false.
    ! kMC-specific parameters
    control_parameters_init%start_traj       =  1
    control_parameters_init%n_trajs          =  -1
    control_parameters_init%rate_file_name   = 'none'
    control_parameters_init%log_scale_t1     = -1.0_dp
    control_parameters_init%power_scale_t1   = -1.0_dp

    control_parameters_init%file_name_base = file_name_base
    !  read control parameters from the input file
    call open_for_read(inp_unit, trim(file_name_base)//'.control' )

    ios = 0
    do while (ios == 0)

      read(inp_unit, '(A)', iostat=ios) buffer
        ! ios < 0: end of record condition encountered or endfile condition detected
        ! ios > 0: an error is detected
        ! ios = 0  otherwise!

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
                if ( .not. read_int(words(2), control_parameters_init%n_rows) )&
                  stop err // "nlat parameter must be integer"
                control_parameters_init%n_cols = control_parameters_init%n_rows

              case (3)

                if ( .not. read_int(words(2), control_parameters_init%n_rows) .or. &
                     .not. read_int(words(3), control_parameters_init%n_cols)     )&
                  stop err // "nlat parameters must be integers"

              case default
                stop err // "nlat must have 1 or 2 parameters."

            end select

          case('step_period')

            if (nwords/=2) stop err // "step_period must have 1 parameter."
            if ( .not. read_int(words(2), control_parameters_init%step_period) )&
              stop err // "step_period must be integer"

          case('adsorbates')

            if (nwords==1) stop err // "adsorbates must have at least 1 parameter."
            allocate(control_parameters_init%ads_names(nwords - 1))
            do i=1,nwords-1
              read(words(i+1),*) control_parameters_init%ads_names(i)
            end do

          case('coverages')

            if (nwords==1) stop err // "coverages must have at least 1 parameter."
            allocate(coverages(nwords - 1))
            do i=1,nwords-1
              if (.not. read_num(words(i+1), coverages(i)) ) &
                stop err // "coverages must be numeric"
            end do

          case('temperature')

            if (nwords/=2) stop err // "temperature must have 1 parameter."
            if ( .not. read_num(words(2), control_parameters_init%temperature) )&
              stop err // "temperature must be a number"

            read(words(2),*) control_parameters_init%temperature

          case('energy')
            if (nwords/=2) stop err // "energy must have 1 parameter."
            control_parameters_init%energy_file_name = words(2)

          case('start_conf')
            if (nwords/=2) stop err // "start_conf must have 1 parameter."
            control_parameters_init%cfg_file_name = words(2)

          case('rdf')
            if (nwords/=4) stop err // "rdf must have 3 parameters."
            if ( .not. (read_int(words(2), control_parameters_init%rdf_period) .and. &
                        read_num(words(3), control_parameters_init%rdf_bin_size) .and. &
                        read_int(words(4), control_parameters_init%rdf_n_bins ) )) &
              stop err // "rdf's period, bin size, and number of bins must be integer, real, and integer, respectively."

          case('show_progress')
            if (nwords/=2) stop err // "show_progress must have 1 parameter."
            if (.not. read_logical(words(2), control_parameters_init%show_progress) ) &
              stop err // "show_progress must have type logical"

          case('gc_period')
            if (nwords/=2) stop err // "gc_period must have 1 parameter."
            if ( .not. read_int(words(2), control_parameters_init%gc_period) ) &
              stop err // " gc_period must be integer"

          case('gc_chempots')
            if (nwords==1) stop err // "gc_chempot must have at least 1 parameter."
            allocate(control_parameters_init%chem_pots(nwords - 1))
            do i=1,nwords-1
              if ( .not. read_num(words(i+1), control_parameters_init%chem_pots(i)) )&
                stop err // "gc_chempots must be numeric"
            end do

          case('mmc_save_period')
            if (nwords/=2) stop err // "save_period must have 1 parameter."
            if ( .not. read_int(words(2), control_parameters_init%save_period) ) &
              stop err // "mmc_save_period must be integer"

          case('mmc_nsteps')
            if (nwords/=2) stop err // "mmc_nsteps must have 1 parameter."
            if ( .not. read_int(words(2), control_parameters_init%n_mmc_steps) ) &
              stop err // "mmc_nsteps must be integer"

          case('mmc_hist_period')
            if (nwords/=2) stop err // "mmc_hist_period must have 1 parameter."
            if ( .not. read_int(words(2), control_parameters_init%hist_period) )&
              stop err // "mmc_hist_period must be integer"

          case('mmc_conf_save')
            if (nwords/=2) stop err // "mmc_conf_save must have 1 parameter."
            if ( .not. read_logical(words(2), control_parameters_init%conf_save) )&
              stop err // "mmc_conf_save must have type logical"

          case('mmc_running_avgs_save')
            if (nwords/=2) stop err // "mmc_running_avgs_save must have 1 parameter."
            if ( .not. read_logical(words(2), control_parameters_init%running_avgs_save) )&
              stop err // "mmc_running_avgs_save must have type logical"

          case('kmc_ntrajs')
            if ( nwords/=2 .and. nwords/=3 ) stop err // "kmc_ntrajs_period must have 1 parameter."
            if ( .not. read_int(words(2), control_parameters_init%n_trajs) )&
              stop err // "kmc_ntrajs must be integer."
            if (nwords==3) then
              if ( .not. read_int(words(3), control_parameters_init%start_traj) )&
                stop err // "kmc starting trajectory must be integer"
            endif

          case('kmc_time')
            if (nwords==1) stop err // "kmc_time must have at least 1 parameter."
            allocate(control_parameters_init%t_end(nwords - 1))
            do i=1,nwords-1
              if ( .not. read_num(words(i+1), control_parameters_init%t_end(i)) )&
                stop err // "kmc_time must be numeric"
            end do

          case('kmc_nbins')
            if (nwords==1) stop err // "kmc_nbins must have at least 1 parameter."
            allocate(control_parameters_init%n_bins(nwords - 1))
            do i=1,nwords-1
              if ( .not. read_int(words(i+1), control_parameters_init%n_bins(i)) )&
                stop err // "kmc_nbins must be integer"
            end do

          case('kmc_log_bin')
            if (nwords /= 2) stop err // "kmc_log_bin must have 1 parameter."
            if ( .not. read_num(words(2), control_parameters_init%log_scale_t1) )&
              stop err // "kmc_log_bin must be numeric"

          case('kmc_power_bin')
            if (nwords /= 2) stop err // "kmc_power_bin must have 1 parameter."
            if ( .not. read_num(words(2), control_parameters_init%power_scale_t1) )&
              stop err // "kmc_power_bin must be numeric"

          case('kmc_rates')
            if (nwords/=2) stop err // "kmc_rates must have 1 parameter."
            control_parameters_init%rate_file_name = words(2)

          case('')

          case default
            print *, warning // 'Skipping invalid key ', trim(words(1))
        end select

      end if

    end do ! ios

    close(inp_unit)

    ! Check the input consistency

    if (control_parameters_init%n_rows <= 0 .or. control_parameters_init%n_cols <= 0) then
      print*, " n_rows: ", control_parameters_init%n_rows
      print*, " n_cols: ", control_parameters_init%n_cols
      stop err // "lattice size parameters must be positive"
    end if

    if (size(coverages) /= size(control_parameters_init%ads_names))&
      stop err // "adsorbates and coverages are inconsistent"
    ! set the number of adsorbates
    control_parameters_init%n_species = size(coverages)

    do i=1,control_parameters_init%n_species
      if (coverages(i) < 0.0_dp ) stop err // " coverages have to be non-negative."
    end do

    if (control_parameters_init%step_period > 0 &
        .and. mod(control_parameters_init%n_cols,&
                  control_parameters_init%step_period) /= 0)&
      stop err // "inconsistent number of columns (nlat's second value) and step_period."

    if (control_parameters_init%temperature < 0) &
      stop err // " negative temperatures are not supported."

    if (control_parameters_init%rdf_period /= -1) then
      if (control_parameters_init%rdf_bin_size <= 0.0_dp .or.  &
          control_parameters_init%rdf_n_bins <= 0 ) &
        stop err // " rdf's bin size and number of bins must be positive."
    end if

    allocate(control_parameters_init%n_ads(control_parameters_init%n_species))
    ! Calculate number of adsorbate particles
    control_parameters_init%n_ads = nint(coverages*control_parameters_init%n_rows&
                                                 *control_parameters_init%n_cols)

    select case (control_parameters_init%algorithm)

    case('mmc')

      ! Set mmc_save_period to the proper value
      if (control_parameters_init%save_period == 0 )&
        control_parameters_init%save_period = control_parameters_init%n_mmc_steps + 1
      if (control_parameters_init%save_period == default_int )&
        stop err // "mmc_save_period is undefined."
      if (control_parameters_init%save_period < 0 )&
        stop err // "mmc_save_period must be non-negative."

      if (control_parameters_init%gc_period > 0) then
        ! TODO: implement gc for the mixture of species
        if (control_parameters_init%n_species > 1)&
            stop err // "grand canonical mmc is not implemented for more than 1 species. Consult the experts"
        ! check gc_coverages consistency
        if (size(coverages) /= size(control_parameters_init%chem_pots))&
            stop err // "gc_chempot and coverages are inconsistent"
      end if

    case('bkl')

      if (control_parameters_init%n_trajs <= 0)&
        stop err // "kmc_ntrajs must be positive or can be missing."
      if (control_parameters_init%start_traj <= 0)&
        stop err // "kmc starting trajectory must be positive."

      if ( .not.allocated(control_parameters_init%n_bins) ) stop err // "kmc_nbins is undefined."
      if ( .not.allocated(control_parameters_init%t_end)  ) stop err // "kmc_time is undefined."

      if ( size(control_parameters_init%n_bins) /=  size(control_parameters_init%t_end) ) &
        stop err // "kmc_bins and kmc_time are inconsistent."

      do i=1,size(control_parameters_init%n_bins)
        ! if (control_parameters_init%n_bins(i) == 0 ) control_parameters_init%n_bins(i) = 1
        if (control_parameters_init%n_bins(i) <= 0 ) stop err // "kmc_nbins must be positive."
      end do

      do i=1,size(control_parameters_init%t_end)
        if (control_parameters_init%t_end(i) <= 0 ) stop err // "kmc_time must be positive."
      end do

      if (control_parameters_init%rate_file_name == 'none')&
        stop err // "kmc_rates file name is not defined"

    end select

  end function

end module control_parameters_class
