module control_parameters_class

  use constants
  use open_file
  use utilities

  implicit none

  type, public :: control_parameters

    character(len=3) :: algorithm   ! MC algorithm to use (kmc or mmc)
    integer :: n_rows               ! number of rows
    integer :: n_cols               ! number of columns
    integer :: step_period      	! = step_density^-1, 0 means no steps
    integer :: n_species            ! number of the adsorbate types
    integer :: save_period          ! period for the output
    character(len=10),&
        allocatable :: ads_names(:) ! adsorbate names
    integer, allocatable :: n_ads(:)! initial number of adsorbates
    real(dp):: temperature          ! in K
    character(len=max_string_length) ::&
                   energy_file_name,& ! name of the file with adsorption and interaction energies
                   cfg_file_name      ! name of the file with initial configuration

    ! MMC-specific parameters

    integer :: n_mmc_steps          ! number of mmc steps
    integer :: hist_period          ! period for histogram  calculation

    ! kMC-specific parameters

    integer  :: n_trajs	            ! number of kMC trajectories
    integer  :: n_bins              ! number of time intervals in kmc histogram
    real(dp) :: t_end               ! kmc simulation time
    character(len=max_string_length) ::&
                   rate_file_name   ! name of the file with rate parameters
  contains

!    procedure :: read  => control_parameters_read

  end type


  interface control_parameters

    module procedure :: control_parameters_init

  end interface

contains

  function control_parameters_init(file_name_base)

    type(control_parameters) control_parameters_init

    character(len=max_string_length), intent(in) :: file_name_base

    character(len=*), parameter :: err = "Error in the control file: "

    integer :: ios, nwords
    character(len=max_string_length) :: buffer
    character(len=max_string_length) :: words(100)




    control_parameters_init%algorithm        = ''
    control_parameters_init%n_rows           = -1
    control_parameters_init%n_cols           = -1
    control_parameters_init%step_period      = -1
    control_parameters_init%n_species        = -1
    control_parameters_init%save_period      = -1
    control_parameters_init%temperature      = -1.0_dp
    control_parameters_init%energy_file_name = 'none'
    control_parameters_init%cfg_file_name    = 'none'
    ! MMC-specific parameters
    control_parameters_init%n_mmc_steps      = -1
    control_parameters_init%hist_period      = -1
    ! kMC-specific parameters
    control_parameters_init%n_trajs          = -1
    control_parameters_init%n_bins           = -1
    control_parameters_init%t_end            = -1.0_dp
    control_parameters_init%rate_file_name   = 'none'

    !  read control parameters from the input file

    call open_for_read(inp_unit, trim(file_name_base)//'.in' )

    ios = 0
    do while (ios == 0)

      read(inp_unit, '(A)', iostat=ios) buffer
        ! ios < 0: end of record condition encountered or endfile condition detected
        ! ios > 0: an error is detected
        ! ios = 0  otherwise

      if (ios == 0) then

        ! Split an input string
        call split_string(buffer, words, nwords)

        select case (words(1)) ! take a keyword

          case('algorithm')

            if (control_parameters_init%algorithm /= '')&
              stop err // "Multiple use of the algorithm key"

            read(words(2),'(A)') control_parameters_init%algorithm

            print*, 'algorithm is .',control_parameters_init%algorithm,'.'

!          case('nlat')
!              read(buffer,*,iostat=ios) nlat
!          case('step_period')
!              read(buffer,*,iostat=ios) step_period
!          case('temperature')
!              read(buffer,*,iostat=ios) temperature
!              beta = 1.0d0/temperature ! thermodynamic temperature
!          case('coverage')
!              read(buffer,*,iostat=ios) coverage
!          case('energy')
!              read(buffer,*,iostat=ios) energy_file
!              energy_file = trim(energy_file)
!          case('save_period')
!              read(buffer,*,iostat=ios) save_period
!          case('ini_conf')
!              read(buffer,*,iostat=ios) cfg_fname
!          case('mmc_hist_period')
!              read(buffer,*,iostat=ios) hist_period
!          case('mmc_nsteps')
!              read(buffer,*,iostat=ios) nsteps
!          case('kmc_ntrajs')
!              read(buffer,*,iostat=ios) ntrajs
!          case('kmc_time')
!              read(buffer,*,iostat=ios) t_end
!          case('kmc_nbins')
!              read(buffer,*,iostat=ios) n_bins
!          case('kmc_rates')
!              read(buffer,*,iostat=ios) rate_file
!              rate_file = trim(rate_file)
!          case default
!              if (label(1:1) /= '!')&
!                  print *, 'Skipping invalid label at line', label
        end select

      end if

    end do ! ios

    close(inp_unit)


  end function

end module control_parameters_class
