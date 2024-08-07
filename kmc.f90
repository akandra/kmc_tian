program kmc_tian

use constants
use mc_lat_class
use control_parameters_class
use energy_parameters_class
use mmc
use kmc
use checks

implicit none

type(mc_lat) :: lat     ! Declare a variable of type mc_lat.
type(control_parameters) :: c_pars
type( energy_parameters) :: e_pars

character(len=max_string_length) file_name_base

! initialize vector of conditions for debug trap
debug = .false.
debug(10) = .false. ! check energy correction consistency in rc calculations

! Print the banner
print'(2A)', "kmc_tian. Release ",version

! Take a file name base
select case (command_argument_count())
    case(1)
      call get_command_argument(1,file_name_base)
    case default
      stop "Wrong number of arguments"
end select

! Read simulation parameters
c_pars = control_parameters_init(file_name_base)
e_pars = energy_parameters_init(c_pars)
!   initialize lattice
lat = mc_lat_init(c_pars, e_pars)

select case (c_pars%algorithm)

  case ('mmc')
    call metropolis(lat, c_pars, e_pars) ! Warning: mmc is implemented for the hops only

  case ('bkl')
    call Bortz_Kalos_Lebowitz(lat, c_pars, e_pars)

  case ('chk')
    call code_checks(lat, c_pars)

  case default
    stop 'Error: mc algorithm not defined'

end select

end program
