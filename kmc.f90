program kmc_tian

use mc_lat_class
use control_parameters_class
use energy_parameters_class
use rates_class
use mmc
use kmc
use utilities

implicit none

type(mc_lat) :: lat     ! Declare a variable of type mc_lat.
type(control_parameters) :: c_pars
type( energy_parameters) :: e_pars

character(len=max_string_length) file_name_base

! Take a file name base
select case (command_argument_count())
    case(1)
      call get_command_argument(1,file_name_base)
    case default
      stop "Wrong number of arguments"
end select

! Read simulation parameters
c_pars = control_parameters_init(file_name_base)
!print*, control_pars
e_pars = energy_parameters_init(c_pars)
!print*, energy_pars
!   initialize lattice
lat = mc_lat_init(c_pars, e_pars)

!  call lat%print_st
!  call lat%print_ocs
!  call lat%print_ads

select case (c_pars%algorithm)

  case ('mmc')
    call metropolis(lat, c_pars, e_pars)

  case ('bkl')
    call Bortz_Kalos_Lebowitz(lat, c_pars, e_pars)

  case default
    stop 'Error: mc algorithm not defined'

end select

call open_for_write(outcfg_unit,trim(c_pars%file_name_base)//'.out')

write(outcfg_unit,'(A10,A10,A15)') "# rows","# cols","step_period"
write(outcfg_unit,'(3i10)') lat%n_rows, lat%n_cols, c_pars%step_period
write(outcfg_unit,'(100A10)') adjustr(c_pars%ads_names)
write(outcfg_unit,'(100i10)') lat%n_ads
write(outcfg_unit,'(5A10)') "#","row","col","ads_site", "species"
call lat%print_ads(outcfg_unit)


close(outcfg_unit)

end program
