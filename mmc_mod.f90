module mmc

  use mc_lat_class
  use energy_parameters_class
  use energy_mod

  implicit none

contains

subroutine metropolis(lat, c_pars, e_pars)

  type(mc_lat), intent(in) :: lat
  type(control_parameters), intent(in) :: c_pars
  type(energy_parameters ), intent(in) :: e_pars


  integer :: i, istep, n_ads_total
  real(dp) :: energy_old

!
!        call open_for_write(outeng_unit,trim(fname)//'.en')
!        write(outeng_unit,*) 0,&
!            total_energy(nlat, nads, nnn, occupations, site_type, &
!                         ads_list, nn_list, ads_energy, int_energy)/eV2K
!
!!        write(*,cfg_fmt) transpose(site_type)
!!        print*
!!        write(*,cfg_fmt) transpose(occupations)
!!        print*
!!        print*,energy_file
!!        write(*,'(3f12.3)') int_energy/eV2K
!!        write(*,*) total_energy(nlat, nads, nnn, occupations, site_type, &
!!                         ads_list, nn_list, ads_energy, int_energy)/eV2K
!!        stop 33
!
  !loop over mmc steps
  do istep=1, c_pars%n_mmc_steps

    n_ads_total = 0
    do i=1,c_pars%n_species
      n_ads_total = n_ads_total + lat%n_ads(i)
    end do

    do i=1, n_ads_total

      energy_old = energy(i, lat, e_pars)
      print *, energy_old
      pause
!
!                ihop = floor(nnn*ran1()) + 1
!
!                i_new = modulo(ads_list(i,1) + nn_list(ihop,1)-1,nlat) + 1
!                j_new = modulo(ads_list(i,2) + nn_list(ihop,2)-1,nlat) + 1
!
!                if (occupations(i_new, j_new) == 0) then
!
!                    i_old = ads_list(i,1)
!                    j_old = ads_list(i,2)
!
!                    occupations(i_new,j_new) = i
!                    occupations(i_old,j_old) = 0
!                    ads_list(i,:) = (/i_new,j_new/)
!
!                    delta_E = energy(i, nlat, nads, nnn, occupations, site_type, &
!                                     ads_list, nn_list, ads_energy, int_energy) &
!                            - energy_old
!
!                    if (exp(- delta_E/temperature) < ran1()) then
!
!                        occupations(i_old,j_old) = i
!                        occupations(i_new,j_new) = 0
!                        ads_list(i,:) = (/i_old,j_old/)
!
!                    end if
!
!                end if
!
    enddo
!
!!print*,istep, hist_period, mod(istep, hist_period)
!
!        if (hist_period > 0 .and. mod(istep, hist_period) == 0) then
!            call hoshen_kopelman(cluster_label, largest_label, occupations, ads_list, &
!                                                nn_list, nads, nlat, nnn)
!
!            call count_cluster_sizes(cluster_sizes, cluster_label,&
!                                                                ads_list, nads, nlat)
!
!            hist_counter = hist_counter + 1
!            do i=1,largest_label
!                if (cluster_sizes(i) > 0) &
!                    hist(cluster_sizes(i)) = hist(cluster_sizes(i)) + 1
!            end do
!
!        !print*,largest_label
!        !print*,cluster_sizes
!        !print*
!        !write(*,cfg_fmt) (cluster_label(m,:), m=1,nlat)
!        !print*
!        !write(*,cfg_fmt) hist
!        !
!
!        end if
!
!            if (mod(istep, save_period) == 0) then
!                print*, istep
!                write(outcfg_unit,cfg_fmt) transpose(occupations)
!                write(outeng_unit,*) istep-1,&
!                    total_energy(nlat, nads, nnn, occupations, site_type, &
!                                 ads_list, nn_list, ads_energy, int_energy)/eV2K
!
!                if (hist_period > 0) then
!                    write(outhst_unit,*) hist_counter
!                    write(outhst_unit,ads_fmt) hist
!                end if
!           end if
!
    enddo
!
!        close(outeng_unit)
!


end subroutine metropolis


end module mmc
