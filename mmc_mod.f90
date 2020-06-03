module mmc

  use constants
  use utilities
  use mc_lat_class
  use control_parameters_class
  use energy_parameters_class
  use energy_mod

  implicit none

  private
  public::  metropolis

contains

subroutine metropolis(lat, c_pars, e_pars)

  type(mc_lat), intent(inout) :: lat
  type(control_parameters), intent(in) :: c_pars
  type(energy_parameters ), intent(in) :: e_pars

  integer :: i, istep, ihop, species
  integer :: new_row, new_col, old_row, old_col
  real(dp) :: energy_old, beta, delta_E
  character(len=max_string_length) :: n_ads_fmt
  integer, dimension(lat%n_rows,lat%n_cols) :: cluster_label
  integer, dimension(maxval(lat%n_ads)) :: cluster_sizes
  integer :: largest_label, hist_counter
  integer, dimension(c_pars%n_species,maxval(lat%n_ads)) :: hist

  ! inverse thermodynamic temperature
  beta = 1.0_dp/(kB*c_pars%temperature)

  ! Output formats
  write(n_ads_fmt,'(i10)') lat%n_ads_tot()
  n_ads_fmt = '('//trim(adjustl(n_ads_fmt))//'i8)'

  ! write initial state of the lattice to a file
  call open_for_write(outcfg_unit,trim(c_pars%file_name_base)//'.confs')

  write(outcfg_unit,'(A10,A10,A15)') &
                "# rows","# cols","step_period"
  write(outcfg_unit,'(3i10)') lat%n_rows, lat%n_cols, c_pars%step_period
  write(outcfg_unit,'(100A10)') adjustr(c_pars%ads_names)
  write(outcfg_unit,'(A10,i0)') "mmc step ", 0
  write(outcfg_unit,'(100i10)') lat%n_ads
!  write(outcfg_unit,'(5A10)') "#","row","col","ads_site", "species"
  call lat%print_ads(outcfg_unit)

  ! write initial total energy of the system
  call open_for_write(outeng_unit,trim(c_pars%file_name_base)//'.en')
  write(outeng_unit,*) 0, total_energy(lat,e_pars)

  ! open file for saving cluster size histogram
  if (c_pars%hist_period > 0) &
    call open_for_write(outhst_unit,trim(c_pars%file_name_base)//'.hist')
  ! Initialize cluster histogram counters
  hist_counter = 0
  hist = 0

  !loop over mmc steps
  do istep=1, c_pars%n_mmc_steps

    do i=1, lat%n_ads_tot()

      energy_old = energy(i, lat, e_pars)

      ! We consider hops to the nearest-neighbor cells only
      ! And we doubt that we ever need something else
      ihop = floor(lat%n_nn(1)*ran1()) + 1

      call lat%hop(i,ihop,new_row,new_col)

      if (lat%occupations(new_row,new_col) == 0) then

        old_row = lat%ads_list(i)%row
        old_col = lat%ads_list(i)%col

        lat%occupations(new_row,new_col) = i
        lat%occupations(old_row,old_col) = 0
        lat%ads_list(i)%row = new_row
        lat%ads_list(i)%col = new_col

        delta_E = energy(i, lat, e_pars) - energy_old

        if (exp(- beta*delta_E) < ran1()) then
          ! reject the hop
          lat%occupations(new_row,new_col) = 0
          lat%occupations(old_row,old_col) = i
          lat%ads_list(i)%row = old_row
          lat%ads_list(i)%col = old_col
        end if

      end if

    enddo

    if (c_pars%hist_period > 0 .and. mod(istep, c_pars%hist_period) == 0) then
      hist_counter = hist_counter + 1
      do species=1,c_pars%n_species
        call lat%hoshen_kopelman(species, cluster_label, largest_label)
        call lat%cluster_size(species, cluster_label, cluster_sizes)
        do i=1,largest_label
          if (cluster_sizes(i) > 0) &
            hist(species,cluster_sizes(i)) = hist(species,cluster_sizes(i)) + 1
        end do
!        print*,'species =',species, 'largest_label = ',largest_label
!        print*,'cluster_sizes = ',cluster_sizes
!        print*,'hist:',hist(species,:)
!        print*
      end do
    end if

    if (mod(istep, c_pars%save_period) == 0) then
      write(*,'(A,i0)'), 'mmc step = ', istep
      ! Save configuration
      write(outcfg_unit,'(A10,i0)') "mmc step ",istep
      write(outcfg_unit,'(100i10)') lat%n_ads
      call lat%print_ads(outcfg_unit)
      ! Save energy
      write(outeng_unit,*) istep, total_energy(lat,e_pars)
      ! Save cluster size histogram
      if (c_pars%hist_period > 0) then
        write(outhst_unit,*) 'counts ', hist_counter
        do species=1,c_pars%n_species
          write(outhst_unit,*) species
          write(outhst_unit,n_ads_fmt) hist(species,:)
        end do
! Warning: needs discussion
!        hist_counter = 0
!        hist = 0
      end if
    end if

  enddo ! over mmc steps

  close(outcfg_unit)
  close(outeng_unit)
  if (c_pars%hist_period > 0) close(outhst_unit)

end subroutine metropolis


end module mmc
