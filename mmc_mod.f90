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

  integer :: i, istep, ihop, species, species1, species2, j
  integer :: new_row, new_col, new_ads_site, old_row, old_col, old_ads_site, old_lst
  real(dp) :: energy_old, beta, delta_E
  character(len=max_string_length) :: n_ads_fmt, rdf_fmt, version_header
  integer, dimension(lat%n_rows,lat%n_cols) :: cluster_label
  integer, dimension(lat%n_rows*lat%n_cols) :: cluster_sizes
  integer :: largest_label, hist_counter, rdf_counter
  integer, dimension(c_pars%n_species,lat%n_rows*lat%n_cols) :: hist
  integer, dimension(c_pars%n_species,c_pars%n_species,c_pars%rdf_n_bins) :: rdf_hist
  integer :: n_sites, row, col, new_n_ads, n_ads_sites, i_ads, n_ads_tot_old, i_rand, i_lst
  real(dp) :: delta, probability, beta_delta, cov_factor
  logical :: remove
  ! Running averages
  ! lst-specific coverage in ML
  integer,  dimension(c_pars%n_species,n_max_lat_site_types) :: lst_pops
  ! replica-specific coverage in ML
  integer,  dimension(c_pars%n_species,lat%n_cols/c_pars%step_period) :: replica_pops
  ! vector of zigzag-related coverage in ML
  integer :: ns(4), zigzag(2,2,2,2)

!  real(dp), dimension(c_pars%rdf_n_bins) :: dr2
!  real(dp), dimension(c_pars%n_species) :: coverage
!  real(dp) :: unit_cell_area

  n_sites = lat%n_rows*lat%n_cols
  cov_factor = 1.0_dp/(n_sites*c_pars%save_period)

  ! inverse thermodynamic temperature
  beta = 1.0_dp/(kB*c_pars%temperature)
  ! Factors to calculate rdf from adsorbate counts
!  unit_cell_area = abs(lat%lat_vec_1(1)*lat%lat_vec_2(2) - lat%lat_vec_1(2)*lat%lat_vec_2(1))
!  coverage = c_pars%n_ads/float(c_pars%n_rows*c_pars%n_cols)
!  dr2(1) = 1.0_dp
!  do i=2,c_pars%rdf_n_bins
!    dr2(i) = unit_cell_area / (2*pi*(i - 1)*c_pars%rdf_bin_size**2)
!  end do


  ! Output formats
  write(n_ads_fmt,'(i10)') lat%n_ads_tot()
  n_ads_fmt = '('//trim(adjustl(n_ads_fmt))//'i8)'
  write(rdf_fmt,'(i10)') c_pars%rdf_n_bins
  rdf_fmt = '('//trim(adjustl(rdf_fmt))//'i8)'
  version_header = '! kmc_tian Release ' // version

  ! write initial state of the lattice to a file
  call open_for_write(outcfg_unit,trim(c_pars%file_name_base)//'.confs')

  write(outcfg_unit,'(A)') trim(version_header)
  write(outcfg_unit,'(A10,A10,A15)') &
                "# rows","# cols","step_period"
  write(outcfg_unit,'(3i10)') lat%n_rows, lat%n_cols, c_pars%step_period
  write(outcfg_unit,'(100A10)') adjustr(c_pars%ads_names)
  write(outcfg_unit,'(A10,i0)') "mmc_step ", 0
  write(outcfg_unit,'(100i10)') lat%n_ads
!  write(outcfg_unit,'(5A10)') "#","row","col","ads_site", "species"
  call lat%print_ads(outcfg_unit)

  ! write initial total energy of the system
  call open_for_write(outeng_unit,trim(c_pars%file_name_base)//'.en')
  write(outeng_unit,'(3A12)') 'mmc_step', ' energy(eV)', 'n_ads'
  write(outeng_unit,*) 0, total_energy(lat,e_pars), lat%n_ads

  ! open file for saving running averages
  if (c_pars%running_avgs_save > 0) call open_for_write(outrav_unit,trim(c_pars%file_name_base)//'.rav')
  ! Initialize accumulators
  lst_pops = 0
  replica_pops = 0
  zigzag = 0

  ! open file for saving cluster size histogram
  if (c_pars%hist_period > 0) call open_for_write(outhst_unit,trim(c_pars%file_name_base)//'.hist')
  ! Initialize cluster histogram and counters
  hist_counter = 0
  hist = 0

  ! open file for saving rdf
  if (c_pars%rdf_period > 0) call open_for_write(outrdf_unit,trim(c_pars%file_name_base)//'.rdf')
  ! rdf histogram and counters
  rdf_counter = 0
  rdf_hist = 0

  if (c_pars%show_progress > 0) then
    write(*,'(20X,A)') "M.M.C. Code's progress report:"
    call progress_bar(0)
  end if

  !loop over mmc steps
  do istep=1, c_pars%n_mmc_steps

    ! Loop over all adsorbates to do a canonical mmc step
    do i=1, lat%n_ads_tot()

      energy_old = energy(i, lat, e_pars)

      old_row = lat%ads_list(i)%row
      old_col = lat%ads_list(i)%col
      old_ads_site = lat%ads_list(i)%ast
      old_lst = lat%lst(old_row,old_col)

      ! We consider hops to the nearest-neighbor cells only
      ! And we doubt that we ever need something else
      ihop = floor(lat%n_nn(old_lst,1)*ran1()) + 1

      call lat%hop(i,ihop,new_row,new_col,new_ads_site)

!      print*, 'ads. no. = ',i,'row = ', new_row,'col = ', new_col,&
!              'ads. site = ', ads_site_names(new_ads_site)
!      pause

      if (lat%occupations(new_row,new_col) == 0) then

        lat%occupations(new_row,new_col) = i
        lat%occupations(old_row,old_col) = 0
        lat%ads_list(i)%row = new_row
        lat%ads_list(i)%col = new_col
        lat%ads_list(i)%ast = new_ads_site

        delta_E = energy(i, lat, e_pars) - energy_old

        ! Accept if delta_E is non-positive
        ! Calculate probability if delta_E is positive
        if ( delta_E > 0.0_dp ) then
          probability = exp(- beta*delta_E)
          if ( probability < ran1() ) then
            ! reject the hop
            lat%occupations(new_row,new_col) = 0
            lat%occupations(old_row,old_col) = i
            lat%ads_list(i)%row = old_row
            lat%ads_list(i)%col = old_col
            lat%ads_list(i)%ast = old_ads_site
          end if
        end if

      end if

    end do

    ! Loop over  grand-canonical steps
    ! WARNING: implemented for 1 species only
    species = 1
    if (c_pars%gc_period > 0 .and. mod(istep, c_pars%gc_period) == 0) then

      ! select a random site
      i_rand = irand(n_sites)
      row = (i_rand-1)/lat%n_cols + 1
      col = i_rand - (row - 1)*lat%n_cols

      if (lat%occupations(row,col) == 0) then

        ! Do trial addition of an adsorbate

        ! increase the number of adsorbates
        lat%n_ads(species) = lat%n_ads(species)  + 1
        ! update occupations
        new_n_ads = lat%n_ads_tot()
        lat%occupations(row,col) = new_n_ads
        ! update adsorbate list
        lat%ads_list(new_n_ads)%id  = species
        lat%ads_list(new_n_ads)%row = row
        lat%ads_list(new_n_ads)%col = col
        n_ads_sites = size( lat%avail_ads_sites( lat%lst(row,col), species )%list )
        if ( n_ads_sites == 1 ) then
          lat%ads_list(new_n_ads)%ast = lat%avail_ads_sites( lat%lst(row,col), species )%list(1)
        else
          lat%ads_list(new_n_ads)%ast = lat%avail_ads_sites( lat%lst(row,col), species )%list( irand(n_ads_sites) )
        end if

        ! energy change due to the added adsorbate minus chemical potential
        delta = energy(new_n_ads, lat, e_pars) - c_pars%chem_pots(species)

        ! If delta is non-positive, accept the addition of adsorbate , i.e., do nothing
        ! this check also prevents overflow in exp
        if ( delta > 0.0_dp ) then
          probability = exp(- beta*delta)
          if ( probability < ran1() ) then
            ! reject the addition of adsorbate
            lat%occupations(row,col) = 0
            lat%ads_list(new_n_ads)%row = 0
            lat%ads_list(new_n_ads)%col = 0
            lat%ads_list(new_n_ads)%ast = 0
            lat%n_ads(species) = lat%n_ads(species)  - 1
          end if
        end if

      else

        ! Decide whether to remove an adsorbate

        ! Save total number of adsorbates
        n_ads_tot_old = lat%n_ads_tot()

        ! Get number of an adsorbate to be removed
        i_ads = lat%occupations(row,col)

        ! chemical potential minus the energy change due to the removal of the adsorbate
        delta = c_pars%chem_pots(species) - energy(i_ads, lat, e_pars)

        ! Prevents overflow in exp
        beta_delta = - beta*delta
        remove = .false.
        if ( beta_delta > exp_arg_too_big ) then
          remove = .true.
        else
          probability = exp(beta_delta)
          if ( probability >= 1.0_dp ) then
            remove = .true.
          else if ( probability > ran1() ) then
            remove = .true.
          end if
        end if

        if (remove) then
          ! Remove adsorbate from the lattice
          lat%occupations(row, col) = 0
          ! Update the number of adsorbates in the lat structure
          lat%n_ads(species) = lat%n_ads(species) - 1
          ! Rearrange adsorbates except when the last adsorbate removed
          if ( i_ads < n_ads_tot_old ) then
            ! Put the last adsorbate in place of removed adsorbate
            lat%ads_list(i_ads) = lat%ads_list(n_ads_tot_old)
            ! Update the adsorbate number in the lattice
            lat%occupations(lat%ads_list(i_ads)%row, lat%ads_list(i_ads)%col) = i_ads
          end if
        end if

      end if

    end if ! gc_period

    ! Accumulate running averages
    if (c_pars%running_avgs_save > 0) then

      do i=1, lat%n_ads_tot()

        row     = lat%ads_list(i)%row
        col     = lat%ads_list(i)%col
        species = lat%ads_list(i)%id
        i_lst   = lat%lst(row,col)
        lst_pops(species,i_lst) = lst_pops(species,i_lst) + 1
        replica_pops(species, (col-1)/c_pars%step_period + 1) = replica_pops(species, (col-1)/c_pars%step_period + 1) + 1

        ! counting zigzags
        if (i_lst == corner_site) then ! an adsorbate at the corner

          ns(1) = lat%occupations(                             row, modulo(col-1-1,c_pars%n_cols)+1 ) ! { 0,-1}
          ns(2) = lat%occupations( modulo(row+2-1,c_pars%n_rows)+1, modulo(col-1-1,c_pars%n_cols)+1 ) ! { 2,-1}
          ns(3) = lat%occupations( modulo(row+2-1,c_pars%n_rows)+1,                             col ) ! { 2, 0}
          ns(4) = lat%occupations( modulo(row-2-1,c_pars%n_rows)+1,                             col ) ! {-2, 0}

          do j=1,4
            if (ns(j) > 0) then
              ns(j) = 2
            else
              ns(j) = 1
            end if
          end do

          zigzag(ns(1),ns(2),ns(3),ns(4)) = zigzag(ns(1),ns(2),ns(3),ns(4)) + 1

        end if

      end do

    end if

    ! Calculate the cluster size histogramm
    if (c_pars%hist_period > 0 .and. mod(istep, c_pars%hist_period) == 0) then
      hist_counter = hist_counter + 1

!      do species=1,c_pars%n_species
!        call lat%hoshen_kopelman(species, cluster_label, largest_label)
!        call lat%cluster_size(species, cluster_label, cluster_sizes)
!        do i=1,largest_label
!          if (cluster_sizes(i) > 0) &
!            hist(species,cluster_sizes(i)) = hist(species,cluster_sizes(i)) + 1
!        end do
!        print*,'species =',species, 'largest_label = ',largest_label
!        print*,'cluster_sizes = ',cluster_sizes
!        print*,'hist:',hist(species,:)
!        print*
!      end do

    end if

    ! Calculate rdf
    if (c_pars%rdf_period > 0 .and. mod(istep, c_pars%rdf_period) == 0) then

      rdf_counter = rdf_counter + 1
      call lat%rdf_hist(c_pars,rdf_hist)

    end if

    if (mod(istep, c_pars%save_period) == 0) then
      if (c_pars%show_progress > 0) call progress_bar(100*istep/c_pars%n_mmc_steps)

      ! Save energy
      write(outeng_unit,*) istep, total_energy(lat,e_pars), lat%n_ads

      ! Save configuration
      if (c_pars%conf_save > 0 ) then
        write(outcfg_unit,'(A10,i0)') "mmc_step ",istep
        write(outcfg_unit,'(100i10)') lat%n_ads
        call lat%print_ads(outcfg_unit)
      end if

      ! Save running averages
      if (c_pars%running_avgs_save > 0) then
        write(outrav_unit,'(i0,1000f15.8)') istep, cov_factor*lst_pops, cov_factor*zigzag, cov_factor*replica_pops
        lst_pops = 0
        replica_pops = 0
        zigzag = 0
      end if

      ! Save cluster size histogram
      if (c_pars%hist_period > 0) then
        write(outhst_unit,*) 'counts ', hist_counter
        do species=1,c_pars%n_species
          write(outhst_unit,*) species
          write(outhst_unit,n_ads_fmt) hist(species,:)
        end do
        hist_counter = 0
        hist = 0
      end if

      ! Save rdf
      if (c_pars%rdf_period > 0) then
        write(outrdf_unit,*) 'counts ', rdf_counter
        do species1=1,c_pars%n_species
        do species2=1,c_pars%n_species
          write(outrdf_unit,*) species1, species2
          write(outrdf_unit,rdf_fmt) rdf_hist(species1,species2,:)
        end do
        end do
        rdf_counter = 0
        rdf_hist = 0
      end if

    end if

  end do ! over mmc steps

  ! Save the last configuration
  if (mod(c_pars%n_mmc_steps, c_pars%save_period) /= 0) then
    write(outcfg_unit,'(A10,i0)') "mmc_step ",c_pars%n_mmc_steps
    write(outcfg_unit,'(100i10)') lat%n_ads
    call lat%print_ads(outcfg_unit)
  end if

  close(outcfg_unit)
  close(outeng_unit)
  if (c_pars%running_avgs_save > 0) close(outrav_unit)
  if (c_pars%hist_period > 0) close(outhst_unit)
  if (c_pars%rdf_period  > 0) close(outrdf_unit)

  print*
  print*,"MMC done. Goodbye."

end subroutine metropolis


end module mmc
