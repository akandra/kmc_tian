module kmc

  use constants
  use utilities
  use mc_lat_class
  use control_parameters_class
  use energy_parameters_class
  use energy_mod
  use rates_hopping_class

  implicit none

  private
  public:: Bortz_Kalos_Lebowitz

contains

subroutine Bortz_Kalos_Lebowitz(lat, c_pars, e_pars)

  type(mc_lat), intent(inout) :: lat
  type(control_parameters), intent(inout) :: c_pars
  type(energy_parameters ), intent(in) :: e_pars

  type(hopping_rates_type) :: r_hop

  real(dp) :: beta
  character(len=max_string_length) :: buffer, n_ads_fmt
  integer :: itraj, ibin, ibin_new, k_change
  integer :: i,m,ads,m_nn, n_nn, n_nn2, iads, species
  integer :: kmc_nsteps
  integer, dimension(lat%n_rows,lat%n_cols) :: cluster_label
  ! Warning: size of cluster_sizes and hist array are too large
  !          consider a way to decrease them
  integer, dimension(lat%n_rows*lat%n_cols) :: cluster_sizes
  integer :: largest_label
  integer, dimension(c_pars%n_species,lat%n_rows*lat%n_cols) :: hist
  integer :: n_ads_total
  integer :: row, col, row_new, col_new, lst_new, ast_new, id
  real(dp) :: energy_old, energy_new
  real(dp) :: time, total_rate, u, rate_acc
  real(dp) :: delta_t, time_new, step_bin
  integer, dimension(lat%n_nn(1),lat%n_nn(1)/2) :: nn_new
  integer, dimension(2*lat%n_nn(1)) :: change_list

  ! Create a rate structure
  r_hop = hopping_rates_init(c_pars, lat, e_pars)
  call r_hop%print(c_pars)

  ! inverse thermodynamic temperature
  beta = 1.0_dp/(kB*c_pars%temperature)
  ! total number of adsorbates
  n_ads_total = lat%n_ads_tot()
  !
  n_nn  = lat%n_nn(1)
  n_nn2 = n_nn/2
  ! List of additional nn directions to scan after hop
  do m=1,n_nn
  do i=1,n_nn2
    nn_new(m,i) = modulo( m+i-n_nn2, n_nn ) + 1
  end do
  end do

  ! time binning for distributions
  step_bin = c_pars%t_end/c_pars%n_bins

  write(*,'(20X,A)') "B.K.L. Code's progress report:"
  ! Loop over trajectories
  do itraj=1, c_pars%n_trajs

    call progress_bar( 'current trajectory ', 0* 100*itraj/c_pars%n_trajs , '   total', 0)

    ibin = 1
    kmc_nsteps = 0

    ! initialize random number generator
    call random_seed(put=itraj*randseed*42)

    write(buffer,'(i6.6)') itraj

    ! write initial state of the lattice to a file
    call open_for_write(outcfg_unit,trim(c_pars%file_name_base)//trim(buffer)//'.confs')
    write(outcfg_unit,'(A10,A10,A15)') &
                  "# rows","# cols","step_period"
    write(outcfg_unit,'(3i10)') lat%n_rows, lat%n_cols, c_pars%step_period
    write(outcfg_unit,'(100A10)') adjustr(c_pars%ads_names)
    write(outcfg_unit,'(A6,1pe12.3)') "time ",0.0_dp
    write(outcfg_unit,'(100i10)') lat%n_ads
    call lat%print_ads(outcfg_unit)
    ! write initial total energy of the system
    call open_for_write(outeng_unit,trim(c_pars%file_name_base)//trim(buffer)//'.en')
    write(outeng_unit,'(1pe12.3,1pe12.3)') 0.0_dp, total_energy(lat,e_pars)
    ! open file for saving cluster size histogram
    call open_for_write(outhst_unit,trim(c_pars%file_name_base)//trim(buffer)//'.hist')
    ! Initialize cluster histogram
    hist = 0

    ! Construct rates array
    do i=1,n_ads_total
      call r_hop%construct(i, lat, e_pars, beta)
    end do

!    call lat%print_ocs
!    do i = 1,n_ads_total
!    do m = 1,n_nn
!      print*, i,m,(r_hop%rates(i,m)%list(iads)==r1%rates(i,m)%list(iads) ,iads=1,size(r%rates(i,m)%list))
!
!    end do
!    end do
!    stop 110

    ! start time propagation
    time = big_bang
    do while (time<c_pars%t_end)

      ! calculate total rate
      total_rate = 0.0_dp
      do ads=1,n_ads_total
      do m=1,n_nn
        total_rate = total_rate + sum(r_hop%rates(ads,m)%list)
      end do
      end do
      !      print*, total_rate

      ! random number to select a process
      u = ran1()*total_rate

      ! determine reaction channel (adsorbate, direction, ads. site of available ones)
      !                            (      ads,      m_nn,                        iads)
      rate_acc = 0.0_dp
      extloop: do ads=1,n_ads_total
        do m_nn=1,n_nn
        do iads=1,size(r_hop%rates(ads,m_nn)%list) ! Warning: check timing of size calculation
          rate_acc = rate_acc + r_hop%rates(ads,m_nn)%list(iads)
          if (u < rate_acc) exit extloop
        end do
        end do
      end do extloop

      delta_t = -log(ran1())/total_rate   ! when does a hop occur?
      time_new = time + delta_t
      kmc_nsteps = kmc_nsteps + 1

      call progress_bar(                  &
          'current trajectory',           &
          int(100*time_new/c_pars%t_end), &
          '   total',                     &
          ! total =  % completed trajs + contribution form current trajctory
          100*(itraj-1)/c_pars%n_trajs + int(100.*time_new/(c_pars%t_end*c_pars%n_trajs)))


!      print*, 'ran. number is ', u/total_rate, 'rate_acc is ',rate_acc
!      print*, 'reaction channel is:', 'ads = ',ads, 'dir = ',m_nn, 'ads. site is ',ast_new
!      print*, 'rate  is ',rates(ads,m_nn)%list(ast_new)
!      stop 321

      if (time_new > c_pars%t_end) time_new = c_pars%t_end

      ibin_new = int(time_new/step_bin) + 1

      if (ibin_new - ibin > 0 ) then

        !-----Save configuration
        write(outcfg_unit,'(A6,1pe12.3)') "time ",ibin*step_bin
        write(outcfg_unit,'(100i10)') lat%n_ads
        call lat%print_ads(outcfg_unit)

        !-----Save energy
        write(outeng_unit,'(1pe12.3,1pe12.3)') ibin*step_bin, total_energy(lat,e_pars)

        !-----Save histogram with cluster sizes
        ! Output formats
        write(n_ads_fmt,'(i10)') n_ads_total
        n_ads_fmt = '('//trim(adjustl(n_ads_fmt))//'i8)'
        ! Calculate histogram
        do species=1,c_pars%n_species
          call lat%hoshen_kopelman(species, cluster_label, largest_label)
          call lat%cluster_size(species, cluster_label, cluster_sizes)
          do i=1,largest_label
            if (cluster_sizes(i) > 0) &
              hist(species,cluster_sizes(i)) = hist(species,cluster_sizes(i)) + 1
          end do
        end do
        ! Save cluster size histogram
        write(outhst_unit,'(A6,1pe12.3)') "time ",ibin*step_bin
        do species=1,c_pars%n_species
          write(outhst_unit,*) species
          write(outhst_unit,n_ads_fmt) hist(species,:)
        end do
        hist = 0

      end if

      ibin = ibin_new
      time = time_new ! time shift

!------------ Update rate constants
!   Particle (ads) is going to hop in direction (m_nn) to ads. site (site)

      ! create a list of adsorbates affected by hop
      change_list = 0
      k_change = 1
      ! Put the hopping particle iads into the list
      change_list(k_change) = ads

      ! scan over old neighbors
      do m=1,n_nn
        ! position of neighbor m
        call lat%neighbor(ads,m,row,col)
        if (lat%occupations(row,col) > 0) then
          k_change = k_change + 1
          change_list(k_change) = lat%occupations(row,col)
        end if
      end do

      ! a new position of particle (ads) after a hop to a neighbor (m_nn)
      call lat%neighbor(ads,m_nn,row_new,col_new)
      id = lat%ads_list(ads)%id
      lst_new = lat%lst(row_new,col_new)
      ast_new = lat%avail_ads_sites(id,lst_new)%list(iads)

      ! Make a hop:
      ! Delete an adsorbate from its old position
      lat%occupations(lat%ads_list(ads)%row, lat%ads_list(ads)%col) = 0
      ! Put the adsorbate in a new position
      lat%ads_list(ads)%row  = row_new
      lat%ads_list(ads)%col  = col_new
      lat%ads_list(ads)%ast  = ast_new
      ! Update adsorbate position
      lat%occupations(row_new, col_new) = ads

      ! scan over additional new neighbors
      do m=1,n_nn2
        ! position of neighbor nn_new(m_nn,m)
        call lat%neighbor( ads, nn_new(m_nn,m), row,col)
        if (lat%occupations(row,col) > 0) then
          k_change = k_change + 1
          change_list(k_change) = lat%occupations(row,col)
        end if
      end do

!      print*,ads, m_nn, &
!        ads_site_names(lat%avail_ads_sites(lat%ads_list(ads)%id,&
!                            lat%lst(row_new, col_new))%list(ast_new))
!      stop
!      call lat%print_ocs
!      do i=1,n_ads_total
!      do m=1,n_nn
!        print'(2i,10e16.4)', i,m,rates(i,m)%list
!      end do
!      end do
!      print*,k_change
!      print*,change_list
!      stop 15

      ! Update rate array for the affected adsorbates
      do i=1,k_change
        call r_hop%construct(change_list(i), lat, e_pars, beta)
      end do

    end do ! over time
!-------------------------------------------------------------

    close(outcfg_unit)
    close(outeng_unit)
    close(outhst_unit)

  end do ! over trajectories

end subroutine Bortz_Kalos_Lebowitz


end module kmc
