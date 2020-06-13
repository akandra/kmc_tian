module kmc

  use constants
  use utilities
  use mc_lat_class
  use control_parameters_class
  use energy_parameters_class
  use energy_mod
  use rates_class

  implicit none

  private
  public:: Bortz_Kalos_Lebowitz

  type, private :: v_list_dp

    real(dp), dimension(:),allocatable :: list

  end type

contains

subroutine Bortz_Kalos_Lebowitz(lat, c_pars, e_pars)

  type(mc_lat), intent(inout) :: lat
  type(control_parameters), intent(inout) :: c_pars
  type(energy_parameters ), intent(in) :: e_pars

  type( rates_type)        :: r
  ! rates array (n_adsorbates x n_neighbors)
  type(v_list_dp), dimension(:,:), allocatable :: rates

  real(dp) :: beta
  character(len=max_string_length) :: n_ads_fmt
  integer :: itraj, step_bin, ibin, k_change
  integer :: i,m,ads,m_nn, n_nn, n_nn2
  integer :: kmc_nsteps
  integer :: largest_label, hist_counter
  integer, dimension(c_pars%n_species,maxval(lat%n_ads)) :: hist
  integer :: n_ads_total
  integer :: row_old, col_old, st_old, ast_old
  integer :: row_new, col_new, st_new, ast_new
  integer :: row, col, site, id
  integer, dimension(lat%n_nn(1)) :: nn_row, nn_col
  real(dp) :: energy_old, energy_new
  real(dp) :: time, total_rate, u, rate_acc
  real(dp) :: delta_t, time_new
  integer, dimension(lat%n_nn(1),lat%n_nn(1)/2) :: nn_new
  integer, dimension(2*lat%n_nn(1)) :: change_list

  ! Create a rate structure
  r = rates_init(c_pars, lat, e_pars)
  call r%print_r_hop(c_pars)

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

  ! Output formats
  write(n_ads_fmt,'(i10)') lat%n_ads_tot()
  n_ads_fmt = '('//trim(adjustl(n_ads_fmt))//'i8)'

  ! write initial state of the lattice to a file
  call open_for_write(outcfg_unit,trim(c_pars%file_name_base)//'.confs')

  write(outcfg_unit,'(A10,A10,A15)') &
                "# rows","# cols","step_period"
  write(outcfg_unit,'(3i10)') lat%n_rows, lat%n_cols, c_pars%step_period
  write(outcfg_unit,'(100A10)') adjustr(c_pars%ads_names)
  write(outcfg_unit,'(A10,i0)') "bkl step ", 0
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

  write(*,'(20X,A)') "B.K.L. Code's progress report:"
  call progress_bar(0)

  ! Allocate rates array
  allocate( rates(n_ads_total,n_nn) )
  do i=1,n_ads_total
    row = lat%ads_list(i)%row
    col = lat%ads_list(i)%col
    site = lat%site_type(row,col)
    id  = lat%ads_list(i)%id
    do m=1,n_nn
      allocate( rates(i,m)%list(size(lat%avail_ads_sites(id,site)%list)) )
    end do
  end do

  ! time binning for distributions
  if (c_pars%n_bins > 0) step_bin = c_pars%t_end/c_pars%n_bins

  do itraj=1, c_pars%n_trajs

    print*, 'Running trajectory no.',itraj
    print*
    ibin = 1
    kmc_nsteps = 0

    ! Warning: check the meaning of arguments
    ! initialize random number generator
!    call random_seed(size=i)
    call random_seed(put=itraj*randseed)

    !write(buffer,'(i6.6)') itraj
!            if (n_bins > 0) then
!
!                call open_for_write(outeng_unit,trim(fname)//trim(buffer)//'.en')
!                write(outeng_unit,*) big_bang ,&
!                    total_energy(nlat, nads, nnn, occupations, site_type, &
!                                 ads_list, nn_list, ads_energy, int_energy)/eV2K
!
!                call open_for_write(outhst_unit,trim(fname)//trim(buffer)//'.csz')
!                write(outhst_unit,*) t_end, n_bins
!
!            end if

    ! Construct rates array
    do i=1,n_ads_total

      ! energy for particle i in its old position
      energy_old = energy(i, lat, e_pars)

      ! Save the old configuration
      row_old  = lat%ads_list(i)%row
      col_old  = lat%ads_list(i)%col
      st_old   = lat%site_type(row_old,col_old)
      ast_old  = lat%ads_list(i)%site
      id       = lat%ads_list(i)%id

      ! Delete particle i from the old position
      lat%occupations(row_old,col_old) = 0

      ! Loop over possible new positions of particle i
      do m=1,n_nn

        ! Get position and site type of neighbour m
        call lat%neighbor(i, m, row_new, col_new)
        st_new  = lat%site_type(row_new, col_new)

        ! Check if the cell is free
        if (lat%occupations(row_new, col_new) > 0) then

          rates(i,m)%list = 0.0d0

        else

          ! Put particle i to site m
          lat%ads_list(i)%row = row_new
          lat%ads_list(i)%col = col_new
          lat%occupations(row_new,col_new) = i

          ! Loop over adsorption site
          do ads = 1, size(lat%avail_ads_sites(id,st_new)%list)

            ! Move particle i to adsorption site list(iads)
            ast_new = lat%avail_ads_sites(id,st_new)%list(ads)
            lat%ads_list(i)%site = ast_new

            ! Calculate energy of i in new position
            energy_new = energy(i, lat, e_pars)

            ! Apply detailed balance when
            ! energy in the old position < energy in the new position
            if (energy_old < energy_new) then
                rates(i,m)%list(ads) = r%r_hop(id, st_old, ast_old, st_new, ast_new)&
                    *exp( -beta*(energy_new - energy_old) )
!                print*
!                print*, 'id ',id, ' old site ',st_old,' old ads. site ', ast_old
!                print*, ' new site ',st_new,' new ads. site ', ast_new
!                print*, 'rate ',r%r_hop(id, st_old, ast_old, st_new, ast_new)
            else
                rates(i,m)%list(ads) = r%r_hop(id, st_old, ast_old, st_new, ast_new)
!                print*, id,st_old, ast_old, st_new, ast_new,r%r_hop(id, st_old, ast_old, st_new, ast_new)
            end if

!            print '(A,i4,A,i4,A,i4)', 'ads ',i,' neighbor ',m,' ads. site ',lat%ads_list(i)%site
!            print '(A,e18.4,A,e18.4)',' E_old = ', energy_old, ' E_new = ',energy_new
!            print *,' r_hop = ',r%r_hop(id, st_old, ast_old, st_new, ast_new),&
!                                      ' rate = ', rates(i,m)%list(iads)
!            write(*,*) 'pause'
!            read(*,*)

          end do ! iads

          ! Delete particle i from site m
          lat%occupations(row_new,col_new) = 0

        end if

      end do ! m

      ! Return particle i to the old position
      lat%ads_list(i)%row  = row_old
      lat%ads_list(i)%col  = col_old
      lat%ads_list(i)%site = ast_old
      lat%occupations(row_old,col_old) = i

    end do ! i

    call lat%print_ocs
!    do i = 1,n_ads_total
!    do m = 1,n_nn
!      print*, i,m,(rates(i,m)%list(iads),iads=1,size(rates(i,m)%list))
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
        total_rate = total_rate + sum(rates(ads,m)%list)
      end do
      end do
      !      print*, total_rate

      ! random number to select a process
      u = ran1()*total_rate

      ! determine reaction channel (adsorbate, direction, ads. site)
      !                            (      ads,      m_nn,   ast_new)
      rate_acc = 0.0_dp
      extloop: do ads=1,n_ads_total
        do m_nn=1,n_nn
        do ast_new=1,size(rates(ads,m_nn)%list) ! Warning: check timing of size calculation
          rate_acc = rate_acc + rates(ads,m_nn)%list(ast_new)
          if (u < rate_acc) exit extloop
        end do
        end do
      end do extloop

      delta_t = -log(ran1())/total_rate   ! when does a hop occur?
      time_new = time + delta_t
      kmc_nsteps = kmc_nsteps + 1
!      print*, 'ran. number is ', u/total_rate, 'rate_acc is ',rate_acc
!      print*, 'reaction channel is:', 'ads = ',ads, 'dir = ',m_nn, 'ads. site is ',site
!      stop 321

!      if (time_new > t_end) time_new = t_end

!                if(n_bins > 0) then
!
!
!                    ibin_new = int(time_new/step_bin) + 1
!
!                    if (ibin_new - ibin > 0 ) then
!                        call hoshen_kopelman(cluster_label, largest_label, occupations, ads_list, &
!                                                nn_list, nads, nlat, nnn)
!                        call count_cluster_sizes(cluster_sizes, cluster_label,&
!                                                                ads_list, nads, nlat)
!                        hist = 0
!                        do i=1,largest_label
!                            if (cluster_sizes(i) > 0) &
!                                hist(cluster_sizes(i)) = hist(cluster_sizes(i)) + 1
!                        end do
!
!                        write(outhst_unit,*) ibin
!                        write(outhst_unit,ads_fmt) hist
!
!!                        write(6,cfg_fmt) (occupations(i,:), i=1,nlat)
!                        write(outeng_unit,*) time,&
!                            total_energy(nlat, nads, nnn, occupations, site_type, &
!                                         ads_list, nn_list, ads_energy, int_energy)/eV2K
!
!                    end if
!
!                    ibin = ibin_new
!
!                end if
!
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

      ! Make a hop:
      ! Delete an adsorbate from its old position
      lat%occupations(lat%ads_list(ads)%row, lat%ads_list(ads)%col) = 0
      ! Put the adsorbate in a new position
      lat%ads_list(ads)%row  = row_new
      lat%ads_list(ads)%col  = col_new
      lat%ads_list(ads)%site = ast_new
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

      print*,ads, m_nn, &
        ads_site_names(lat%avail_ads_sites(lat%ads_list(ads)%id,&
                            lat%site_type(row_new, col_new))%list(ast_new))
!      call lat%print_ocs
!      do i=1,n_ads_total
!      do m=1,n_nn
!        print'(2i,10e16.4)', i,m,rates(i,m)%list
!      end do
!      end do
!      print*,k_change
!      print*,change_list
!      stop 15

We are here there and everywere

      ! Update rate array
      ! see building up the rate array above
      ! with i replaced by change_list(i)

!                do kk=1,k_change
!
!                    i = change_list(kk)
!                    i_old = ads_list(i,1)
!                    j_old = ads_list(i,2)
!                    st_old = site_type(i_old,j_old)
!
!                    energy_acc_old = ads_energy(st_old)
!
!                    ! Calculate the interaction energy for particle i in its old positions
!                    ! Loop over the neighbors of particle i
!                    do m=1,nnn
!
!                        nn_pos(m,1) = modulo(i_old + nn_list(m,1)-1,nlat) + 1
!                        nn_pos(m,2) = modulo(j_old + nn_list(m,2)-1,nlat) + 1
!
!                        i_new = nn_pos(m,1)
!                        j_new = nn_pos(m,2)
!                        st_new  = site_type(i_new,j_new)
!
!                        if (occupations(i_new, j_new) > 0)&
!                            energy_acc_old = energy_acc_old + int_energy(st_old, st_new)
!                    end do
!
!                    ! Calculate the interaction energy for particle i in its NEW positions
!
!                    ! Loop over possible new positions of particle i
!                    do m=1,nnn
!
!                        i_new = nn_pos(m,1)
!                        j_new = nn_pos(m,2)
!                        st_new  = site_type(i_new,j_new)
!
!                        if (occupations(i_new, j_new) == 0) then
!
!                            ! Excluding self-counting
!                            energy_acc_new = ads_energy(st_new) - int_energy(st_old,st_new)
!
!                            ! Loop over the neighbors of particle i in its possible NEW positions
!                            do m2=1,nnn
!
!                                i_new2 = modulo(i_new + nn_list(m2,1)-1,nlat) + 1
!                                j_new2 = modulo(j_new + nn_list(m2,2)-1,nlat) + 1
!                                st_new2= site_type(i_new2,j_new2)
!
!                !                print*, 'i, i_old, j_old, st_old :',i, i_old, j_old, st_old
!                !                print*, 'm, i_new, j_new, st_new :',m, i_new, j_new, st_new
!                !                print*, 'm2,i_new2,j_new2,st_new2:',m2, i_new2,j_new2, st_new2
!                !                pause
!
!                                if (occupations(i_new2,j_new2) > 0) &
!                                    energy_acc_new = energy_acc_new + int_energy(st_new,st_new2)
!                            end do
!
!                            ! Apply detailed balance when
!                            ! total energy in the old position < total energy in the new position
!                            if (energy_acc_old < energy_acc_new) then
!                                rates(i,m) = r_hop(st_old,st_new)&
!                                    *exp( -beta*(energy_acc_new - energy_acc_old) )
!                            else
!                                rates(i,m) = r_hop(st_old,st_new)
!                            end if
!
!                        else ! if the neighbor site is occupied
!                            rates(i,m) = 0.0d0
!                        end if
!
!!                        write(*,cfg_fmt) transpose(occupations)
!!                        print*
!!                        print*, 'i, i_old, j_old, st_old:',i, i_old, j_old, st_old
!!                        print*, 'm, i_new, j_new, st_new:',m, i_new, j_new, st_new
!!                        print'(A,f18.2)', 'rates(i,m):',rates(i,m)
!!                        print'(A,f18.4)', 'Eacc_old  = ',energy_acc_old/eV2K
!!                        print'(A,f18.4)', 'E_acc_new = ',energy_acc_new/eV2K
!!                        pause
!                    end do
!
!                end do
!
    end do ! over time
!!-------------------------------------------------------------
!
!            if (n_bins > 0) then
!                close(outeng_unit)
!                close(outhst_unit)
!            end if
!
!            print*,'Number of kmc-steps is ', kmc_nsteps
!            print*
!
  end do ! over trajectories
!
!        deallocate(rates)
!
!    case default
!        stop 'Error: mc algorithm not defined'


  close(outcfg_unit)
  close(outeng_unit)
  if (c_pars%hist_period > 0) close(outhst_unit)

end subroutine Bortz_Kalos_Lebowitz


end module kmc
