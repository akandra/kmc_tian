module kmc

  use constants
  use utilities
  use mc_lat_class
  use control_parameters_class
  use energy_parameters_class
  use energy_mod
  use open_file
  use rates_hopping_class
  use rates_desorption_class
  use reaction_class

  implicit none

  private
  public:: Bortz_Kalos_Lebowitz

contains

subroutine Bortz_Kalos_Lebowitz(lat, c_pars, e_pars)

  type(mc_lat), intent(inout) :: lat
  type(control_parameters), intent(inout) :: c_pars
  type(energy_parameters ), intent(in)    :: e_pars

  type(reaction_type) :: r

  character(len=max_string_length) :: buffer, n_ads_fmt, rdf_fmt, version_header
  integer :: itraj, ibin, ibin_new, ibin_current, bin_shift
  integer :: time_segment, time_segment_old, time_segment_new
  integer :: i, species, species1, species2
  integer :: kmc_nsteps
!dja  integer, dimension(lat%n_rows,lat%n_cols) :: cluster_label
!dja  ! Warning: size of cluster_sizes and hist array are large
!dja   !          consider a way to decrease them
!dja  integer, dimension(lat%n_rows*lat%n_cols) :: cluster_sizes
!dja  integer :: largest_label
  integer, dimension(c_pars%n_species,lat%n_rows*lat%n_cols) :: hist
  real(dp) :: time, end_of_time, time_bin, time_shift
  real(dp) :: delta_t
  real(dp), dimension(size(c_pars%t_end)) :: time_limits, step_bin
  logical :: bin_log_scale, bin_power_scale
  real(dp) :: log_t1, bin_power_law

  integer, dimension(c_pars%n_species,c_pars%n_species,c_pars%rdf_n_bins) :: rdf_hist

!dja  integer :: j, k, n_chan  ! only for debug printout
  integer :: j, k  ! only for debug printout

  version_header = '! kmc_tian Release ' // version

  ! Create a rate structure
  r = reaction_init(c_pars, lat, e_pars)
!------------------------------------------------------------------------------
!  debug printing
!------------------------------------------------------------------------------
!print *, 'debug printouts'
!print '(A, i0)', ' n_processes = ', r%association%n_processes
!print *, ' lat%n_ads   = ', lat%n_ads
!
!print '(A)', ' n_chan'
!do i=1, lat%n_ads(1)
!  n_chan = size(r%association%rate_info(i)%list)
!  write(*, '(i4, i6)',advance='no') i, n_chan
!  print*, r%association%rate_info(i)%list(:)%rate
!
!end do
!
!stop 111

  bin_log_scale   = c_pars%log_scale_t1 > 0.0_dp
  bin_power_scale = c_pars%power_scale_t1  > 0.0_dp
  ! do linear time binning for distributions
  if (bin_log_scale) then
    log_t1 = log10(c_pars%log_scale_t1)
    step_bin(1) = ( log10(c_pars%t_end(1)) - log_t1 )/c_pars%n_bins(1)
  elseif (bin_power_scale) then
    bin_power_law = ( log(c_pars%t_end(1)) - log(c_pars%power_scale_t1) )/ log( 1.0_dp*c_pars%n_bins(1) )
    step_bin(1) = c_pars%power_scale_t1**(1.0/bin_power_law)
  else
    step_bin = c_pars%t_end/c_pars%n_bins
  end if

  ! calculate cumulative times for segments
  time_limits(1) = c_pars%t_end(1)
  do i=2, size(c_pars%t_end)
    time_limits(i) = time_limits(i-1) + c_pars%t_end(i)
  end do

  if (c_pars%show_progress) write(*,'(20X,A)') "B.K.L. Code's progress report:"

  ! Loop over trajectories
  do itraj=c_pars%start_traj, c_pars%start_traj - 1 + c_pars%n_trajs

!    debug(1) = (itraj==3)

    if (c_pars%show_progress) call progress_bar( 'current trajectory', 0 , &
                                 '   total', 100*(itraj-c_pars%start_traj+1)/c_pars%n_trajs )
    ibin = 0
    kmc_nsteps = 0

    ! initialize random number generator
    call random_seed(put=itraj*randseed*42)
    ! Set initial lattice based on trajectory number
    call lat%conf_init(c_pars, itraj )
    ! Adjust the respective records in reaction variable
    r%n_ads_total = lat%n_ads_tot()

    ! Initialize reaction counter
    r%counter = 0
    !-------- Construct rates arrays
    call r%construct(lat, e_pars)

    if (r%n_ads_total>0) then

      ! Calculate initial cluster size histogram
      hist = 0
!      do species=1,c_pars%n_species
!        call lat%hoshen_kopelman(species, cluster_label, largest_label)
!        call lat%cluster_size(species, cluster_label, cluster_sizes)
!        do i=1,largest_label
!          if (cluster_sizes(i) > 0) &
!            hist(species,cluster_sizes(i)) = hist(species,cluster_sizes(i)) + 1
!        end do
!      end do

      ! Calculate initial rdf
      if (c_pars%rdf_period > 0) then
        rdf_hist = 0
        call lat%rdf_hist(c_pars, rdf_hist)
      end if

    end if

    ! write trajectory files
    !-----------------------------------------------------------------
    ! Prepare string with traj number for using in filenames
    write(buffer,'(i9.9)') itraj
    ! write initial state of the lattice into a file
    call open_for_write(outcfg_unit,trim(c_pars%file_name_base)//'_'//trim(buffer)//'.confs')
    write(outcfg_unit,'(A)') trim(version_header)
    write(outcfg_unit,'(A20,A20,A20)') "# rows(|| to steps)","# cols","step_period"
    write(outcfg_unit,'(3i20)') lat%n_rows, lat%n_cols, c_pars%step_period
    write(outcfg_unit,'(100A10)') adjustr(c_pars%ads_names)
    write(outcfg_unit,'(A6,1pe19.10)') "time ",0.0_dp
    write(outcfg_unit,'(100i10)') lat%n_ads
    call lat%print_ads(outcfg_unit)

    ! write initial total energy of the system
    call open_for_write(outeng_unit,trim(c_pars%file_name_base)//'_'//trim(buffer)//'.en')
    write(outeng_unit,'(2A12)') 'time(s)', ' energy(eV)'
    write(outeng_unit,'(1pe19.10,1pe19.10)') 0.0_dp, total_energy(lat,e_pars)

    ! open file for saving reaction counts
    call open_for_write(outcnt_unit,trim(c_pars%file_name_base)//'_'//trim(buffer)//'.counts')
    write(outcnt_unit,'(A12,100A20)') 'time(s)', &
    ((' '//trim(c_pars%ads_names(j))//'_'//trim(reaction_names(k)), k=1,n_reaction_types), j=1,c_pars%n_species)
    write(outcnt_unit,'(1pe19.10,100i20)') 0.0_dp, r%counter

    ! open file for saving cluster size histogram
    call open_for_write(outhst_unit,trim(c_pars%file_name_base)//'_'//trim(buffer)//'.hist')

    ! open file for saving rdf and define rdf-specific variables
    if (c_pars%rdf_period > 0) then
      call open_for_write(outrdf_unit,trim(c_pars%file_name_base)//'_'//trim(buffer)//'.rdf')
      ! Output formats
      write(rdf_fmt,'(i10)') c_pars%rdf_n_bins
      rdf_fmt = '('//trim(adjustl(rdf_fmt))//'i8)'
    end if

    if (r%n_ads_total>0) then

    !-----Save histogram with cluster sizes

      ! Output formats
      write(n_ads_fmt,'(i10)') r%n_ads_total
      n_ads_fmt = '('//trim(adjustl(n_ads_fmt))//'i8)'
      ! Save cluster size histogram
      write(outhst_unit,'(A6,1pe19.10)') "time ",0.0_dp
      do species=1,c_pars%n_species
        write(outhst_unit,*) species
        write(outhst_unit,n_ads_fmt) (hist(species,i), i=1,lat%n_ads(species))
      end do

    !-----Save rdf

      if (c_pars%rdf_period > 0) then
        ! Write it out
        write(outrdf_unit,'(A6,1pe12.3)') "time ",0.0_dp
        do species1=1,c_pars%n_species
        do species2=1,c_pars%n_species
          write(outrdf_unit,*) species1, species2
          write(outrdf_unit,rdf_fmt) rdf_hist(species1,species2,:)
        end do
        end do
      end if

    else
      ! No output when no adsorbates
      write(outhst_unit,'(A6,1pe19.10)') "time ",0.0_dp
      write(outrdf_unit,'(A6,1pe19.10)') "time ",0.0_dp

    end if

    ! set the time span
    time = big_bang
    end_of_time = time_limits( size(time_limits) )
    time_segment_old = 1
    time_segment_new = 1
    bin_shift = 0
    time_shift = 0.0_dp

    ! start time propagation
    time_loop: do while ( time < end_of_time )

      ! propagate the time
      delta_t = -log(ran1())/r%acc_rate(n_reaction_types)   ! when does a reaction occur?
      !write(*,*) 'Total rate: ', r%acc_rate(n_reaction_types), 'time step: ', delta_t

      time = time + delta_t

      ! set time to the final value at the end of trajectory
      if (time > end_of_time) time = end_of_time

      ! Conditionally jump over time segments
      do while ( time > time_limits(time_segment_new) )
        time_segment_new = time_segment_new + 1
      end do

      do time_segment = time_segment_old, time_segment_new

        if (time_segment == time_segment_new) then

          if (bin_log_scale) then
            ! Exclude possible negative values of ibin_new for small kmc_times
            if (time <= c_pars%log_scale_t1) then
              ibin_new = 0
            else
              ibin_new = int( ( log10(time) - log_t1 )/step_bin(time_segment) )
            end if
          elseif (bin_power_scale) then
            ibin_new = int( time**(1./bin_power_law)/step_bin(time_segment) )
          else
            ibin_new = bin_shift + int( (time - time_shift)/step_bin(time_segment) )
          end if

        else

          time_shift = time_limits(time_segment)
          bin_shift  = bin_shift + c_pars%n_bins(time_segment)
          ibin_new = bin_shift

        end if

        ! Write when switched to the next bin taking into account jumps over more than 1 bin
        ! by repeating output values
        do ibin_current=ibin+1,ibin_new

          if (c_pars%show_progress) call progress_bar(&
            'current trajectory',&
            int(100*time/end_of_time),&
            '   total',&
            ! total =  % completed trajs + contribution form current trajectory
            100*(itraj-c_pars%start_traj)/c_pars%n_trajs + int(100.*time/(end_of_time*c_pars%n_trajs)))

          if (r%n_ads_total>0) then

            ! Reset the histogram
            hist = 0
            ! Calculate histogram
!            do species=1,c_pars%n_species
!              call lat%hoshen_kopelman(species, cluster_label, largest_label)
!              call lat%cluster_size(species, cluster_label, cluster_sizes)
!              do i=1,largest_label
!                if (cluster_sizes(i) > 0) &
!                  hist(species,cluster_sizes(i)) = hist(species,cluster_sizes(i)) + 1
!              end do
!            end do

            if (c_pars%rdf_period > 0) then
              ! Calculate rdf
              rdf_hist = 0
              call lat%rdf_hist(c_pars, rdf_hist)
            end if

          end if

          ! write into trajectory files
          if (bin_log_scale) then
            time_bin = 10**( log_t1 + ibin_current*step_bin(time_segment) )
          elseif (bin_power_scale) then
            time_bin = ( ibin_current*step_bin(time_segment) )**bin_power_law
          else
            time_bin = time_shift + (ibin_current - bin_shift)*step_bin(time_segment)
          end if

          ! configuration
          write(outcfg_unit,'(A6,1pe19.10)') "time ", time_bin
          write(outcfg_unit,'(100i10)') lat%n_ads
          call lat%print_ads(outcfg_unit)

          ! energy
          write(outeng_unit,'(1pe19.10,1pe19.10)') time_bin, total_energy(lat,e_pars)

          ! reaction counts
          write(outcnt_unit,'(1pe19.10,100i20)') time_bin, r%counter
          ! Reset reaction counters
!          r%counter = 0

          if (r%n_ads_total>0) then

            ! histogram with cluster sizes
            ! Output formats
            write(n_ads_fmt,'(i10)') r%n_ads_total
            n_ads_fmt = '('//trim(adjustl(n_ads_fmt))//'i8)'
            ! Save cluster size histogram
            write(outhst_unit,'(A6,1pe19.10)') "time ",time_bin
            do species=1,c_pars%n_species
              write(outhst_unit,*) species
              write(outhst_unit,n_ads_fmt) (hist(species,i), i=1,lat%n_ads(species))
            end do

            ! rdf
            if (c_pars%rdf_period > 0) then
              write(outrdf_unit,'(A6,1pe19.10)') "time ",time_bin
              do species1=1,c_pars%n_species
              do species2=1,c_pars%n_species
                write(outrdf_unit,*) species1, species2
                write(outrdf_unit,rdf_fmt) rdf_hist(species1,species2,:)
              end do
              end do
            end if

          else

            ! No histogram output when no adsorbates
            write(outhst_unit,'(A6,1pe19.10)') "time ",time_bin
            write(outrdf_unit,'(A6,1pe19.10)') "time ",time_bin

          end if

        end do ! ibin

        ibin = ibin_new

      end do ! time segment

      time_segment_old = time_segment_new

      ! do reaction based on a random number
      call r%do_reaction(ran1(), lat, e_pars)

      ! end trajectory if all essential species are gone
      if (.not. any(e_pars%stopping_trigger)) then
        print*
        print *, 'all ' // trim(stopping_trigger_name) // ' species have left the surface: exiting kMC loop'
!        print*
        exit
      end if

      ! end trajectory if there's no processes left
      if (r%acc_rate(n_reaction_types) == 0.0_dp) then
        print*
        print *, 'nothing to do, since total rate is zero: exiting kMC loop'
!        print*
        exit
      end if

      kmc_nsteps = kmc_nsteps + 1

    end do time_loop
!-------------------------------------------------------------


    close(outcfg_unit)
    close(outeng_unit)
    close(outhst_unit)
    if (c_pars%rdf_period > 0) close(outrdf_unit)
    close(outcnt_unit)

  end do ! over trajectories

  print*
  print*,"kMC done. Goodbye."

end subroutine Bortz_Kalos_Lebowitz


end module kmc
