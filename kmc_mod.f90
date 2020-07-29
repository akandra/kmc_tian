module kmc

  use constants
  use utilities
  use mc_lat_class
  use control_parameters_class
  use energy_parameters_class
  use energy_mod
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

  type(reaction_type) :: r, r_save

  character(len=max_string_length) :: buffer, n_ads_fmt
  integer :: itraj, ibin, ibin_new
  integer :: i, species
  integer :: kmc_nsteps
  integer, dimension(lat%n_rows,lat%n_cols) :: cluster_label
  ! Warning: size of cluster_sizes and hist array are large
  !          consider a way to decrease them
  integer, dimension(lat%n_rows*lat%n_cols) :: cluster_sizes
  integer :: largest_label
  integer, dimension(c_pars%n_species,lat%n_rows*lat%n_cols) :: hist
  real(dp) :: time
  real(dp) :: delta_t, time_new, step_bin

  integer :: j, k, n_chan  ! only for debug printout

  ! initialize vector of conditions for debug trap
  debug =  .false.

  ! Create a rate structure
  r = reaction_init(c_pars, lat, e_pars)

!!------------------------------------------------------------------------------
!!  debug printing
!!------------------------------------------------------------------------------
!print *, 'debug printouts'
!print '(A, i0)', ' n_processes = ', r%dissociation%n_processes
!print '(A, i0)', ' lat%n_ads   = ', lat%n_ads
!
!print '(A)', ' n_chan'
!do i=1, lat%n_ads(1)
!  n_chan = size(r%dissociation%rate_info(i)%list)
!  write(*, '(i4, i6)',advance='no') i, n_chan
!  print*, r%dissociation%rate_info(i)%list(1)%rate
!
!end do
!
!stop 111


  ! time binning for distributions
  step_bin = c_pars%t_end/c_pars%n_bins

  write(*,'(20X,A)') "B.K.L. Code's progress report:"

  ! Loop over trajectories
  do itraj=1, c_pars%n_trajs

    debug(1) = (itraj==2)

    call progress_bar( 'current trajectory ', 0* 100*itraj/c_pars%n_trajs , '   total', 0)

    ibin = 0
    kmc_nsteps = 0

    ! initialize random number generator
    call random_seed(put=itraj*randseed*42)

    write(buffer,'(i6.6)') itraj

    ! Restore initial lattice and adjust the respective records in reaction variable
    call lat%restore
    r%n_ads_total = lat%n_ads_tot()
    r%counter = 0

    !-------- Construct rates arrays
    call r%construct(lat, e_pars)

    ! write initial state of the lattice into a file
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
    write(outeng_unit,'(2A12)') 'time(s)', ' energy(eV)'
    write(outeng_unit,'(1pe12.3,1pe12.3)') 0.0_dp, total_energy(lat,e_pars)
    ! open file for saving reaction counts
    call open_for_write(outcnt_unit,trim(c_pars%file_name_base)//trim(buffer)//'.counts')
    write(outcnt_unit,'(A12,100A20)') 'time(s)', &
      ((' '//trim(c_pars%ads_names(j))//'_'//trim(reaction_names(k)), k=1,n_reaction_types), j=1,c_pars%n_species)
    ! Initialize reaction counters
    r%counter = 0
    ! open file for saving cluster size histogram
    call open_for_write(outhst_unit,trim(c_pars%file_name_base)//trim(buffer)//'.hist')
    !-----Save histogram with cluster sizes
    if (r%n_ads_total>0) then
    ! Output formats
      write(n_ads_fmt,'(i10)') r%n_ads_total
      n_ads_fmt = '('//trim(adjustl(n_ads_fmt))//'i8)'
      ! Calculate histogram
      hist = 0
      do species=1,c_pars%n_species
        call lat%hoshen_kopelman(species, cluster_label, largest_label)
        call lat%cluster_size(species, cluster_label, cluster_sizes)
        do i=1,largest_label
          if (cluster_sizes(i) > 0) &
            hist(species,cluster_sizes(i)) = hist(species,cluster_sizes(i)) + 1
        end do
      end do
      ! Save cluster size histogram
      write(outhst_unit,'(A6,1pe12.3)') "time ",0.0_dp
      do species=1,c_pars%n_species
        write(outhst_unit,*) species
        write(outhst_unit,n_ads_fmt) (hist(species,i), i=1,lat%n_ads(species))
      end do
    else
      ! No histogram output when no adsorbates
      write(outhst_unit,'(A6,1pe12.3)') "time ",0.0_dp
    end if
    ! start time propagation
    time = big_bang
    do while (time<c_pars%t_end)

      ! do reaction based on a random number
      call r%do_reaction(ran1(), lat, e_pars)

      ! End the trajectory if there processes left
      if (r%total_rate == 0.0_dp) then
        print*
        print *, 'total rate is zero: exiting kMC loop'
        print*
        exit
      end if

      ! update time
      delta_t = -log(ran1())/r%total_rate   ! when does a reaction occur?
      time_new = time + delta_t
      kmc_nsteps = kmc_nsteps + 1

      if (time_new > c_pars%t_end) time_new = c_pars%t_end


      call progress_bar(                  &
          'current trajectory',           &
          int(100*time_new/c_pars%t_end), &
          '   total',                     &
          ! total =  % completed trajs + contribution form current trajctory
          100*(itraj-1)/c_pars%n_trajs + int(100.*time_new/(c_pars%t_end*c_pars%n_trajs)))

      ibin_new = int(time_new/step_bin)

      if (ibin_new - ibin > 0 ) then

        !-----Save configuration
        write(outcfg_unit,'(A6,1pe12.3)') "time ",ibin_new*step_bin
        write(outcfg_unit,'(100i10)') lat%n_ads
        call lat%print_ads(outcfg_unit)

        !-----Save energy
        write(outeng_unit,'(1pe12.3,1pe12.3)') ibin_new*step_bin, total_energy(lat,e_pars)

        !-----Save reaction counts
        write(outcnt_unit,'(1pe12.3,100i20)') ibin_new*step_bin, r%counter
        ! Reset reaction counters
        r%counter = 0

        if (r%n_ads_total>0) then
          !-----Save histogram with cluster sizes
          ! Output formats
          write(n_ads_fmt,'(i10)') r%n_ads_total
          n_ads_fmt = '('//trim(adjustl(n_ads_fmt))//'i8)'
          ! Reset the histogram
          hist = 0
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
          write(outhst_unit,'(A6,1pe12.3)') "time ",ibin_new*step_bin
          do species=1,c_pars%n_species
            write(outhst_unit,*) species
            write(outhst_unit,n_ads_fmt) (hist(species,i), i=1,lat%n_ads(species))
          end do

        else
          ! No histogram output when no adsorbates
          write(outhst_unit,'(A6,1pe12.3)') "time ",ibin_new*step_bin

        end if

      end if

      ibin = ibin_new
      time = time_new ! time shift

    end do ! over time
!-------------------------------------------------------------

    close(outcfg_unit)
    close(outeng_unit)
    close(outhst_unit)
    close(outcnt_unit)

  end do ! over trajectories

end subroutine Bortz_Kalos_Lebowitz


end module kmc
