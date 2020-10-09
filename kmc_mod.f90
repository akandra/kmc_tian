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

  type(reaction_type) :: r

  character(len=max_string_length) :: buffer, n_ads_fmt, rdf_fmt
  integer :: itraj, ibin, ibin_new, ibin_current
  integer :: i, species, species1, species2
  integer :: kmc_nsteps
  integer, dimension(lat%n_rows,lat%n_cols) :: cluster_label
  ! Warning: size of cluster_sizes and hist array are large
  !          consider a way to decrease them
  integer, dimension(lat%n_rows*lat%n_cols) :: cluster_sizes
  integer :: largest_label
  integer, dimension(c_pars%n_species,lat%n_rows*lat%n_cols) :: hist
  real(dp) :: time
  real(dp) :: delta_t, time_new, step_bin

  integer, dimension(c_pars%n_species,c_pars%n_species,c_pars%rdf_n_bins) :: rdf_hist

  real(dp), dimension(                                                    c_pars%n_bins+1) :: acc_en
  integer,  dimension(size(c_pars%n_ads),                                 c_pars%n_bins+1) :: acc_nads
  integer,  dimension(n_reaction_types,c_pars%n_species,                  c_pars%n_bins+1) :: acc_rcounts
  integer,  dimension(c_pars%n_species,lat%n_rows*lat%n_cols,             c_pars%n_bins+1) :: acc_hist
  integer,  dimension(c_pars%n_species,c_pars%n_species,c_pars%rdf_n_bins,c_pars%n_bins+1) :: acc_rdf
  real(dp) :: energy, n_trajs_inv
  integer, dimension(c_pars%n_species)  :: n_ads_max

  integer :: j, k, n_chan  ! only for debug printout

  ! initialize vector of conditions for debug trap
  debug =  .false.

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

  n_ads_max = 0
  ! initialize accumulators
  acc_en      = 0.0_dp
  acc_nads    = 0
  acc_rcounts = 0
  acc_hist    = 0
  acc_rdf     = 0

  ! time binning for distributions
  step_bin = c_pars%t_end/c_pars%n_bins

  write(*,'(20X,A)') "B.K.L. Code's progress report:"

  ! Loop over trajectories
  do itraj=1, c_pars%n_trajs

    debug(1) = (itraj==3)

    call progress_bar( 'current trajectory ', 0* 100*itraj/c_pars%n_trajs , '   total', 0)

    ibin = 0
    kmc_nsteps = 0

    ! initialize random number generator
    call random_seed(put=itraj*randseed*42)

    write(buffer,'(i6.6)') itraj

    ! Set initial lattice
    call lat%conf_init(c_pars, mod(itraj-1,lat%n_ini_confs)+1 )
    ! Adjust the respective records in reaction variable
    r%n_ads_total = lat%n_ads_tot()

    ! Initialize reaction counter
    r%counter = 0
    !-------- Construct rates arrays
    call r%construct(lat, e_pars)

    ! Update accumulators for initial state

    ! number of adsorbates
    acc_nads(:,1) = acc_nads(:,1) + lat%n_ads

    ! energy
    energy = total_energy(lat,e_pars)
    acc_en(1) = acc_en(1) + energy

    ! reaction counts
    acc_rcounts(:,:,1) = acc_rcounts(:,:,1) + r%counter

    if (r%n_ads_total>0) then

    ! cluster size histogram
    ! Calculate initial cluster size histogram
      hist = 0
      do species=1,c_pars%n_species
        call lat%hoshen_kopelman(species, cluster_label, largest_label)
        call lat%cluster_size(species, cluster_label, cluster_sizes)
        do i=1,largest_label
          if (cluster_sizes(i) > 0) &
            hist(species,cluster_sizes(i)) = hist(species,cluster_sizes(i)) + 1
        end do
      end do
      ! update
      acc_hist(:,:,1) = acc_hist(:,:,1) + hist

      ! rdf
      ! Calculate initial rdf
      if (c_pars%rdf_period > 0) then
        rdf_hist = 0
        call lat%rdf_hist(c_pars, rdf_hist)
        ! update
        acc_rdf(:,:,:,1) = acc_rdf(:,:,:,1) + rdf_hist
      end if

    end if

    ! write trajectory files if needed
    if (c_pars%output_key /= output_key_avg) then

      ! write initial state of the lattice into a file
      call open_for_write(outcfg_unit,trim(c_pars%file_name_base)//'_'//trim(buffer)//'.confs')
      write(outcfg_unit,'(A10,A10,A15)') "# rows","# cols","step_period"
      write(outcfg_unit,'(3i10)') lat%n_rows, lat%n_cols, c_pars%step_period
      write(outcfg_unit,'(100A10)') adjustr(c_pars%ads_names)
      write(outcfg_unit,'(A6,1pe12.3)') "time ",0.0_dp
      write(outcfg_unit,'(100i10)') lat%n_ads
      call lat%print_ads(outcfg_unit)

      ! write initial total energy of the system
      call open_for_write(outeng_unit,trim(c_pars%file_name_base)//'_'//trim(buffer)//'.en')
      write(outeng_unit,'(2A12)') 'time(s)', ' energy(eV)'
      write(outeng_unit,'(1pe12.3,1pe12.3)') 0.0_dp, energy

      ! open file for saving reaction counts
      call open_for_write(outcnt_unit,trim(c_pars%file_name_base)//'_'//trim(buffer)//'.counts')
      write(outcnt_unit,'(A12,100A20)') 'time(s)', &
      ((' '//trim(c_pars%ads_names(j))//'_'//trim(reaction_names(k)), k=1,n_reaction_types), j=1,c_pars%n_species)
      write(outcnt_unit,'(1pe12.3,100i20)') 0.0_dp, r%counter

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
        write(outhst_unit,'(A6,1pe12.3)') "time ",0.0_dp
        do species=1,c_pars%n_species
          write(outhst_unit,*) species
          write(outhst_unit,n_ads_fmt) (hist(species,i), i=1,lat%n_ads(species))
          if (n_ads_max(species) < lat%n_ads(species)) n_ads_max(species) = lat%n_ads(species)
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
        write(outhst_unit,'(A6,1pe12.3)') "time ",0.0_dp
        write(outrdf_unit,'(A6,1pe12.3)') "time ",0.0_dp

      end if

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

      ibin_new = int(time_new/step_bin)

      ! Do write when switched to the next bin taking into account jumps over more than 1 by repeating output values
      do ibin_current=ibin+1,ibin_new

        call progress_bar(                  &
          'current trajectory',           &
          int(100*time_new/c_pars%t_end), &
          '   total',                     &
          ! total =  % completed trajs + contribution form current trajctory
          100*(itraj-1)/c_pars%n_trajs + int(100.*time_new/(c_pars%t_end*c_pars%n_trajs)))

        ! update accumulators
        ! nads
        acc_nads(:,ibin_current+1) = acc_nads(:,ibin_current+1) + lat%n_ads
        ! energy
        energy = total_energy(lat,e_pars)
        acc_en(ibin_current+1) = acc_en(ibin_current+1) + energy
        ! reaction counts
        acc_rcounts(:,:,ibin_current+1) = acc_rcounts(:,:,ibin_current+1) + r%counter

        if (r%n_ads_total>0) then

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
          ! update
          acc_hist(:,:,ibin_current+1) = acc_hist(:,:,ibin_current+1) + hist

          if (c_pars%rdf_period > 0) then
            ! Calculate rdf
            rdf_hist = 0
            call lat%rdf_hist(c_pars, rdf_hist)
            ! update
            acc_rdf(:,:,:,ibin_current+1) = acc_rdf(:,:,:,ibin_current+1) + rdf_hist
          end if

        end if

        ! write into trajectory files if needed
        if (c_pars%output_key /= output_key_avg) then

          ! configuration
          write(outcfg_unit,'(A6,1pe12.3)') "time ",ibin_current*step_bin
          write(outcfg_unit,'(100i10)') lat%n_ads
          call lat%print_ads(outcfg_unit)

          ! energy
          write(outeng_unit,'(1pe12.3,1pe12.3)') ibin_current*step_bin, energy

          ! reaction counts
          write(outcnt_unit,'(1pe12.3,100i20)') ibin_current*step_bin, r%counter

          if (r%n_ads_total>0) then

            ! histogram with cluster sizes

            ! Output formats
            write(n_ads_fmt,'(i10)') r%n_ads_total
            n_ads_fmt = '('//trim(adjustl(n_ads_fmt))//'i8)'
            ! Save cluster size histogram
            write(outhst_unit,'(A6,1pe12.3)') "time ",ibin_current*step_bin
            do species=1,c_pars%n_species
              write(outhst_unit,*) species
              write(outhst_unit,n_ads_fmt) (hist(species,i), i=1,lat%n_ads(species))
              if (n_ads_max(species) < lat%n_ads(species)) n_ads_max(species) = lat%n_ads(species)
            end do

            ! rdf

            if (c_pars%rdf_period > 0) then
              write(outrdf_unit,'(A6,1pe12.3)') "time ",ibin_current*step_bin
              do species1=1,c_pars%n_species
              do species2=1,c_pars%n_species
                write(outrdf_unit,*) species1, species2
                write(outrdf_unit,rdf_fmt) rdf_hist(species1,species2,:)
              end do
              end do
            end if

          else
            ! No histogram output when no adsorbates
            write(outhst_unit,'(A6,1pe12.3)') "time ",ibin_current*step_bin
            write(outrdf_unit,'(A6,1pe12.3)') "time ",ibin_current*step_bin

          end if

        end if

      end do

      ! Reset reaction counters
      r%counter = 0

      ibin = ibin_new
      time = time_new ! time shift

    end do ! over time
!-------------------------------------------------------------

    if (c_pars%output_key /= outcfg_unit) then
      close(outcfg_unit)
      close(outeng_unit)
      close(outhst_unit)
      if (c_pars%rdf_period > 0) close(outrdf_unit)
      close(outcnt_unit)
    end if

  end do ! over trajectories

  ! average accumulated values and save the result
  n_trajs_inv = 1.0_dp/c_pars%n_trajs
  ! n_ads
  call open_for_write(outcfg_unit,trim(c_pars%file_name_base)//'-nads.avg')
  write(outcfg_unit,'(100A12)') 'time(s)', c_pars%ads_names
  do i=1, c_pars%n_bins+1
    write(outcfg_unit,'(100(1pe12.3))') (i-1)*step_bin, acc_nads(:,i)*n_trajs_inv
  end do
  close(outcfg_unit)

  ! energy
  call open_for_write(outeng_unit,trim(c_pars%file_name_base)//'-en.avg')
  write(outeng_unit,'(2A12)') 'time(s)', ' energy(eV)'
  do i=1, c_pars%n_bins+1
    write(outeng_unit,'(1pe12.3,1pe12.3)') (i-1)*step_bin, acc_en(i)*n_trajs_inv
  end do
  close(outeng_unit)
  close(outcfg_unit)

  ! reaction counts
  call open_for_write(outcnt_unit,trim(c_pars%file_name_base)//'-counts.avg')
  write(outcnt_unit,'(A20,100A20)') 'time(s)', &
      ((' '//trim(c_pars%ads_names(j))//'_'//trim(reaction_names(k)), k=1,n_reaction_types), j=1,c_pars%n_species)
  do i=1, c_pars%n_bins+1
    write(outcnt_unit,'(100(1pe20.3))') (i-1)*step_bin, acc_rcounts(:,:,i)*n_trajs_inv
  end do
  close(outcnt_unit)

  ! cluster size histogram
  call open_for_write(outhst_unit,trim(c_pars%file_name_base)//'-hist.avg')
  ! Output formats
  write(n_ads_fmt,'(i10)') c_pars%n_cols*c_pars%n_rows
  n_ads_fmt = '('//trim(adjustl(n_ads_fmt))//'(1pe12.3))'
  do i=1, c_pars%n_bins+1
    ! Save cluster size histogram
    write(outhst_unit,'(A6,1pe12.3)') "time ", (i-1)*step_bin
    do species=1,c_pars%n_species
      write(outhst_unit,*) species
      write(outhst_unit,n_ads_fmt) (acc_hist(species,j,i)*n_trajs_inv, j=1,n_ads_max(species))
    end do
  end do
  close(outhst_unit)

  ! rdf
  if (c_pars%rdf_period > 0) then
    call open_for_write(outrdf_unit,trim(c_pars%file_name_base)//'-rdf.avg')
    ! Output formats
    write(rdf_fmt,'(i10)') c_pars%rdf_n_bins
    rdf_fmt = '('//trim(adjustl(rdf_fmt))//'(1pe12.3))'
    do i=1, c_pars%n_bins+1
      ! Save rdf
      write(outrdf_unit,'(A6,1pe12.3)') "time ",(i-1)*step_bin
      do species1=1,c_pars%n_species
      do species2=1,c_pars%n_species
        write(outrdf_unit,*) species1, species2
        write(outrdf_unit,rdf_fmt) acc_rdf(species1,species2,:,i)*n_trajs_inv
      end do
      end do
    end do
    close(outrdf_unit)
  end if

  ! compress trajectory files
  if (c_pars%output_key == output_key_gz) then
    call compress(c_pars%file_name_base, "confs")
    call compress(c_pars%file_name_base, "en")
    call compress(c_pars%file_name_base, "counts")
    call compress(c_pars%file_name_base, "hist")
    if (c_pars%rdf_period > 0) call compress(c_pars%file_name_base, "rdf")
  end if


end subroutine Bortz_Kalos_Lebowitz


end module kmc
