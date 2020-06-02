program kmc_tian

use mc_lat_class
use control_parameters_class
use energy_parameters_class
use rates_class
use mmc
!use utilities

implicit none

type(mc_lat) :: lattice     ! Declare a variable of type mc_lat.
type(control_parameters) :: control_pars
type( energy_parameters) :: energy_pars
type( rates_type)        :: rates

character(len=max_string_length) file_name_base


! Take a file name base
select case (command_argument_count())
    case(1)
      call get_command_argument(1,file_name_base)
    case default
      stop "Wrong number of arguments"
end select

! Read simulation parameters
control_pars = control_parameters_init(file_name_base)
!print*, control_pars
energy_pars = energy_parameters_init(control_pars)
!print*, energy_pars
rates = rates_init(control_pars)
!call rates%print_r_hop(control_pars)

!   initialize lattice
lattice = mc_lat_init(control_pars)

  call lattice%print_st
  call lattice%print_ocs
  call lattice%print_ads

!
!! configuration input/output format
!write(cfg_fmt,'(i6)') nlat
!cfg_fmt = '('//trim(adjustl(cfg_fmt))//'i8)'
!write(ads_fmt,'(i6)') nads
!ads_fmt = '('//trim(adjustl(ads_fmt))//'i8)'
!
!if (hist_period > 0) call open_for_write(outhst_unit,trim(fname)//'.hist')
!

!! List of opposite directions
!do m=1,nnn
!    nn_opps(m) = modulo(m-1+nnn/2,nnn)+1
!end do
!
!! List of additional nn directions to scan after hop
!do m=1,nnn
!do n=1,nnn/2
!    nn_new(m,n) = modulo(m+n-nnn/2,nnn) + 1
!end do
!end do
!
!
!hist = 0
!hist_counter = 0
!
!write(outcfg_unit,*) nlat, nads
!write(outcfg_unit,cfg_fmt) transpose(occupations)
!
!
select case (control_pars%algorithm)
!
  case ('mmc')

    call metropolis(lattice, control_pars, energy_pars)


  case ('kmc')
!
!        ! Read in the rates and
!        ! construct an array of free-particle hopping rates
!        ! with account for detailed balance
!
!        r_hop = -1.0d0
!
!        call open_for_read(inp_unit,rate_file)
!
!
!        ios = 0
!        do while (ios == 0)
!
!           read(inp_unit, '(A)', iostat=ios) buffer
!
!           if (ios == 0) then
!
!                ! Find the first instance of whitespace.  Split label and data.
!                pos1 = scan(buffer, ' ')
!                label = buffer(1:pos1)
!                buffer = buffer(pos1+1:)
!
!                select case (trim(label))
!
!                case('Arrhenius')
!
!                    read(buffer,*) st_tag, rate_par1, rate_par2
!                    rate_par2 = rate_par2 * eV2K
!
!                    select case (st_tag)
!
!                    case('tt')
!                        r_hop(terrace_site,terrace_site) = &
!                            arrhenius(temperature, rate_par1, rate_par2)
!                    case('ss')
!                        r_hop(step_site,step_site) = &
!                            arrhenius(temperature, rate_par1, rate_par2)
!                    case('ts')
!                        r_hop(terrace_site,step_site) = &
!                            arrhenius(temperature, rate_par1, rate_par2)
!                        r_hop(step_site,terrace_site) = r_hop(terrace_site,step_site)
!                    case('cc')
!                        r_hop(corner_site,corner_site) = &
!                            arrhenius(temperature, rate_par1, rate_par2)
!                    case('tc')
!                        r_hop(terrace_site,corner_site) = &
!                            arrhenius(temperature, rate_par1, rate_par2)
!                        r_hop(corner_site,terrace_site) = r_hop(terrace_site,corner_site)
!                    case('sc')
!                        r_hop(step_site,corner_site) = &
!                            arrhenius(temperature, rate_par1, rate_par2)
!                        r_hop(corner_site,step_site) = r_hop(step_site,corner_site)
!
!                    case default
!                        print*, 'Error in the rate file: unknown site-type tag: ', st_tag
!                        stop
!                    end select
!
!                case('extArrhenius')
!                    stop 'Extended Arrhenius is not yet implemented.'
!
!                case default
!
!                    if (label(1:1) /= '!') stop 'Error in the rate file: unknown rate form.'
!
!                end select
!
!            end if
!
!        end do
!
!        close(inp_unit)
!
!        if (any(r_hop < 0)) stop 'Error in the rate file: not all the rates are defined.'
!
!!        print'(4i6)', transpose(site_type)
!        print*,"Free-particle hopping rate is: "
!        print'(3e16.6)', r_hop
!        print*
!!        stop 44
!
!        allocate(rates(nads,nnn))
!
!        ! time binning for distributions
!        if (n_bins > 0) step_bin = t_end/n_bins
!
!        do itraj=1, ntrajs
!
!            print*, 'Running trajectory no.',itraj
!            print*
!            ibin = 1
!            kmc_nsteps = 0
!
!            ! initialize random number generator
!            call random_seed(size=nseed)
!            call random_seed(put=itraj*seed)
!
!            write(buffer,'(i6.6)') itraj
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
!
!
!            ! Construct rate array
!            do i=1,nads
!
!                i_old = ads_list(i,1)
!                j_old = ads_list(i,2)
!                st_old = site_type(i_old,j_old)
!
!                energy_acc_old = ads_energy(st_old)
!
!                ! Calculate the interaction energy for particle i in its old positions
!                ! Loop over the neighbors of particle i
!                do m=1,nnn
!
!                    nn_pos(m,1) = modulo(i_old + nn_list(m,1)-1,nlat) + 1
!                    nn_pos(m,2) = modulo(j_old + nn_list(m,2)-1,nlat) + 1
!
!                    i_new = nn_pos(m,1)
!                    j_new = nn_pos(m,2)
!                    st_new  = site_type(i_new,j_new)
!
!                    if (occupations(i_new, j_new) > 0)&
!                        energy_acc_old = energy_acc_old + int_energy(st_old, st_new)
!                end do
!
!                ! Calculate the interaction energy for particle i in its NEW positions
!
!                ! Loop over possible new positions of particle i
!                do m=1,nnn
!
!                    i_new = nn_pos(m,1)
!                    j_new = nn_pos(m,2)
!                    st_new  = site_type(i_new,j_new)
!
!                    if (occupations(i_new, j_new) == 0) then
!
!                        ! Excluding self-counting
!                        energy_acc_new = ads_energy(st_new) - int_energy(st_old,st_new)
!
!                        ! Loop over the neighbors of particle i in its possible NEW positions
!                        do m2=1,nnn
!
!                            i_new2 = modulo(i_new + nn_list(m2,1)-1,nlat) + 1
!                            j_new2 = modulo(j_new + nn_list(m2,2)-1,nlat) + 1
!                            st_new2= site_type(i_new2,j_new2)
!
!            !                print*, 'i, i_old, j_old, st_old :',i, i_old, j_old, st_old
!            !                print*, 'm, i_new, j_new, st_new :',m, i_new, j_new, st_new
!            !                print*, 'm2,i_new2,j_new2,st_new2:',m2, i_new2,j_new2, st_new2
!            !                pause
!
!                            if (occupations(i_new2,j_new2) > 0) &
!                                energy_acc_new = energy_acc_new + int_energy(st_new,st_new2)
!                        end do
!
!                        ! Apply detailed balance when
!                        ! total energy in the old position < total energy in the new position
!                        if (energy_acc_old < energy_acc_new) then
!                            rates(i,m) = r_hop(st_old,st_new)&
!                                *exp( -beta*(energy_acc_new - energy_acc_old) )
!                        else
!                            rates(i,m) = r_hop(st_old,st_new)
!                        end if
!
!                    else ! if the neighbor site is occupied
!                        rates(i,m) = 0.0d0
!                    end if
!
!!                            write(*,cfg_fmt) transpose(occupations)
!!                            print*
!!                            print*, 'i, i_old, j_old, st_old:',i, i_old, j_old, st_old
!!                            print*, 'm, i_new, j_new, st_new:',m, i_new, j_new, st_new
!!                            print'(A,f18.2)', 'rates(i,m):',rates(i,m)
!!                            print'(A,f18.4)', 'Eacc_old  = ',energy_acc_old/eV2K
!!                            print'(A,f18.4)', 'E_acc_new = ',energy_acc_new/eV2K
!!                            pause
!                end do
!
!            end do
!
!            ! start time propagation
!            time = big_bang
!            do while (time<t_end)
!
!                total_rate = sum(rates)  ! calculate total rate
!                u = ran1()*total_rate    ! random number to select a process
!
!                rate_acc = 0.d0
!                extloop: do i_ads=1,nads    ! determine reaction channel
!                         do i_nn=1,nnn
!                            rate_acc = rate_acc + rates(i_ads,i_nn)
!                            if (u < rate_acc) exit extloop
!                        end do
!                end do extloop
!
!                delta_t = -log(ran1())/total_rate   ! when does a hop occur?
!                time_new = time + delta_t
!                kmc_nsteps = kmc_nsteps + 1
!
!                if (time_new > t_end) time_new = t_end
!
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
!                time = time_new ! time shift
!
!                ! Update rate constants
!                ! Particle i_ads is going to hop in direction i_nn
!
!                change_list = 0
!                k_change = 1
!                ! Put the hopping particle into the list
!                change_list(k_change) = i_ads
!
!                ! scan over old neighbors
!                do m=1,nnn
!                    ! (i,j) is a position of neighbor m
!                    i = modulo(ads_list(i_ads,1) + nn_list(m,1)-1,nlat) + 1
!                    j = modulo(ads_list(i_ads,2) + nn_list(m,2)-1,nlat) + 1
!
!                    if (occupations(i,j) > 0) then
!                        k_change = k_change + 1
!                        change_list(k_change) = occupations(i,j)
!                    end if
!
!                end do
!
!                ! a new position of particle i_ads after a hop to a neighbor i_nn
!                i_new = modulo(ads_list(i_ads,1) + nn_list(i_nn,1)-1,nlat) + 1
!                j_new = modulo(ads_list(i_ads,2) + nn_list(i_nn,2)-1,nlat) + 1
!                ! scan over additional new neighbors
!                do m=1,nnn/2
!
!                    ! (i,j) is a position of neighbor m
!                    i = modulo(i_new + nn_list(nn_new(i_nn,m),1)-1,nlat) + 1
!                    j = modulo(j_new + nn_list(nn_new(i_nn,m),2)-1,nlat) + 1
!
!                    if (occupations(i,j) > 0) then
!                        k_change = k_change + 1
!                        change_list(k_change) = occupations(i,j)
!                    end if
!
!                end do
!
!!                print*,i_ads, i_nn
!!                print(cfg_fmt), transpose(occupations)
!!                print*
!!                print'(6e16.4)', (rates(i,:)/r_hop,i=1,nads)
!                !print*,k
!!                print*,change_list
!                !stop 15
!
!                ! Update occupations and ads. list
!                occupations(ads_list(i_ads,1),ads_list(i_ads,2)) = 0
!                occupations(i_new,j_new) = i_ads
!                ads_list(i_ads,:) = (/i_new,j_new/)
!
!                ! Update rate array
!                ! see building up the rate array above
!                ! with i replaced by change_list(i)
!
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
!            end do ! over time
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
!        end do ! over trajectories
!
!        deallocate(rates)
!
!    case default
!        stop 'Error: mc algorithm not defined'
end select
!
!
!close(outcfg_unit)
!
!if (hist_period > 0) close(outhst_unit)
!
!call system_clock(icount2)
!print*
!print*, 'Number of clock ticks = ', icount2-icount1
!
call open_for_write(outcfg_unit,trim(control_pars%file_name_base)//'.out')

write(outcfg_unit,'(A10,A10,A15)') "# rows","# cols","step_period"
write(outcfg_unit,'(3i10)') lattice%n_rows, lattice%n_cols, &
                            control_pars%step_period
write(outcfg_unit,'(100A10)') adjustr(control_pars%ads_names)
write(outcfg_unit,'(100i10)') lattice%n_ads
write(outcfg_unit,'(5A10)') "#","row","col","ads_site", "species"
call lattice%print_ads(outcfg_unit)


close(outcfg_unit)

!deallocate(cluster_label, cluster_sizes, hist, change_list, nn_new)
!deallocate(ads_list,nn_pos,nn_list,nn_opps,temp1D,site_type,occupations)


end program
