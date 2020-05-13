program mmc

! load module with random number and absorption rate functions
use utilities
use open_file

implicit none

! Assignment for the surface site types
integer(1), parameter :: terrace_site = 1
integer(1), parameter :: step_site    = 2
integer(1), parameter :: corner_site  = 3

integer, parameter :: inp_unit  = 1
integer, parameter :: outcfg_unit  = 10
integer, parameter :: outeng_unit  = 11
integer, parameter :: outhst_unit  = 12

real(8), parameter :: eV2K = 11604.52
real(8), parameter :: big_bang = 0.0d0 ! in fortnights
real(8), parameter :: pi = 4.0d0*atan(1.0d0)

character(len=3) :: algorithm   ! MC algorithm to use (kmc or mmc)
integer :: nlat         ! size of 2D lattice (nlat x nlat)
integer :: step_period  ! = step_density^-1, 0 means no steps
real(8) :: temperature  ! temperature in K
real(8) :: coverage     ! coverage in ML per 2D MC lattice (nlat x nlat)
real(8) :: eps          ! O-O interaction energy in eV
integer :: save_period  ! period for conf. output
integer :: hist_period  ! period for histogram  calculation
integer :: nsteps   ! number of Metropolis MC steps
integer :: ntrajs	! number of kmc trajectories (ignored for mmc)
real(8) :: t_end     ! kmc simulation time (s)
integer :: n_bins  	! number of time intervals in kmc histogram
character(len=120) :: rate_file='' ! name of the file with rate parameters
character(len=120) :: energy_file=''! name of the file with adsorption and interaction energy
character(len=120) :: cfg_fname='' ! name of the file with initial configuration

real(8) :: energy, total_energy, arrhenius

integer :: nseed = 8
integer, parameter :: seed(8) = (/1,6,3,5,7,3,3,7/)
integer :: icount1, icount2

integer :: nads, nnn, nn_counter, nlat_old, nads_old
integer :: i, j, k, m, m2, n, ihop, istep, itraj, kk, k_change
integer :: i_old, j_old, i_new, j_new, i_new2, j_new2
integer :: st_old, st_new, st_new2

real(8), dimension(2) :: b1, b2
integer, dimension(:,:), allocatable   :: occupations, ads_list
integer(1), dimension(:,:), allocatable   :: site_type
integer, dimension(:,:), allocatable   :: nn_list, nn_pos, nn_new
integer, dimension(:),   allocatable   :: nn_opps, temp1D, change_list

integer, dimension(:),   allocatable   :: cluster_sizes, hist
integer, dimension(:,:), allocatable   :: cluster_label

real(8) :: time, delta_t, time_new, step_bin
real(8) :: rate_acc, u, total_rate
integer :: i_nn, i_ads, ibin_new, ibin, kmc_nsteps
real(8), dimension(:,:),   allocatable :: rates

real(8) :: energy_old, delta_E, beta
real(8) :: energy_acc_old, energy_acc_new

integer :: ios, pos1, largest_label, hist_counter
character(len=120) buffer, label, fname

real(8), dimension(3) :: ads_energy
real(8), dimension(3,3) :: int_energy

real(8), dimension(3,3) :: r_hop
character(len=2) :: st_tag
real(8) :: rate_par1, rate_par2, rtemp

character(len=10) cfg_fmt
character(len=13) ads_fmt

! Read in simulation parameters

select case (iargc())

case(1)
    call getarg(1,fname)
case default
    stop "Wrong number of arguments"
end select


call open_for_read(inp_unit, trim(fname)//'.inp' )
ios = 0
do while (ios == 0)

        read(inp_unit, '(A)', iostat=ios) buffer
        if (ios == 0) then

        ! Find the first instance of whitespace.  Split label and data.
            pos1 = scan(buffer, ' ')
            label = buffer(1:pos1)
            buffer = buffer(pos1+1:)

            select case (label)
            case('algorithm')
                read(buffer,*,iostat=ios) algorithm
            case('nlat')
                read(buffer,*,iostat=ios) nlat
            case('step_period')
                read(buffer,*,iostat=ios) step_period
            case('temperature')
                read(buffer,*,iostat=ios) temperature
                beta = 1.0d0/temperature ! thermodynamic temperature
            case('coverage')
                read(buffer,*,iostat=ios) coverage
            case('energy')
                read(buffer,*,iostat=ios) energy_file
                energy_file = trim(energy_file)
!            case('eps')
!                read(buffer,*,iostat=ios) eps
!                eps = eps * eV2K
            case('save_period')
                read(buffer,*,iostat=ios) save_period
            case('ini_conf')
                read(buffer,*,iostat=ios) cfg_fname
            case('mmc_hist_period')
                read(buffer,*,iostat=ios) hist_period
            case('mmc_nsteps')
                read(buffer,*,iostat=ios) nsteps
            case('kmc_ntrajs')
                read(buffer,*,iostat=ios) ntrajs
            case('kmc_time')
                read(buffer,*,iostat=ios) t_end
            case('kmc_nbins')
                read(buffer,*,iostat=ios) n_bins
            case('kmc_rates')
                read(buffer,*,iostat=ios) rate_file
                rate_file = trim(rate_file)
            case default
                if (label(1:1) /= '!')&
                    print *, 'Skipping invalid label at line', label
            end select
        end if

end do ! ios

close(inp_unit)

! Number of adsorbate particles
nads = nint(coverage*nlat*nlat)

! Consistency checks

if (step_period > 0 .and. mod(nlat,step_period) /= 0)&
    stop "Error in the inp file: inconsistent nlat and step_period values"

call open_for_read(inp_unit,energy_file)

ios = 0
do while (ios == 0)

        read(inp_unit, '(A)', iostat=ios) buffer
        if (ios == 0) then

        ! Find the first instance of whitespace.  Split label and data.
            pos1 = scan(buffer, ' ')
            label = buffer(1:pos1)
            buffer = buffer(pos1+1:)

            select case (label)

            case('adsorption')

                read(buffer,*,iostat=ios) &
                    ads_energy(terrace_site),&
                    ads_energy(   step_site),&
                    ads_energy( corner_site)
                if (ios /= 0) stop 'Error in energy file: adsorption'
                ads_energy = ads_energy*eV2K

            case('interaction')
                read(buffer,*,iostat=ios)&
                    int_energy(terrace_site,terrace_site),&
                    int_energy(   step_site,   step_site),&
                    int_energy(terrace_site,   step_site),&
                    int_energy( corner_site, corner_site),&
                    int_energy(terrace_site, corner_site),&
                    int_energy(   step_site, corner_site)
                if (ios /= 0) stop 'Error in energy file: interaction'

                int_energy(  step_site,terrace_site) = int_energy(terrace_site,  step_site)
                int_energy(corner_site,terrace_site) = int_energy(terrace_site,corner_site)
                int_energy(corner_site,   step_site) = int_energy(   step_site,corner_site)

                int_energy = int_energy*eV2K

            case default
                if (label(1:1) /= '!')&
                    print *, 'Skipping invalid label at line', label

            end select

        end if

end do ! ios

close(inp_unit)


! Number of neighbors for the hexagonal structure
nnn = 6

! configuration input/output format
write(cfg_fmt,'(i6)') nlat
cfg_fmt = '('//trim(adjustl(cfg_fmt))//'i8)'
write(ads_fmt,'(i6)') nads
ads_fmt = '('//trim(adjustl(ads_fmt))//'i8)'

call open_for_write(outcfg_unit,trim(fname)//'.confs')
if (hist_period > 0) call open_for_write(outhst_unit,trim(fname)//'.hist')

! allocate memory for arrays
allocate(occupations(nlat,nlat), site_type(nlat,nlat))
allocate(ads_list(nads,2), change_list(2*nnn),nn_new(nnn,nnn/2))
allocate(nn_list(nnn,2), nn_pos(nnn,2), nn_opps(nnn), temp1D(nlat*nlat))
allocate(cluster_label(nlat,nlat), cluster_sizes(nads), hist(nads))

! NN list for the hexagonal structure
!  11    12*   13*   14
!
!     21*   22*   23*   24
!
!        31*   32*   33    34
!
!           41    42    43    44

nn_list(1,:) = (/ 0, 1/)
nn_list(2,:) = (/ 1, 0/)
nn_list(3,:) = (/ 1,-1/)
nn_list(4,:) = (/ 0,-1/)
nn_list(5,:) = (/-1, 0/)
nn_list(6,:) = (/-1, 1/)

! List of opposite directions
do m=1,nnn
    nn_opps(m) = modulo(m-1+nnn/2,nnn)+1
end do

! List of additional nn directions to scan after hop
do m=1,nnn
do n=1,nnn/2
    nn_new(m,n) = modulo(m+n-nnn/2,nnn) + 1
end do
end do

! Hexagonal transformation matrix

b1 = (/         1.0d0,       0.0d0/)
b2 = (/ cos(pi/3.0d0), sin(pi/3.0d0) /)

site_type = terrace_site
! Define where steps are
if (step_period > 0) then
    do j=1,nlat,step_period
        site_type(:,j)   = step_site
        site_type(:,j+1) = corner_site
    end do
end if

!print*
!write(*,cfg_fmt) transpose(site_type)
!stop

if (cfg_fname=='') then

    temp1D = 0
    do i=1,nads
        temp1D(i) = i
    end do
    occupations = reshape(temp1D,(/nlat,nlat/))

else

    call open_for_read(inp_unit,trim(cfg_fname))
    read(inp_unit,*) nlat_old, nads_old
    if (nlat_old /= nlat .OR. nads_old /= nads ) &
        stop 'Error: inconsistent input and configuration files'
    read(inp_unit,cfg_fmt) (occupations(i,:), i=1,nlat)
    close(inp_unit)
    print*
    print*,'    Dear Sir,'
    print*,'I would like to humbly notify you that the initial configuration is read from:'
    print*, trim(cfg_fname)
    print*,'Remaining your humble servant, K.M.C. Code'
    print*

end if

do j=1,nlat
do i=1,nlat
    k = occupations(i,j)
    if (k > 0) then
        ads_list(k,1) = i
        ads_list(k,2) = j
    end if
end do
end do

hist = 0
hist_counter = 0

write(outcfg_unit,*) nlat, nads
write(outcfg_unit,cfg_fmt) transpose(occupations)


select case (algorithm)

    case ('mmc')

        call open_for_write(outeng_unit,trim(fname)//'.en')
        write(outeng_unit,*) 0,&
            total_energy(nlat, nads, nnn, occupations, site_type, &
                         ads_list, nn_list, ads_energy, int_energy)/eV2K

!        write(*,cfg_fmt) transpose(site_type)
!        print*
!        write(*,cfg_fmt) transpose(occupations)
!        print*
!        print*,energy_file
!        write(*,'(3f12.3)') int_energy/eV2K
!        write(*,*) total_energy(nlat, nads, nnn, occupations, site_type, &
!                         ads_list, nn_list, ads_energy, int_energy)/eV2K
!        stop 33

        !loop over mmc steps
        do istep=1, nsteps

            do i=1, nads

                energy_old = energy(i, nlat, nads, nnn, occupations, site_type, &
                                       ads_list, nn_list, ads_energy, int_energy)

                ihop = floor(nnn*ran1()) + 1

                i_new = modulo(ads_list(i,1) + nn_list(ihop,1)-1,nlat) + 1
                j_new = modulo(ads_list(i,2) + nn_list(ihop,2)-1,nlat) + 1

                if (occupations(i_new, j_new) == 0) then

                    i_old = ads_list(i,1)
                    j_old = ads_list(i,2)

                    occupations(i_new,j_new) = i
                    occupations(i_old,j_old) = 0
                    ads_list(i,:) = (/i_new,j_new/)

                    delta_E = energy(i, nlat, nads, nnn, occupations, site_type, &
                                     ads_list, nn_list, ads_energy, int_energy) &
                            - energy_old

                    if (exp(- delta_E/temperature) < ran1()) then

                        occupations(i_old,j_old) = i
                        occupations(i_new,j_new) = 0
                        ads_list(i,:) = (/i_old,j_old/)

                    end if

                end if

            enddo

!print*,istep, hist_period, mod(istep, hist_period)

        if (hist_period > 0 .and. mod(istep, hist_period) == 0) then
            call hoshen_kopelman(cluster_label, largest_label, occupations, ads_list, &
                                                nn_list, nads, nlat, nnn)

            call count_cluster_sizes(cluster_sizes, cluster_label,&
                                                                ads_list, nads, nlat)

            hist_counter = hist_counter + 1
            do i=1,largest_label
                if (cluster_sizes(i) > 0) &
                    hist(cluster_sizes(i)) = hist(cluster_sizes(i)) + 1
            end do

        !print*,largest_label
        !print*,cluster_sizes
        !print*
        !write(*,cfg_fmt) (cluster_label(m,:), m=1,nlat)
        !print*
        !write(*,cfg_fmt) hist
        !

        end if

            if (mod(istep, save_period) == 0) then
                print*, istep
                write(outcfg_unit,cfg_fmt) transpose(occupations)
                write(outeng_unit,*) istep-1,&
                    total_energy(nlat, nads, nnn, occupations, site_type, &
                                 ads_list, nn_list, ads_energy, int_energy)/eV2K

                if (hist_period > 0) then
                    write(outhst_unit,*) hist_counter
                    write(outhst_unit,ads_fmt) hist
                end if
           end if

        enddo

        close(outeng_unit)

    case ('kmc')

        ! Read in the rates and
        ! construct an array of free-particle hopping rates
        ! with account for detailed balance

        r_hop = -1.0d0

        call open_for_read(inp_unit,rate_file)


        ios = 0
        do while (ios == 0)

           read(inp_unit, '(A)', iostat=ios) buffer

           if (ios == 0) then

                ! Find the first instance of whitespace.  Split label and data.
                pos1 = scan(buffer, ' ')
                label = buffer(1:pos1)
                buffer = buffer(pos1+1:)

                select case (trim(label))

                case('Arrhenius')

                    read(buffer,*) st_tag, rate_par1, rate_par2
                    rate_par2 = rate_par2 * eV2K

                    select case (st_tag)

                    case('tt')
                        r_hop(terrace_site,terrace_site) = &
                            arrhenius(temperature, rate_par1, rate_par2)
                    case('ss')
                        r_hop(step_site,step_site) = &
                            arrhenius(temperature, rate_par1, rate_par2)
                    case('ts')
                        r_hop(terrace_site,step_site) = &
                            arrhenius(temperature, rate_par1, rate_par2)
                        r_hop(step_site,terrace_site) = r_hop(terrace_site,step_site)
                    case('cc')
                        r_hop(corner_site,corner_site) = &
                            arrhenius(temperature, rate_par1, rate_par2)
                    case('tc')
                        r_hop(terrace_site,corner_site) = &
                            arrhenius(temperature, rate_par1, rate_par2)
                        r_hop(corner_site,terrace_site) = r_hop(terrace_site,corner_site)
                    case('sc')
                        r_hop(step_site,corner_site) = &
                            arrhenius(temperature, rate_par1, rate_par2)
                        r_hop(corner_site,step_site) = r_hop(step_site,corner_site)

                    case default
                        print*, 'Error in the rate file: unknown site-type tag: ', st_tag
                        stop
                    end select

                case('extArrhenius')
                    stop 'Extended Arrhenius is not yet implemented.'

                case default

                    if (label(1:1) /= '!') stop 'Error in the rate file: unknown rate form.'

                end select

            end if

        end do

        close(inp_unit)

        if (any(r_hop < 0)) stop 'Error in the rate file: not all the rates are defined.'

!        print'(4i6)', transpose(site_type)
        print*,"Free-particle hopping rate is: "
        print'(3e16.6)', r_hop
        print*
!        stop 44

        allocate(rates(nads,nnn))

        ! time binning for distributions
        if (n_bins > 0) step_bin = t_end/n_bins

        do itraj=1, ntrajs

            print*, 'Running trajectory no.',itraj
            print*
            ibin = 1
            kmc_nsteps = 0

            ! initialize random number generator
            call random_seed(size=nseed)
            call random_seed(put=itraj*seed)

            write(buffer,'(i6.6)') itraj
            if (n_bins > 0) then

                call open_for_write(outeng_unit,trim(fname)//trim(buffer)//'.en')
                write(outeng_unit,*) big_bang ,&
                    total_energy(nlat, nads, nnn, occupations, site_type, &
                                 ads_list, nn_list, ads_energy, int_energy)/eV2K

                call open_for_write(outhst_unit,trim(fname)//trim(buffer)//'.csz')
                write(outhst_unit,*) t_end, n_bins

            end if


            ! Construct rate array
            do i=1,nads

                i_old = ads_list(i,1)
                j_old = ads_list(i,2)
                st_old = site_type(i_old,j_old)

                energy_acc_old = ads_energy(st_old)

                ! Calculate the interaction energy for particle i in its old positions
                ! Loop over the neighbors of particle i
                do m=1,nnn

                    nn_pos(m,1) = modulo(i_old + nn_list(m,1)-1,nlat) + 1
                    nn_pos(m,2) = modulo(j_old + nn_list(m,2)-1,nlat) + 1

                    i_new = nn_pos(m,1)
                    j_new = nn_pos(m,2)
                    st_new  = site_type(i_new,j_new)

                    if (occupations(i_new, j_new) > 0)&
                        energy_acc_old = energy_acc_old + int_energy(st_old, st_new)
                end do

                ! Calculate the interaction energy for particle i in its NEW positions

                ! Loop over possible new positions of particle i
                do m=1,nnn

                    i_new = nn_pos(m,1)
                    j_new = nn_pos(m,2)
                    st_new  = site_type(i_new,j_new)

                    if (occupations(i_new, j_new) == 0) then

                        ! Excluding self-counting
                        energy_acc_new = ads_energy(st_new) - int_energy(st_old,st_new)

                        ! Loop over the neighbors of particle i in its possible NEW positions
                        do m2=1,nnn

                            i_new2 = modulo(i_new + nn_list(m2,1)-1,nlat) + 1
                            j_new2 = modulo(j_new + nn_list(m2,2)-1,nlat) + 1
                            st_new2= site_type(i_new2,j_new2)

            !                print*, 'i, i_old, j_old, st_old :',i, i_old, j_old, st_old
            !                print*, 'm, i_new, j_new, st_new :',m, i_new, j_new, st_new
            !                print*, 'm2,i_new2,j_new2,st_new2:',m2, i_new2,j_new2, st_new2
            !                pause

                            if (occupations(i_new2,j_new2) > 0) &
                                energy_acc_new = energy_acc_new + int_energy(st_new,st_new2)
                        end do

                        ! Apply detailed balance when
                        ! total energy in the old position < total energy in the new position
                        if (energy_acc_old < energy_acc_new) then
                            rates(i,m) = r_hop(st_old,st_new)&
                                *exp( -beta*(energy_acc_new - energy_acc_old) )
                        else
                            rates(i,m) = r_hop(st_old,st_new)
                        end if

                    else ! if the neighbor site is occupied
                        rates(i,m) = 0.0d0
                    end if

!                            write(*,cfg_fmt) transpose(occupations)
!                            print*
!                            print*, 'i, i_old, j_old, st_old:',i, i_old, j_old, st_old
!                            print*, 'm, i_new, j_new, st_new:',m, i_new, j_new, st_new
!                            print'(A,f18.2)', 'rates(i,m):',rates(i,m)
!                            print'(A,f18.4)', 'Eacc_old  = ',energy_acc_old/eV2K
!                            print'(A,f18.4)', 'E_acc_new = ',energy_acc_new/eV2K
!                            pause
                end do

            end do

            ! start time propagation
            time = big_bang
            do while (time<t_end)

                total_rate = sum(rates)  ! calculate total rate
                u = ran1()*total_rate    ! random number to select a process

                rate_acc = 0.d0
                extloop: do i_ads=1,nads    ! determine reaction channel
                         do i_nn=1,nnn
                            rate_acc = rate_acc + rates(i_ads,i_nn)
                            if (u < rate_acc) exit extloop
                        end do
                end do extloop

                delta_t = -log(ran1())/total_rate   ! when does a hop occur?
                time_new = time + delta_t
                kmc_nsteps = kmc_nsteps + 1

                if (time_new > t_end) time_new = t_end

                if(n_bins > 0) then


                    ibin_new = int(time_new/step_bin) + 1

                    if (ibin_new - ibin > 0 ) then
                        call hoshen_kopelman(cluster_label, largest_label, occupations, ads_list, &
                                                nn_list, nads, nlat, nnn)
                        call count_cluster_sizes(cluster_sizes, cluster_label,&
                                                                ads_list, nads, nlat)
                        hist = 0
                        do i=1,largest_label
                            if (cluster_sizes(i) > 0) &
                                hist(cluster_sizes(i)) = hist(cluster_sizes(i)) + 1
                        end do

                        write(outhst_unit,*) ibin
                        write(outhst_unit,ads_fmt) hist

!                        write(6,cfg_fmt) (occupations(i,:), i=1,nlat)
                        write(outeng_unit,*) time,&
                            total_energy(nlat, nads, nnn, occupations, site_type, &
                                         ads_list, nn_list, ads_energy, int_energy)/eV2K

                    end if

                    ibin = ibin_new

                end if

                time = time_new ! time shift

                ! Update rate constants
                ! Particle i_ads is going to hop in direction i_nn

                change_list = 0
                k_change = 1
                ! Put the hopping particle into the list
                change_list(k_change) = i_ads

                ! scan over old neighbors
                do m=1,nnn
                    ! (i,j) is a position of neighbor m
                    i = modulo(ads_list(i_ads,1) + nn_list(m,1)-1,nlat) + 1
                    j = modulo(ads_list(i_ads,2) + nn_list(m,2)-1,nlat) + 1

                    if (occupations(i,j) > 0) then
                        k_change = k_change + 1
                        change_list(k_change) = occupations(i,j)
                    end if

                end do

                ! a new position of particle i_ads after a hop to a neighbor i_nn
                i_new = modulo(ads_list(i_ads,1) + nn_list(i_nn,1)-1,nlat) + 1
                j_new = modulo(ads_list(i_ads,2) + nn_list(i_nn,2)-1,nlat) + 1
                ! scan over additional new neighbors
                do m=1,nnn/2

                    ! (i,j) is a position of neighbor m
                    i = modulo(i_new + nn_list(nn_new(i_nn,m),1)-1,nlat) + 1
                    j = modulo(j_new + nn_list(nn_new(i_nn,m),2)-1,nlat) + 1

                    if (occupations(i,j) > 0) then
                        k_change = k_change + 1
                        change_list(k_change) = occupations(i,j)
                    end if

                end do

!                print*,i_ads, i_nn
!                print(cfg_fmt), transpose(occupations)
!                print*
!                print'(6e16.4)', (rates(i,:)/r_hop,i=1,nads)
                !print*,k
!                print*,change_list
                !stop 15

                ! Update occupations and ads. list
                occupations(ads_list(i_ads,1),ads_list(i_ads,2)) = 0
                occupations(i_new,j_new) = i_ads
                ads_list(i_ads,:) = (/i_new,j_new/)

                ! Update rate array
                ! see building up the rate array above
                ! with i replaced by change_list(i)

                do kk=1,k_change

                    i = change_list(kk)
                    i_old = ads_list(i,1)
                    j_old = ads_list(i,2)
                    st_old = site_type(i_old,j_old)

                    energy_acc_old = ads_energy(st_old)

                    ! Calculate the interaction energy for particle i in its old positions
                    ! Loop over the neighbors of particle i
                    do m=1,nnn

                        nn_pos(m,1) = modulo(i_old + nn_list(m,1)-1,nlat) + 1
                        nn_pos(m,2) = modulo(j_old + nn_list(m,2)-1,nlat) + 1

                        i_new = nn_pos(m,1)
                        j_new = nn_pos(m,2)
                        st_new  = site_type(i_new,j_new)

                        if (occupations(i_new, j_new) > 0)&
                            energy_acc_old = energy_acc_old + int_energy(st_old, st_new)
                    end do

                    ! Calculate the interaction energy for particle i in its NEW positions

                    ! Loop over possible new positions of particle i
                    do m=1,nnn

                        i_new = nn_pos(m,1)
                        j_new = nn_pos(m,2)
                        st_new  = site_type(i_new,j_new)

                        if (occupations(i_new, j_new) == 0) then

                            ! Excluding self-counting
                            energy_acc_new = ads_energy(st_new) - int_energy(st_old,st_new)

                            ! Loop over the neighbors of particle i in its possible NEW positions
                            do m2=1,nnn

                                i_new2 = modulo(i_new + nn_list(m2,1)-1,nlat) + 1
                                j_new2 = modulo(j_new + nn_list(m2,2)-1,nlat) + 1
                                st_new2= site_type(i_new2,j_new2)

                !                print*, 'i, i_old, j_old, st_old :',i, i_old, j_old, st_old
                !                print*, 'm, i_new, j_new, st_new :',m, i_new, j_new, st_new
                !                print*, 'm2,i_new2,j_new2,st_new2:',m2, i_new2,j_new2, st_new2
                !                pause

                                if (occupations(i_new2,j_new2) > 0) &
                                    energy_acc_new = energy_acc_new + int_energy(st_new,st_new2)
                            end do

                            ! Apply detailed balance when
                            ! total energy in the old position < total energy in the new position
                            if (energy_acc_old < energy_acc_new) then
                                rates(i,m) = r_hop(st_old,st_new)&
                                    *exp( -beta*(energy_acc_new - energy_acc_old) )
                            else
                                rates(i,m) = r_hop(st_old,st_new)
                            end if

                        else ! if the neighbor site is occupied
                            rates(i,m) = 0.0d0
                        end if

!                        write(*,cfg_fmt) transpose(occupations)
!                        print*
!                        print*, 'i, i_old, j_old, st_old:',i, i_old, j_old, st_old
!                        print*, 'm, i_new, j_new, st_new:',m, i_new, j_new, st_new
!                        print'(A,f18.2)', 'rates(i,m):',rates(i,m)
!                        print'(A,f18.4)', 'Eacc_old  = ',energy_acc_old/eV2K
!                        print'(A,f18.4)', 'E_acc_new = ',energy_acc_new/eV2K
!                        pause
                    end do

                end do

            end do ! over time
!-------------------------------------------------------------

            if (n_bins > 0) then
                close(outeng_unit)
                close(outhst_unit)
            end if

            print*,'Number of kmc-steps is ', kmc_nsteps
            print*

        end do ! over trajectories

        deallocate(rates)

    case default
        stop 'Error: mc algorithm not defined'
end select


close(outcfg_unit)

if (hist_period > 0) close(outhst_unit)

call system_clock(icount2)
print*
print*, 'Number of clock ticks = ', icount2-icount1

call open_for_write(outcfg_unit,trim(fname)//'.out')

write(outcfg_unit,*) nlat, nads
write(outcfg_unit,cfg_fmt) transpose(occupations)

close(outcfg_unit)

deallocate(cluster_label, cluster_sizes, hist, change_list, nn_new)
deallocate(ads_list,nn_pos,nn_list,nn_opps,temp1D,site_type,occupations)


end program

real(8) function energy(inx, nlat, nads, nnn, occupations, site_type, &
                        ads_list, nn_list, ads_energy, int_energy)

integer, intent(in) :: inx, nlat, nads, nnn
integer, dimension(nlat,nlat), intent(in) :: occupations
integer(1), dimension(nlat,nlat), intent(in) :: site_type
integer, dimension(nads,2), intent(in) :: ads_list
integer, dimension(nnn,2), intent(in) :: nn_list
real(8), dimension(3), intent(in) :: ads_energy
real(8), dimension(3,3), intent(in) :: int_energy

integer :: i, ic, jc, iads, jads
real(8) :: energy_acc

        iads = ads_list(inx,1)
        jads = ads_list(inx,2)

        energy_acc = ads_energy(site_type(iads,jads))

        do i=1, nnn ! Sum up the int. energy over all nearest neighbors

            ic = modulo(iads+nn_list(i,1)-1,nlat) + 1
            jc = modulo(jads+nn_list(i,2)-1,nlat) + 1
            if (occupations(ic,jc) > 0) &
                energy_acc = energy_acc +&
                        int_energy(site_type(iads,jads),site_type(ic,jc))

        end do

    energy = energy_acc

end function

real(8) function total_energy(nlat, nads, nnn, occupations, site_type, &
                              ads_list, nn_list, ads_energy, int_energy)

integer, intent(in) :: nlat, nads, nnn
integer, dimension(nlat,nlat), intent(in) :: occupations
integer(1), dimension(nlat,nlat), intent(in) :: site_type
integer, dimension(nads,2), intent(in) :: ads_list
integer, dimension(nnn,2), intent(in) :: nn_list
real(8), dimension(3), intent(in) :: ads_energy
real(8), dimension(3,3), intent(in) :: int_energy

integer :: i, ic, jc, iads, jads
real(8) :: energy_acc

    energy_acc = 0.0d0
    do inx=1,nads

        iads = ads_list(inx,1)
        jads = ads_list(inx,2)

        energy_acc = energy_acc + ads_energy(site_type(iads,jads))

        do i=1, nnn/2 ! count interaction with neighbors once

            ic = modulo(iads+nn_list(i,1)-1,nlat) + 1
            jc = modulo(jads+nn_list(i,2)-1,nlat) + 1
            if (occupations(ic,jc) > 0) &
                energy_acc = energy_acc +&
                        int_energy(site_type(iads,jads),site_type(ic,jc))

        end do

    end do

    total_energy = energy_acc

end function

integer function lfind(x, labels, nads)

integer, intent(inout) :: x
integer, dimension(nads), intent(inout) :: labels

integer :: y, z

    y = x
    do while (labels(y) /= y)
        y = labels(y)
    end do

!    do while (labels(x) /= x)
!        z = labels(x)
!        labels(x) = y
!        x = z
!    end do
!
    lfind = y

end function

subroutine lunion(x, y, labels, nads)

integer, intent(in) :: x, y, nads
integer, dimension(nads), intent(inout) :: labels

    if (x > y) then
        labels( lfind(x, labels, nads) ) = lfind(y, labels, nads)
    else if ( x < y) then
        labels( lfind(y, labels, nads) ) = lfind(x, labels, nads)
    endif

end subroutine

subroutine hoshen_kopelman(cluster_label, largest_label, occupations, ads_list, nn_list, nads, nlat, nnn)

integer, dimension(nlat,nlat), intent(in) :: occupations
integer, intent(in) :: nads, nlat, nnn
integer, dimension(nads,2), intent(in) :: ads_list
integer, dimension(nnn,2), intent(in) :: nn_list
integer, dimension(nlat,nlat), intent(out) :: cluster_label
integer, intent(out) :: largest_label

integer :: i, j, m, nnn2
integer :: itemp, ip, jp
integer :: lfind

integer, dimension(nads) :: labels
integer, dimension(nnn/2) :: scanned_nn_occs, i_nn, j_nn

nnn2 = nnn/2
scanned_nn_occs = 0
cluster_label   = 0
do i=1,nads
    labels(i) = i
end do

largest_label = 0

!print'(4i4)', transpose(occupations)

do i=1,nlat
do j=1,nlat


    if (occupations(i,j) > 0) then

        do m=1,nnn2

            i_nn(m) = modulo(i + nn_list(nnn2+m,1)-1,nlat) + 1
            j_nn(m) = modulo(j + nn_list(nnn2+m,2)-1,nlat) + 1

            if (occupations(i_nn(m), j_nn(m)) > 0) then
                scanned_nn_occs(m) = 1
            else
                scanned_nn_occs(m) = 0
            end if

            if (i==1) scanned_nn_occs(2:3) = 0
            if (j==1) scanned_nn_occs(1) = 0
            if (j==nlat) scanned_nn_occs(3) = 0

        end do

!print*,scanned_nn_occs

        select case (sum(scanned_nn_occs))

        case (0)

!print*
!print*,"Case 0"
            largest_label = largest_label + 1
            cluster_label(i,j) = largest_label


        case (1)

!print*
!print*,"Case 1"
            do m=1,nnn2
                if (scanned_nn_occs(m) == 1) then
                    cluster_label(i,j) = &
                        lfind(cluster_label(i_nn(m),j_nn(m)), labels, nads)
                end if
            end do

       case (2)

!print*
!print*,"Case 2"

            do m=1,nnn2-1
                itemp = cluster_label(i_nn(m),j_nn(m))
            do n=m+1,nnn2
                if (scanned_nn_occs(m) == 1 .and. scanned_nn_occs(n) == 1) then

                    call lunion(itemp, cluster_label(i_nn(n),j_nn(n)),&
                                                        labels, nads )
                    cluster_label(i,j) = lfind(itemp, labels, nads )

                end if
            end do
            end do

        case (3)

!print*
!print*,"Case 3"
            itemp = cluster_label(i_nn(1),j_nn(1))
            call lunion(itemp, cluster_label(i_nn(2),j_nn(2)),&
                                                        labels, nads )
            call lunion(itemp, cluster_label(i_nn(3),j_nn(3)),&
                                                        labels, nads )
            cluster_label(i,j) = lfind(itemp, labels, nads )

        case default

            stop 'Hoshen-Kopelman: too many neighbors'

        end select

    end if

!print*,i,j, occupations(i,j)
!print*
!print'(4i4)', transpose(cluster_label)
!print*

end do
end do

!print*
!write(*,*) occupations
!print*
!write(*,*) cluster_label
!print*
!do i=1,largest_label
!    print('(i3,a5,i3)'), i, ' --> ' , labels(i)
!end do

! Apply PBC

do i=1,nlat
    ip = cluster_label(i,1)
    if (ip > 0) then
        do m=1,nnn2
            jp = cluster_label( &
                modulo(i + nn_list(nnn2+m,1)-1,nlat) + 1,&
                modulo(1 + nn_list(nnn2+m,2)-1,nlat) + 1  )
            if (jp > 0) call lunion(ip,jp,labels,nads)
        end do
    end if
end do

do i=2,nlat-1
    ip = cluster_label(1,i)
    if (ip > 0) then
        do m=1,nnn2
            jp = cluster_label( &
                modulo(1 + nn_list(nnn2+m,1)-1,nlat) + 1,&
                modulo(i + nn_list(nnn2+m,2)-1,nlat) + 1  )
            if (jp > 0) call lunion(ip,jp,labels,nads)
        end do
    end if
end do

!print*
!do i=1,largest_label
!    print('(i3,a5,i3)'), i, ' --> ' , labels(i)
!end do

! Going down to roots

do i=1,nlat
do j=1,nlat
    if (cluster_label(i,j) > 0) &
     cluster_label(i,j) = lfind(cluster_label(i,j), labels, nads)
end do
end do

end subroutine


subroutine count_cluster_sizes(cluster_sizes, cluster_label, ads_list, nads, nlat)

integer, dimension(nads), intent(out) :: cluster_sizes
integer, intent(in) :: nads, nlat
integer, dimension(nads,2), intent(in) :: ads_list
integer, dimension(nlat,nlat), intent(in) :: cluster_label


integer :: i, ic

    cluster_sizes   = 0

   do i=1, nads
        ic = cluster_label(ads_list(i,1),ads_list(i,2))
        cluster_sizes(ic) = cluster_sizes(ic) + 1
    end do

end subroutine

real(8) function arrhenius(temperature, prefactor, act_energy)

real(8), intent(in) :: temperature, prefactor, act_energy

    arrhenius = prefactor*exp(-act_energy/temperature)


end function
