program mmc

! load module with random number and absorption rate functions
use utilities
use open_file

implicit none

real(8), parameter :: eV2K = 11604.52
real(8), parameter :: big_bang = 0.0d0 ! in fortnights
real(8), parameter :: pi = 4.0d0*atan(1.0d0)

character(len=3) :: algorithm   ! MC algorithm to use (kmc or mmc)
integer :: nlat         ! size of 2D lattice (nlat x nlat)
real(8) :: temperature  ! temperature in K
real(8) :: coverage     ! coverage in ML per 2D MC lattice (nlat x nlat)
real(8) :: eps          ! O-O interaction energy in eV
integer :: save_period  ! period for conf. output
integer :: hist_period  ! period for histogram  calculation
integer :: nsteps   ! number of Metropolis MC steps
integer :: ntrajs	! number of kmc trajectories (ignored for mmc)
real(8) :: t_end     ! kmc simulation time (s)
integer :: n_bins  	! number of time intervals in kmc histogram
character(len=120) :: rate_file ! file name with rate parameters

real(8) :: energy, total_energy, arrhenius

integer :: nseed = 8
integer, parameter :: seed(8) = (/1,6,3,5,7,3,3,7/)
integer :: icount1, icount2

integer :: nads, nnn, nn_counter, nlat_old, nads_old
integer :: i, j, k, m, n, ihop, istep, i_old, j_old, i_new, j_new, itraj

real(8), dimension(2) :: b1, b2
integer, dimension(:,:), allocatable   :: occupations, ads_list, nn_list
integer, dimension(:),   allocatable   :: nn_opps, temp1D

integer, dimension(:),   allocatable   :: cluster_sizes, hist
integer, dimension(:,:), allocatable   :: cluster_label

real(8) :: time, delta_t, time_new, step_bin
real(8) :: r_hop, rate_acc, u, total_rate
integer :: i_nn, i_ads, ibin_new, ibin, k_change, kmc_nsteps
real(8), dimension(:,:),   allocatable :: rates

real(8) :: energy_old, delta_E, bexp

integer :: ios, pos1, largest_label, hist_counter
character(len=120) buffer, label, fname, cfg_fname

character(len=120) :: rate_key
real(8) :: rate_par1, rate_par2, rtemp

character(len=10) cfg_fmt
character(len=13) ads_fmt

! Read in simulation parameters

select case (iargc())

case(1)
    call getarg(1,fname)
case(2)
    call getarg(1,fname)
    call getarg(2,cfg_fname)
case default
    stop "Wrong number of arguments"
end select

call open_for_read(5, trim(fname)//'.inp' )
ios = 0
do while (ios == 0)

        read(5, '(A)', iostat=ios) buffer
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
            case('temperature')
                read(buffer,*,iostat=ios) temperature
            case('coverage')
                read(buffer,*,iostat=ios) coverage
            case('eps')
                read(buffer,*,iostat=ios) eps
                eps = eps * eV2K
            case('save_period')
                read(buffer,*,iostat=ios) save_period
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
                print *, 'Skipping invalid label at line', label
            end select
        end if

end do ! ios

close(5)

call open_for_read(6,rate_file)

read(6,*) rate_key, rate_par1, rate_par2
rate_par2 = rate_par2 * eV2K

close(6)

! Number of adsorbate particles
nads = nint(coverage*nlat*nlat)
! Number of neighbors for the hexagonal structure
nnn = 6

! configuration input/output format
write(cfg_fmt,'(i6)') nlat
cfg_fmt = '('//trim(adjustl(cfg_fmt))//'i8)'
write(ads_fmt,'(i6)') nads
ads_fmt = '('//trim(adjustl(ads_fmt))//'i8)'

call open_for_write(6,trim(fname)//'.confs')
call open_for_write(7,trim(fname)//'.en')
if (hist_period > 0) call open_for_write(8,trim(fname)//'.hist')

! allocate memory for arrays
allocate(occupations(nlat,nlat), ads_list(nads,2))
allocate(nn_list(nnn,2), nn_opps(nnn), temp1D(nlat*nlat))
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

! Hexagonal transformation matrix

b1 = (/         1.0d0,       0.0d0/)
b2 = (/ cos(pi/3.0d0), sin(pi/3.0d0) /)


do m=1,nnn
    nn_opps(m) = modulo(m-1+nnn/2,nnn)+1
end do

select case (iargc())

case(1)

    temp1D = 0
    do i=1,nads
        temp1D(i) = i
    end do
    occupations = reshape(temp1D,(/nlat,nlat/))

case(2)

    call open_for_read(5,trim(cfg_fname))
    read(5,*) nlat_old, nads_old
    if (nlat_old /= nlat .OR. nads_old /= nads ) &
        stop 'Error: inconsistent input and configuration files'
    read(5,cfg_fmt) (occupations(i,:), i=1,nlat)
    close(5)

case default
    stop "Wrong number of arguments"

end select

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

write(6,*) nlat, nads
write(6,cfg_fmt) (occupations(i,:), i=1,nlat)

write(7,*) total_energy(nlat, nads, nnn, occupations, ads_list, nn_list, eps)

select case (algorithm)

    case ('mmc')

        !call system_clock(icount1)
        do istep=2, nsteps

            do i=1, nads

                energy_old = energy(i, nlat, nads, nnn, occupations, ads_list, nn_list, eps)

                ihop = floor(nnn*ran1()) + 1

                i_new = modulo(ads_list(i,1) + nn_list(ihop,1)-1,nlat) + 1
                j_new = modulo(ads_list(i,2) + nn_list(ihop,2)-1,nlat) + 1

                if (occupations(i_new, j_new) == 0) then

                    i_old = ads_list(i,1)
                    j_old = ads_list(i,2)

                    occupations(i_new,j_new) = i
                    occupations(i_old,j_old) = 0
                    ads_list(i,:) = (/i_new,j_new/)

                    delta_E = energy(i, nlat, nads, nnn, occupations,&
                                      ads_list, nn_list, eps) - energy_old

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
                write(6,cfg_fmt) (occupations(i,:), i=1,nlat)
                write(7,*) total_energy(nlat, nads, nnn, occupations, ads_list, nn_list, eps)
                if (hist_period > 0) then
                    write(8,*) hist_counter
                    write(8,ads_fmt) hist
                end if
           end if

        enddo

    case ('kmc')

        allocate(rates(nads,nnn))

        ! time binning for distributions
        if (n_bins > 0) step_bin = t_end/n_bins

        ! Provisional definition of rates
        r_hop = arrhenius(temperature, rate_par1, rate_par2)/nnn
        print*,"Hopping rate is ", r_hop
        print*

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
                call open_for_write(10,trim(fname)//trim(buffer)//'.csz')
                write(10,*) t_end, n_bins
            end if


            ! Construct rate array
            do i=1,nads
            do m=1,nnn

                i_new = modulo(ads_list(i,1) + nn_list(m,1)-1,nlat) + 1
                j_new = modulo(ads_list(i,2) + nn_list(m,2)-1,nlat) + 1

                if (occupations(i_new, j_new) == 0) then
                    rates(i,m) = r_hop
                else
                    rates(i,m) = 0.0d0
                end if

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

                        write(10,*) ibin
                        write(10,ads_fmt) hist

                    end if

                    ibin = ibin_new

                end if

                time = time_new ! time shift


                ! Update rate constants
                ! Particle i_ads is going to hop in direction i_nn

                ! scan over old neighbors
                do m=1,nnn
                    ! (i,j) is a position of neighbor m
                    i = modulo(ads_list(i_ads,1) + nn_list(m,1)-1,nlat) + 1
                    j = modulo(ads_list(i_ads,2) + nn_list(m,2)-1,nlat) + 1
                    ! update rate for neighbor m in direction nn_opps(m)
                    if (occupations(i,j) > 0) rates(occupations(i,j),nn_opps(m)) = r_hop

                end do

                ! Update adsorbate position after hop
                i_new = modulo(ads_list(i_ads,1) + nn_list(i_nn,1)-1,nlat) + 1
                j_new = modulo(ads_list(i_ads,2) + nn_list(i_nn,2)-1,nlat) + 1
                ! Update occupations and ads. list
                occupations(ads_list(i_ads,1),ads_list(i_ads,2)) = 0
                occupations(i_new,j_new) = i_ads
                ads_list(i_ads,:) = (/i_new,j_new/)

                ! scan over new neighbors
                do m=1,nnn
                    ! (i,j) is a position of neighbor m
                    i = modulo(i_new + nn_list(m,1)-1,nlat) + 1
                    j = modulo(j_new + nn_list(m,2)-1,nlat) + 1
                    ! update rate for neighbor m in direction nn_opps(m)
                    ! as well as rates for particle i_ads which committed a hop
                    if (occupations(i,j) > 0) then
                        rates(occupations(i,j),nn_opps(m)) = 0.0d0
                        rates(i_ads,m) = 0.0d0
                    else
                        rates(i_ads,m) = r_hop
                    end if

                end do

            end do ! over time
!-------------------------------------------------------------

            if (n_bins > 0) close(10)

            print*,'Number of kmc-steps is ', kmc_nsteps
            print*

        end do ! over trajectories

        deallocate(rates)

    case default
        stop 'Error: mc algorithm not defined'
end select


close(6)
close(7)
if (hist_period > 0) close(8)

call system_clock(icount2)
print*
print*, 'Number of clock ticks = ', icount2-icount1

call open_for_write(6,trim(fname)//'.out')

write(6,*) nlat, nads
write(6,cfg_fmt) (occupations(i,:), i=1,nlat)

close(6)

deallocate(cluster_label, cluster_sizes, hist)
deallocate(ads_list,nn_list,nn_opps,temp1D,occupations)


end program

real(8) function energy(inx, nlat, nads, nnn, occupations, ads_list, nn_list, eps)

integer, intent(in) :: inx, nlat, nads, nnn
integer, dimension(nlat,nlat), intent(in) :: occupations
integer, dimension(nads,2), intent(in) :: ads_list
integer, dimension(nnn,2), intent(in) :: nn_list
real(8), intent(in) :: eps

integer :: i, ic, jc, counter

    counter = 0
    do i=1, nnn

        ic = modulo(ads_list(inx,1)+nn_list(i,1)-1,nlat) + 1
        jc = modulo(ads_list(inx,2)+nn_list(i,2)-1,nlat) + 1
        if (occupations(ic,jc) > 0) counter = counter + 1

    end do

    energy = eps*counter

end function

real(8) function total_energy(nlat, nads, nnn, occupations, ads_list, nn_list, eps)

integer, intent(in) :: nlat, nads, nnn
integer, dimension(nlat,nlat), intent(in) :: occupations
integer, dimension(nads,2), intent(in) :: ads_list
integer, dimension(nnn,2), intent(in) :: nn_list
real(8), intent(in) :: eps

integer :: i, ic, jc, counter

    counter = 0
    do inx=1,nads
    do i=1, nnn/2 ! count interaction with neighbors once

        ic = modulo(ads_list(inx,1)+nn_list(i,1)-1,nlat) + 1
        jc = modulo(ads_list(inx,2)+nn_list(i,2)-1,nlat) + 1
        if (occupations(ic,jc) > 0) counter = counter + 1

    end do
    end do

    total_energy = counter*eps

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
