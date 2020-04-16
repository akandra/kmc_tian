program mmc

! load module with random number and absorption rate functions
use utilities
use open_file

implicit none

real(8), parameter :: eV2K = 11604.52
real(8), parameter :: big_bang = 0.0d0 ! in fortnights

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

real(8) :: energy, total_energy

integer :: nseed = 8
integer, parameter :: seed(8) = (/1,6,3,5,7,3,3,7/)
integer :: icount1, icount2

integer :: nads, nnn, nn_counter, nlat_old, nads_old
integer :: i, j, k, m, n, ihop, istep, i_old, j_old, i_new, j_new, itraj
integer, dimension(:,:), allocatable   :: occupations, ads_list, nn_list
integer, dimension(:),   allocatable   :: temp1D

integer, dimension(:),   allocatable   :: cluster_sizes, hist
integer, dimension(:,:), allocatable   :: cluster_label

real(8) :: time, delta_t, time_new, step_bin
real(8) :: r_hop, rate_acc, u, total_rate
integer :: i_nn, i_ads, ibin_new, ibin, k_change
real(8), dimension(:,:),   allocatable :: rates
integer, dimension(:), allocatable :: change_list

real(8) :: energy_old, delta_E, bexp

integer :: ios, pos1, largest_label, hist_counter
character(len=120) buffer, label, fname, cfg_fname
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
            case('hist_period')
                read(buffer,*,iostat=ios) hist_period
            case('mmc_nsteps')
                read(buffer,*,iostat=ios) nsteps
            case('kmc_ntrajs')
                read(buffer,*,iostat=ios) ntrajs
            case('kmc_time')
                read(buffer,*,iostat=ios) t_end
            case('kmc_nbins')
                read(buffer,*,iostat=ios) n_bins
            case default
                print *, 'Skipping invalid label at line', label
            end select
        end if

end do ! ios

close(5)

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
allocate(nn_list(nnn,2), temp1D(nlat*nlat))
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

! Provisional declaration block
! Put it later to the beginning

allocate(rates(nads,nnn), change_list(2*nnn))

        ! time binning for distributions
        step_bin = t_end/n_bins

        ! Provisional definition of rates
        r_hop = 1

        do itraj=1, ntrajs

            ! initialize random number generator
            call random_seed(size=nseed)
            call random_seed(put=itraj*seed)

            ibin = 0

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

                ! histogram
                if (time_new > t_end) time_new = t_end

!                ibin_new = int(time_new/step_bin)   ! do we pass a time bin boundary?
!                if (ibin_new == ibin) then          ! if not
!
!                    fract = delta_t/step_bin
!                    do i=1,nlat
!                    do j=1,nlat
!                        vib_d(ibin,vib_state(i,j)) = vib_d(ibin,vib_state(i,j)) + fract
!                    end do
!                    end do
!
!                else                                 ! if so
!
!                    fract = ibin - time/step_bin
!                    do i=1,nlat
!                    do j=1,nlat
!                        vib_d(ibin,vib_state(i,j)) = vib_d(ibin,vib_state(i,j)) + fract
!                    end do
!                    end do
!
!                    fract = time_new/step_bin - ibin_new
!                    do i=1,nlat
!                    do j=1,nlat
!                        do k=ibin,ibin_new-1
!                            vib_d(k,vib_state(i,j)) = vib_d(k,vib_state(i,j)) + 1
!                        end do
!                        vib_d(ibin_new,vib_state(i,j)) = &
!                                        vib_d(ibin_new,vib_state(i,j)) + fract
!                    end do
!                    end do
!                    ibin = ibin_new
!
!            endif

    time = time_new ! time shift


        ! update the state after hop
        change_list = 0
        ! accounting for old neighbors
        do m=1,nnn
            i = modulo(ads_list(i_ads,1) + nn_list(m,1)-1,nlat) + 1
            j = modulo(ads_list(i_ads,2) + nn_list(m,2)-1,nlat) + 1
            change_list(m) = occupations(i,j)
        end do

        ! Adsorbate position after hop
        i_new = modulo(ads_list(i_ads,1) + nn_list(i_nn,1)-1,nlat) + 1
        j_new = modulo(ads_list(i_ads,2) + nn_list(i_nn,2)-1,nlat) + 1
        ! Update occupations and ads. list
        occupations(ads_list(i_ads,1),ads_list(i_ads,2)) = 0
        occupations(i_new,j_new) = i_ads
        ads_list(i_ads,:) = (/i_new,j_new/)

        ! accounting for new neighbors
        do m=1,nnn
            i = modulo(i_new + nn_list(m,1)-1,nlat) + 1
            j = modulo(j_new + nn_list(m,2)-1,nlat) + 1
            change_list(m+nnn) = occupations(i,j)
        end do


        ! Update rate constants for the new state of the system

        do k=1,2*nnn

            i = change_list(k)

            if (i > 0) then

                do m=1,nnn

                    i_new = modulo(ads_list(i,1) + nn_list(m,1)-1,nlat) + 1
                    j_new = modulo(ads_list(i,2) + nn_list(m,2)-1,nlat) + 1

                    if (occupations(i_new, j_new) == 0) then
                        rates(i,m) = r_hop
                    else
                        rates(i,m) = 0.0d0
                    end if

                end do

            end if

        end do

    end do ! over time
!-------------------------------------------------------------

    enddo ! over trajectories

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

! rate constants

!kdiff = 5.0d9 ! in s-1
!Ediff = 0.43*eV2K

deallocate(rates,change_list)
deallocate(cluster_label, cluster_sizes, hist)
deallocate(ads_list,nn_list,temp1D,occupations)


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

!write(*,cfg_fmt) (occupations(m,:), m=1,nlat)
largest_label = 0
do i=1,nlat
do j=1,nlat

    if (occupations(i,j) > 0) then

        do m=1,nnn2

            i_nn(m) = modulo(i + nn_list(nnn2+m,1)-1,nlat) + 1
            j_nn(m) = modulo(j + nn_list(nnn2+m,2)-1,nlat) + 1

            if (occupations(i_nn(m), j_nn(m)) > 0) scanned_nn_occs(m) = 1

            if (i==1) scanned_nn_occs(2:3) = 0
            if (j==1) scanned_nn_occs(1) = 0
            if (j==nlat) scanned_nn_occs(3) = 0

        end do

        select case (sum(scanned_nn_occs))

        case (0)

            largest_label = largest_label + 1
            cluster_label(i,j) = largest_label

        case (1)

            do m=1,nnn2
                if (scanned_nn_occs(m) == 1) then
                    cluster_label(i,j) = &
                        lfind(cluster_label(i_nn(m),j_nn(m)), labels, nads)
                end if
            end do

        case (2)

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
end do
end do

!print*
!write(*,cfg_fmt) (cluster_label(m,:), m=1,nlat)
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
