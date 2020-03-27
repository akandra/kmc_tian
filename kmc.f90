program kmc

! load module with random number and absorption rate functions
use utilities
use open_file

implicit none

integer :: nlat ! size of 2D lattice (nlat x nlat)
integer :: vmax        ! highest vibrational level considered
real(8) :: Temp        ! temperature in K
real(8) :: t_end       ! length of trajectory in seconds
integer :: start_traj, ntrajs     ! 1st traj and number of trajectories
integer :: n_bins      ! number of time bins the trajectory is devided in

integer :: nseed = 8
integer, parameter :: seed(8) = (/1,6,3,5,7,3,3,7/)

integer :: i,j,k,iproc,i_hop,j_hop, i_hop2, j_hop2, itraj,ibin,ibin_new

integer, dimension(:,:), allocatable :: vib_state     ! for each lattice point a vibrational state is defined
integer, dimension(2) :: inlist, jnlist

real(8), dimension(:),   allocatable :: wn, wn1, kn, freqs, energy
real(8), dimension(:,:), allocatable :: wnm, wnmr
real(8) :: rate, time, time_new, u, total_rate, delta_t
real(8), dimension(:,:,:), allocatable   :: Rint

real(8), dimension(:,:), allocatable :: vib_d
real(8) :: step_bin,fract

integer :: ios, pos1
character(len=120) buffer, label, distfname, inputfname


! Read in simulation parameters

if (iargc() == 0) stop "I need an input file with parameters"
call getarg(1,inputfname)

open(38, file=inputfname, status='old', action='read', iostat = ios)

if (ios == 0) then

    do while (ios == 0)
        read(38, '(A)', iostat=ios) buffer
        if (ios == 0) then

        ! Find the first instance of whitespace.  Split label and data.
            pos1 = scan(buffer, ' ')
            label = buffer(1:pos1)
            buffer = buffer(pos1+1:)

            select case (label)
            case('start_traj')
                read(buffer,*,iostat=ios) start_traj
            case('ntrajs')
                read(buffer,*,iostat=ios) ntrajs
            case('nlat')
                read(buffer,*,iostat=ios) nlat
            case('vmax')
                read(buffer,*,iostat=ios) vmax
            case('temperature')
                read(buffer,*,iostat=ios) Temp
            case('t_end')
                read(buffer,*,iostat=ios) t_end
            case('n_bins')
                read(buffer,*,iostat=ios) n_bins
            case default
                print *, 'Skipping invalid label at line', label
            end select
        end if
    end do ! ios

    close(38)

else
    print *, 'Input file does not exist'
    stop
end if

! allocate memory for arrays
allocate(vib_state(nlat,nlat),Rint(nlat,nlat,7), vib_d(n_bins,0:vmax))
allocate(wn(0:vmax), wn1(0:vmax), kn(0:vmax), freqs(0:vmax), energy(0:vmax))
allocate(wnm(0:vmax,0:vmax), wnmr(0:vmax,0:vmax))

! read in rate constants
open(11,file='Ratematrix.dat', status='old', action='read', iostat = ios)
if (ios==0) then
    read(11,*) (wn(i),wn1(i),kn(i),freqs(i),energy(i), i=0,vmax)
    close(11)
else
    print *, 'I need an input file named Ratematrix.dat'
    stop
end if

! test value for resonant V-V rates to speed up calculations
! should be commented out after test phase
!wn1(0)=6.0d7

! define V-V rate matrix

freqs = exp(-1.43878631*freqs/Temp) ! detailed balance Boltzmann factor

wnm = 0.d0
wnm(0:vmax,1) = wn1                 ! V-V energy transfer rates
do i=1,vmax-1
    wnmr(0,i+1) = wnm(i,1)*freqs(i) ! backward rates
end do
wnm = wnm + wnmr                    ! both are stored at different positions
                                    ! in the V-V rate matrix
! time binning for distributions
step_bin = t_end/n_bins

! loop over trajectories
do itraj=start_traj,start_traj+ntrajs-1

    ! initialize random number generator
    call random_seed(size=nseed)
    call random_seed(put=itraj*seed)

!    write (*,*) itraj

    ! set initial values
    time = 0.d0
    vib_state = 0       ! start with no excitation
    ibin = 0            ! binning interval counter
    fract = 0.d0        ! fraction of time in time bin
    vib_d = 0.d0        ! vibrational population distribution

! debug output
!    write(*,*) itraj
!    write(*,'(f18.15)')time
!    write(*,'(6i3)') iproc,i_hop,j_hop,i_hop2,j_hop2
!    write(*,*)
!    write(*,'(4i3)')vib_state
!    pause

    ! initialize all possible rate processes for all molecules of the lattice
    !   1st order rate processes
    do i=1,nlat
    do j=1,nlat
        Rint(i,j,1) = kabs(vib_state(i,j),time) ! excitation by light absorption
        Rint(i,j,2) = kn(vib_state(i,j))        ! fluorescence
        Rint(i,j,3) = wn(vib_state(i,j))        ! VER
    end do
    end do

    !   2nd order rate processes
    do j=1,nlat
        do i=1,nlat-1
            Rint(i,j,4) = wnm(vib_state(i,j),vib_state(i+1,j))  ! V-V
            Rint(i,j,5) = wnm(vib_state(i+1,j),vib_state(i,j))
            Rint(j,i,6) = wnm(vib_state(j,i),vib_state(j,i+1))
            Rint(j,i,7) = wnm(vib_state(j,i+1),vib_state(j,i))
        end do
        Rint(nlat,j,4) = wnm(vib_state(nlat,j),vib_state(1,j))  ! V-V at the boundary
        Rint(nlat,j,5) = wnm(vib_state(1,j),vib_state(nlat,j))
        Rint(j,nlat,6) = wnm(vib_state(j,nlat),vib_state(j,1))
        Rint(j,nlat,7) = wnm(vib_state(j,1),vib_state(j,nlat))
    end do

    ! start time propagation
    do while (time<t_end)

        total_rate = sum(Rint)  ! calculate total rate
        u = ran1()*total_rate   ! random number to select a process

        rate = 0.d0
        extloop: do i=1,nlat    ! determine reaction channel
            do j=1,nlat
                do iproc=1,7
                    rate = rate + Rint(i,j,iproc)
                    if (u < rate) then  ! site (i_hop,j_hop) where
                        i_hop=i         ! the process (iproc) occurs
                        j_hop=j
                        exit extloop
                    end if
                end do
            end do
        end do extloop

        delta_t = -log(ran1())/total_rate   ! when does a hop occur?
        time_new = time + delta_t
        if (time_new > t_end) time_new = t_end

        ! histogram a vibrational population distribution
        ibin_new = int(time_new/step_bin)   ! do we pass a time bin boundary?
        if (ibin_new == ibin) then          ! if not

            fract = delta_t/step_bin
            do i=1,nlat
            do j=1,nlat
                vib_d(ibin,vib_state(i,j)) = vib_d(ibin,vib_state(i,j)) + fract
            end do
            end do

       else                                 ! if so

            fract = ibin - time/step_bin
            do i=1,nlat
            do j=1,nlat
                vib_d(ibin,vib_state(i,j)) = vib_d(ibin,vib_state(i,j)) + fract
            end do
            end do

            fract = time_new/step_bin - ibin_new
            do i=1,nlat
            do j=1,nlat
                do k=ibin,ibin_new-1
                    vib_d(k,vib_state(i,j)) = vib_d(k,vib_state(i,j)) + 1
                end do
                vib_d(ibin_new,vib_state(i,j)) = &
                                vib_d(ibin_new,vib_state(i,j)) + fract
            end do
            end do
            ibin = ibin_new

        endif

    time = time_new ! time shift

! debug output
!        write(*,*) time, ibin
!        write(*,'(11f10.5)') ,vib_d(1,0:10)
!        pause

        ! update vibrational state population after hop
        select case (iproc)
            case (1) ! absorption
                vib_state(i_hop,j_hop) = vib_state(i_hop,j_hop) + 1
            case (2) ! fluoresence
                vib_state(i_hop,j_hop) = vib_state(i_hop,j_hop) - 1
            case (3) ! VER
                vib_state(i_hop,j_hop) = vib_state(i_hop,j_hop) - 1
            case (4) ! horizontal pooling + resonant V-V
                if (i_hop == nlat) then
                    i_hop2 = 1
                else
                    i_hop2 = i_hop + 1
                end if
                j_hop2 = j_hop
                vib_state(i_hop,j_hop)   =   vib_state(i_hop,j_hop) + 1
                vib_state(i_hop2,j_hop2) = vib_state(i_hop2,j_hop2) - 1

            case (5) ! neighbour hor/ntrajsizontal pooling + resonant V-V
                if (i_hop == nlat) then
                    i_hop2 = 1
                else
                    i_hop2 = i_hop + 1
                end if
                j_hop2 = j_hop
                vib_state(i_hop,j_hop)   =   vib_state(i_hop,j_hop) - 1
                vib_state(i_hop2,j_hop2) = vib_state(i_hop2,j_hop2) + 1

            case (6) ! vertical pooling + resonant V-V
                if (j_hop == nlat) then
                    j_hop2 = 1
                else
                    j_hop2 = j_hop + 1
                end if
                i_hop2 = i_hop
                vib_state(i_hop,j_hop)   =   vib_state(i_hop,j_hop) + 1
                vib_state(i_hop2,j_hop2) = vib_state(i_hop2,j_hop2) - 1

            case (7) ! neighbour vertical pooling + resonant V-V
                if (j_hop == nlat) then
                    j_hop2 = 1
                else
                    j_hop2 = j_hop + 1
                end if
                i_hop2 = i_hop
                vib_state(i_hop,j_hop) =  vib_state(i_hop,j_hop) - 1
                vib_state(i_hop2,j_hop2) = vib_state(i_hop2,j_hop2) + 1

        end select

    ! debug output
!    write(*,'(f18.15)')time
!    write(*,'(6i3)') iproc,i_hop,j_hop,i_hop2,j_hop2
!    write(*,*)
!    write(*,'(4i3)')vib_state
!    pause

        ! Update rate constants for the new state of the system
        ! site (i_hop, j_hop)

        !   1st order rate processes
        Rint(i_hop,j_hop,1) = kabs(vib_state(i_hop,j_hop),time)
        Rint(i_hop,j_hop,2) = kn(vib_state(i_hop,j_hop))
        Rint(i_hop,j_hop,3) = wn(vib_state(i_hop,j_hop))

        ! Periodic boundary conditions
        if (i_hop==1) then
            inlist=(/nlat,i_hop+1/)
        elseif (i_hop==nlat) then
            inlist=(/i_hop-1,1/)
        else
            inlist=(/i_hop-1,i_hop+1/)
        end if
        if (j_hop==1) then
            jnlist=(/nlat,j_hop+1/)
        elseif (j_hop==nlat) then
            jnlist=(/j_hop-1,1/)
        else
            jnlist=(/j_hop-1,j_hop+1/)
        end if

        ! 2nd order rate processes

        Rint(i_hop,j_hop,4) = &
            wnm(vib_state(i_hop,j_hop),vib_state(inlist(2),j_hop))
        Rint(i_hop,j_hop,5) = &
            wnm(vib_state(inlist(2),j_hop),vib_state(i_hop,j_hop))

        Rint(inlist(1),j_hop,4) = &
            wnm(vib_state(inlist(1),j_hop),vib_state(i_hop,j_hop))
        Rint(inlist(1),j_hop,5) = &
            wnm(vib_state(i_hop,j_hop),vib_state(inlist(1),j_hop))


        Rint(i_hop,j_hop,6) = &
            wnm(vib_state(i_hop,j_hop),vib_state(i_hop,jnlist(2)))
        Rint(i_hop,j_hop,7) = &
            wnm(vib_state(i_hop,jnlist(2)),vib_state(i_hop,j_hop))

        Rint(i_hop,jnlist(1),6) = &
            wnm(vib_state(i_hop,jnlist(1)),vib_state(i_hop,j_hop))
        Rint(i_hop,jnlist(1),7) = &
            wnm(vib_state(i_hop,j_hop),vib_state(i_hop,jnlist(1)))

        ! Update rate constants for the new state of the system
        ! site (i_hop2, j_hop2)
        if (iproc > 3) then

            !   1st order rate processes
            Rint(i_hop2,j_hop2,1) = kabs(vib_state(i_hop2,j_hop2),time)
            Rint(i_hop2,j_hop2,2) = kn(vib_state(i_hop2,j_hop2))
            Rint(i_hop2,j_hop2,3) = wn(vib_state(i_hop2,j_hop2))

            !   Periodic boundary conditions
            if (i_hop2==1) then
                inlist=(/nlat,i_hop2+1/)
            elseif (i_hop2==nlat) then
                inlist=(/i_hop2-1,1/)
            else
                inlist=(/i_hop2-1,i_hop2+1/)
            end if
            if (j_hop2==1) then
                jnlist=(/nlat,j_hop2+1/)
            elseif (j_hop2==nlat) then
                jnlist=(/j_hop2-1,1/)
            else
                jnlist=(/j_hop2-1,j_hop2+1/)
            end if

            !   2nd order rate processes

            Rint(i_hop2,j_hop2,4) = &
                wnm(vib_state(i_hop2,j_hop2),vib_state(inlist(2),j_hop2))
            Rint(i_hop2,j_hop2,5) = &
                wnm(vib_state(inlist(2),j_hop2),vib_state(i_hop2,j_hop2))

            Rint(inlist(1),j_hop2,4) = &
                wnm(vib_state(inlist(1),j_hop2),vib_state(i_hop2,j_hop2))
            Rint(inlist(1),j_hop2,5) = &
                wnm(vib_state(i_hop2,j_hop2),vib_state(inlist(1),j_hop2))


            Rint(i_hop2,j_hop2,6) = &
                wnm(vib_state(i_hop2,j_hop2),vib_state(i_hop2,jnlist(2)))
            Rint(i_hop2,j_hop2,7) = &
                wnm(vib_state(i_hop2,jnlist(2)),vib_state(i_hop2,j_hop2))

            Rint(i_hop2,jnlist(1),6) = &
                wnm(vib_state(i_hop2,jnlist(1)),vib_state(i_hop2,j_hop2))
            Rint(i_hop2,jnlist(1),7) = &
                wnm(vib_state(i_hop2,j_hop2),vib_state(i_hop2,jnlist(1)))

        end if

    end do ! over time

    ! save data to file
    write(distfname,'(A,I6.6,A)') "vdist",itraj,".dat"
    open(11,file=distfname)
    write(11,'(a18,10000i18)') "time_s", (i, i=0,vmax)
    do i=1, n_bins
        write(11,'(10000e18.5)') (i-0.5)*step_bin,vib_d(i,:)
    end do
    close(11)

end do ! over trajs

deallocate(wnm,wnmr)
deallocate(wn,wn1,kn,freqs,energy)
deallocate(vib_state)
deallocate(Rint)
deallocate(vib_d)

end program


