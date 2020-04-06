program mmc

! load module with random number and absorption rate functions
use utilities
use open_file

implicit none

real(8), parameter :: eV2K = 11604.52

integer :: nlat         ! size of 2D lattice (nlat x nlat)
real(8) :: temperature  ! temperature in K
real(8) :: coverage     ! coverage in ML
real(8) :: eps          ! O-O interaction energy in eV
integer :: nsteps       ! number of Metropolis MC steps
integer :: nsave        ! period for for conf. output

real(8) :: energy, total_energy

integer :: nseed = 8
integer, parameter :: seed(8) = (/1,6,3,5,7,3,3,7/)
integer :: icount1, icount2

integer :: nads, nnn, nn_counter, nlat_old, nads_old
integer :: i, j, k, ihop, istep, i_old, j_old, i_new, j_new
integer, dimension(:,:), allocatable   :: occupations, occupations_new
integer, dimension(:),   allocatable   :: temp1D
integer, dimension(:,:), allocatable   :: ads_list,ads_list_new, nn_list

real(8) :: energy_old, delta_E, bexp

integer :: ios, pos1
character(len=120) buffer, label, fname, cfg_fname
character(len=10) cfg_fmt


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
            case('nlat')
                read(buffer,*,iostat=ios) nlat
            case('temperature')
                read(buffer,*,iostat=ios) temperature
            case('coverage')
                read(buffer,*,iostat=ios) coverage
            case('eps')
                read(buffer,*,iostat=ios) eps
                eps = eps * eV2K
            case('nsteps')
                read(buffer,*,iostat=ios) nsteps
            case('nsave')
                read(buffer,*,iostat=ios) nsave
            case default
                print *, 'Skipping invalid label at line', label
            end select
        end if

end do ! ios

close(5)

! configuration input/output format
write(cfg_fmt,'(i6)') nlat
cfg_fmt = '('//trim(adjustl(cfg_fmt))//'i3)'

call open_for_write(6,trim(fname)//'.confs')
call open_for_write(7,trim(fname)//'.en')

nads = nint(coverage*nlat*nlat)

! Number of neighbors for the hexagonal structure
nnn = 6

! allocate memory for arrays
allocate(occupations(nlat,nlat), temp1D(nlat*nlat), ads_list(nads,2))
allocate(occupations_new(nlat,nlat), ads_list_new(nads,2))
allocate(nn_list(nnn,2))

! NN list for the hexagonal structure
!  11*   12*   13*   14
!
!     21*   22*   23*   24
!
!        31*   32*   33    34
!
!           41*   42*   43    44

nn_list(1,:) = (/ 0, 1/)
nn_list(2,:) = (/ 1, 0/)
nn_list(3,:) = (/ 1,-1/)
nn_list(4,:) = (/ 0,-1/)
nn_list(5,:) = (/-1, 0/)
nn_list(6,:) = (/-1, 1/)

select case (iargc())

case(1)

    temp1D = 0
    temp1D(1:nads) = 1
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

k = 0
do j=1,nlat
do i=1,nlat
    if (occupations(i,j) == 1) then
        k = k + 1
        ads_list(k,1) = i
        ads_list(k,2) = j
    end if
end do
end do

write(6,*) nlat, nads
write(6,cfg_fmt) (occupations(i,:), i=1,nlat)

write(7,*) total_energy(nlat, nads, nnn, occupations, ads_list, nn_list, eps)


call system_clock(icount1)
do istep=2, nsteps

    do i=1, nads

        energy_old = energy(i, nlat, nads, nnn, occupations, ads_list, nn_list, eps)

        ihop = floor(nnn*ran1()) + 1

        i_new = modulo(ads_list(i,1) + nn_list(ihop,1)-1,nlat) + 1
        j_new = modulo(ads_list(i,2) + nn_list(ihop,2)-1,nlat) + 1

        if (occupations(i_new, j_new) == 0) then

            i_old = ads_list(i,1)
            j_old = ads_list(i,2)

            occupations(i_old,j_old) = 0
            occupations(i_new,j_new) = 1
            ads_list(i,:) = (/i_new,j_new/)

            delta_E = energy(i, nlat, nads, nnn, occupations,&
                              ads_list, nn_list, eps) - energy_old

            if (exp(- delta_E/temperature) < ran1()) then

                occupations(i_old,j_old) = 1
                occupations(i_new,j_new) = 0
                ads_list(i,:) = (/i_old,j_old/)

            end if

        end if

    enddo

    if (mod(istep, nsave) == 0) then
        write(6,cfg_fmt) (occupations(i,:), i=1,nlat)
        write(7,*) total_energy(nlat, nads, nnn, occupations, ads_list, nn_list, eps)
    end if

    if (mod(istep, 1000) == 0) then
        print*, istep
!        do i=1,nlat
!        write(*,*) occupations(i,:)
!        enddo
    end if

enddo

close(6)
close(7)

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

deallocate(occupations_new, ads_list_new)
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
        counter = counter + occupations(ic,jc)

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
        counter = counter + occupations(ic,jc)

    end do
    end do

    total_energy = counter*eps

end function
