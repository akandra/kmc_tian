program mmc

! load module with random number and absorption rate functions
use utilities
use open_file

implicit none

real(8), parameter :: eV2K = 11604.52

integer :: nlat         ! size of 2D lattice (nlat x nlat)
real(8) :: temperature  ! temperature in K
real(8) :: coverage     ! coverage in ML

integer :: nseed = 8
integer, parameter :: seed(8) = (/1,6,3,5,7,3,3,7/)

integer :: nads
integer :: i, j, k
integer, dimension(:,:), allocatable   :: occupations
integer, dimension(:),   allocatable   :: temp1D

integer :: ios, pos1
character(len=120) buffer, label, inputfname


! Read in simulation parameters

if (iargc() == 0) stop "I need an input file with parameters"
call getarg(1,inputfname)

call open_for_read(5, inputfname)
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
            case default
                print *, 'Skipping invalid label at line', label
            end select
        end if

end do ! ios

close(5)

! allocate memory for arrays
allocate(occupations(nlat,nlat), temp1D(nlat*nlat))

nads = nint(coverage*nlat*nlat)

temp1D = 0
temp1D(1:nads) = 1
occupations = reshape(temp1D,(/nlat,nlat/))

do i=1, nlat
    print*, occupations(:,i)
enddo

stop

!k=0
!loop1:do i=1, nlat
!do j=1, nlat
!
!    occupations(i,j) = 1
!    k = k + 1
!    if (k == nads) exit loop1
!
!end do
!end do loop1
! rate constants

!kdiff = 5.0d9 ! in s-1
!Ediff = 0.43*eV2K

deallocate(temp1D,occupations)

end program


