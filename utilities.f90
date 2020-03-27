module utilities

contains

function ran1()  !returns random number between 0 - 1
 implicit none
        real(8) ran1,x
        call random_number(x) ! built in fortran 90 random number function
        ran1 = x
end function ran1

function kabs(v,t)
implicit none
real(8) :: kabs
integer :: v
real(8) :: t

    if ((t<=5.d-6).and.(v==0)) then
        kabs = 9.d4
    else
        kabs = 0.d0
    end if

end function kabs

end module utilities
