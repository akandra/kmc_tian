module utilities

  implicit none

contains

function ran1()  !returns random number between 0 - 1

    real(8) ran1,x
    call random_number(x) ! built in fortran 90 random number function
    ran1 = x

end function ran1

subroutine lower_case(str)

  character(*), intent(in out) :: str
  integer :: i

  do i = 1, len(str)
      select case(str(i:i))
          case("A":"Z")
              str(i:i) = achar(iachar(str(i:i))+32)
      end select
  end do

end subroutine lower_case

subroutine split_string ( line, words, nw )

  character(*), intent(in)  :: line
  character(*), intent(out) :: words(:)
  integer,      intent(out) :: nw
  character(len(words)) :: buf( size(words) )
  integer :: i, ios

  nw = 0 ; words(:) = ""
  do i = 1, size(words)
      read( line, *, iostat=ios ) buf( 1 : i )
      if ( ios /= 0 ) exit
      nw = i
      words( 1 : nw ) = buf( 1 : nw )
  enddo

end subroutine

end module utilities
