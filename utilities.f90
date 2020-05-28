module utilities

  use constants

  implicit none

contains

!-----------------------------------------------------
! Random numbers generators
!-----------------------------------------------------

function ran1()  !returns random number between 0 - 1

    real(8) ran1,x
    call random_number(x) ! built in fortran 90 random number function
    ran1 = x

end function ran1

!-----------------------------------------------------
! Subs and function dealing with strings
!-----------------------------------------------------

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

subroutine split_string ( line, words, nw, comment_character )

  character(*), intent(in)  :: line
  character(*), intent(out) :: words(:)
  integer,      intent(out) :: nw
  character, optional       :: comment_character

  character(len(words)) :: buf( size(words) )
  character :: cc
  integer :: i, ios

  if (present(comment_character)) then
    cc = comment_character
  else
    cc = '!'
  end if

  nw = 0 ; words(:) = ""
  do i = 1, size(words)
      read( line, *, iostat=ios ) buf( 1 : i )
      if ( ios /= 0 ) exit
      if ( buf(i)(1:1) == cc ) exit
      nw = i
      words( 1 : nw ) = buf( 1 : nw )
  enddo

end subroutine

logical function read_num(string, x)
  implicit none

  character(len=*), intent(in) :: string
  real(dp), intent(out) :: x
  integer :: error_code

  read(string,*,iostat=error_code) x
  read_num = error_code == 0

end function read_num

!-----------------------------------------------------
! Miscellaneous
!-----------------------------------------------------

integer function get_index (key, keylist)
  character(*), intent(in)  :: key
  character(*), intent(in)  :: keylist(:)
  integer                   :: i, index

  index = 0
  do i = 1,size(keylist)
    if(key == keylist(i)) index=i
  end do
  get_index = index

end function get_index

subroutine error(file_name, line_number, line, message)

  character(*), intent(in) :: file_name, line, message
  character(100) :: dummy
  integer, intent(in) :: line_number

  print*
  print*, '---Dear Sir/Madam, '
  print*, 'It is my duty to inform you that you have an error'
  print*, '  in file ', trim(file_name)
  write(dummy,*) line_number
  print*, '  line ' // trim(adjustl(dummy)) // ': ' // trim(line)
  print*, '  reason: ',trim(message)
  print*, 'Remaining your humble servant, K.M.C. Code'
  stop

end subroutine

end module utilities
