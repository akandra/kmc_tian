module utilities

  use constants
  use iso_fortran_env

  implicit none

  public

contains

  !-----------------------------------------------------
  ! Random numbers generators
  !-----------------------------------------------------

  function ran1()  !returns random number between 0 - 1

      real(dp) ran1,x
      call random_number(x) ! built in fortran 90 random number function
      ran1 = x

  end function ran1

  function irand(i)  !returns integer random number between 1 and i

      integer irand
      integer i
      real(dp) x

      call random_number(x) ! built in fortran 90 random number function
      irand = floor( i*x ) + 1

  end function irand

  !-----------------------------------------------------
  ! Subs and function dealing with strings
  !-----------------------------------------------------

  function lower_case(input_string) result(output_string)

    character(*), intent(in)      :: input_string
    character(len(input_string))  :: output_string
    integer :: i

    do i = 1, len(input_string)

      select case(input_string(i:i))

        case("A":"Z")
          output_string(i:i) = achar(iachar(input_string(i:i))+32)

        case default
          output_string(i:i) = input_string(i:i)

      end select

    end do

  end function lower_case

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

  subroutine progress_bar( percent_done, symbol )

    integer, intent(in)                   :: percent_done
    character(len=1), intent(in),optional :: symbol

    integer             :: i
    character(len=1)    :: cr = char(13)
    character(len=1)    :: symb
    character(len=56)   :: hdr = "    0   10   20   30   40   50   60   70   80   90  100 "
    character(len=56)   :: bar = "???%|                                                  |"
    logical             :: first = .true.

    symb = '*'
    if(present(symbol)) symb = symbol

    ! construct the percent_done bar
    write(unit=bar(1:3),fmt="(i3)") percent_done
    do i=1, 50
      if(percent_done >= i*2 ) then
        bar(i+5:i+5)= symb
      else
        bar(i+5:i+5)= ' '
      end if
    enddo

  !    ! output the percent_done bar to the screen
  !    if (first) then
  !      write(output_unit, "(10x, a56)") hdr
  !      first=.false.
  !    end if

    write(output_unit, "(a, 10x, a56)", advance='no') cr, bar
    flush(output_unit)

  end subroutine progress_bar


!------------------------------------------------------------------------------
!                 Miscellaneous
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  integer function get_index (key, keylist)
!------------------------------------------------------------------------------
    character(*), intent(in)  :: key
    character(*), intent(in)  :: keylist(:)
    integer                   :: i, index

    index = 0
    do i = 1,size(keylist)
      if(key == keylist(i)) index=i
    end do
    get_index = index

  end function get_index

!------------------------------------------------------------------------------
  subroutine error_message(file_name, line_number, line, message, warning, stop)
!------------------------------------------------------------------------------
    character(*), intent(in)      :: file_name
    integer,      intent(in)      :: line_number
    character(*), intent(in)      :: line
    character(*), intent(in)      :: message
    logical, intent(in), optional :: warning
    logical, intent(in), optional :: stop

    logical                   :: st, warn
    character(100)            :: dummy
    character(:), allocatable :: msg

    ! process optional arguments
    st = .true.
    if(present(stop)) st=stop

    warn = .false.
    if(present(warning)) warn = warning

    if (warn) then
      msg = '      It is my duty to give you gentle warning concerning'
    else
      msg = '      It is my duty to inform you that you have an error'
    end if

    ! convert line number to string
    write(dummy,*) line_number

    print*
    print '(A)',  '--- Dear Sir/Madam: '
    print '(A)',  msg
    print '(2A)', '         file: ', trim(file_name)
    print '(2A)', '         line: ' // trim(adjustl(dummy)) // ': ' // trim(line)
    print '(2A)', '         reason: ',trim(message)
    print '(/A)', '      As always, I remain your humble servant, kMC Code'

    if(st) stop 999

  end subroutine

end module utilities
