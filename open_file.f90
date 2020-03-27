module open_file
    !
    ! Purpose:
    !           This module contains routines to open files smoothly
    !
    ! Date          	Author          	History of Revison
    ! ====          	======          	==================
    ! 18.02.2014    	Svenja M. Janke		Original
    !			        Sascha Kandratsenka
    !			        Dan J. Auerbach

    interface
        subroutine open_for_read(lun, file_name)
            integer, intent(in)     :: lun
            character, intent(in)   :: file_name
        end subroutine

        subroutine open_for_write(lun, file_name)
            integer, intent(in)     :: lun
            character, intent(in)   :: file_name
        end subroutine

        subroutine open_for_append(lun, file_name)
            integer, intent(in)     :: lun
            character, intent(in)   :: file_name
        end subroutine

 !       contains module procedure

    endinterface

contains

end module

subroutine open_for_read(lun,file_name)
    implicit none
    integer, intent(in)           :: lun
    character(len=*), intent(in)  :: file_name

    integer                       :: ios
    character(len=120)            :: error_message

    open(unit=lun, file=file_name, status='old', action='read', iostat=ios, iomsg=error_message)
    if (ios==0) return
    print '(/ "Error on open file ", (a), " for read i/o status=", i4 )', TRIM(file_name), ios
    print '( "error message=", (a) )', error_message

    STOP 101

end subroutine open_for_read

subroutine open_for_write(lun,file_name)
    implicit none
    integer, intent(in)           :: lun
    character(len=*), intent(in)  :: file_name

    integer                       :: ios
    character(120)                :: error_message

    open(unit=lun, file=file_name, status='new', action='write', iostat=ios, iomsg=error_message)
    if (ios==0) return

    open(unit=lun, file=file_name, status='replace', action='write', iostat=ios, iomsg=error_message)

    if (ios==0) return
    print *, 'failed to open file ', file_name, ' for write with status=replace. i/o status =',ios
    print *, 'error message: ', error_message
    STOP 103

end subroutine open_for_write

subroutine open_for_append(lun,file_name)
    implicit none
    integer, intent(in)           :: lun
    character(len=*), intent(in)  :: file_name

    integer                       :: ios
    character(120)                :: error_message

    open(unit=lun, file=file_name, status='new', action='write', iostat=ios, iomsg=error_message)
    if (ios==0) return

    open(unit=lun, file=file_name, status='old', access='append', action='write', &
         iostat=ios, iomsg=error_message)

    if (ios==0) return
    print *, 'failed to open file ', file_name, ' for write with status=old. i/o status =',ios
    print *, 'error message: ', error_message
    STOP 103

end subroutine open_for_append

