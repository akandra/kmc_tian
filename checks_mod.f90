module checks

use constants
use mc_lat_class
use control_parameters_class
use energy_parameters_class

implicit none

private
public:: code_checks

contains

subroutine code_checks(lat, c_pars, e_pars)

  type(mc_lat), intent(inout) :: lat
  type(control_parameters), intent(in) :: c_pars
  type(energy_parameters ), intent(in)    :: e_pars


  call check_nn_new_list(lat)

end subroutine

subroutine check_nn_new_list(lat)

  type(mc_lat), intent(inout) :: lat

  integer:: i,j,k
  integer:: row, col, lst, row_new, col_new, lst_new
  integer, dimension(:,:),allocatable :: pos_nn

  ! Put adsorbate on the lattice
  row = 3
  col = 4
  lst = lat%lst(row, col)
  lat%occupations(row,col) = 1
  lat%ads_list(1)%id  = 1
  lat%ads_list(1)%row = row
  lat%ads_list(1)%col = col
  lat%ads_list(1)%ast = 1
  lat%n_ads(1) = lat%n_ads(1) + 1

  allocate(pos_nn(lat%n_nn(lst,1)+1,2))

  pos_nn(1,:) = [row,col]
  do i=1, lat%n_nn(lst,1)
    call lat%neighbor(1,i,row,col)
    pos_nn(1+i,:) = [row,col]
  end do

  print *, 'nn list:'
  do i=1, size(pos_nn,1)
    write(6,'(2i4)') pos_nn(i,:)
  end do
  print*

  do i=1, lat%n_nn(lst,1)

    ! Put another adsorbate into a new site
    lat%occupations(pos_nn(i+1,1),pos_nn(i+1,2)) = i+1
    lat%ads_list(i+1)%id  = 1
    lat%ads_list(i+1)%row = pos_nn(i+1,1)
    lat%ads_list(i+1)%col = pos_nn(i+1,2)
    lat%ads_list(i+1)%ast = 1
    lat%n_ads(1) = lat%n_ads(1) + 1

    lst_new = lat%lst(pos_nn(i+1,1),pos_nn(i+1,2))
    write(6,'(A3,i3,A5,i3)') 'nn:',i, ' lst:',lst_new

!    write(6,'(A)') 'dublicates:'
    do j=1,lat%n_nn(lst_new,1)

      call lat%neighbor(i+1,j,row_new,col_new)
!      do k=1,lat%n_nn(lst,1)+1
!        if (pos_nn(k,1)==row_new .and. pos_nn(k,2)==col_new) then
          write(6,'(3i4)') j,row_new,col_new
!        end if
!      end do

    end do
    write(6,'(A,10i4)') 'nn_new list: ',(lat%nn_new(lst,i,k), k=1,lat%n_nn_new(lst,i))
    print *

  end do


  print *
  call lat%print_ocs()
  call lat%print_st()

  deallocate(pos_nn)

end subroutine

end module
