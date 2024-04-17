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


  call check_nn_new_list(lat,c_pars)

end subroutine

subroutine check_nn_new_list(lat,c_pars)

  type(mc_lat), intent(inout) :: lat
  type(control_parameters), intent(in) :: c_pars

  integer:: i,j,k, ads
  integer:: row, col, lst, row_new, col_new, lst_new
  integer, dimension(:,:),allocatable :: pos_nn
  logical, dimension(:), allocatable :: nn_mask
  logical :: checked = .true.

  write(6,*)

  do ads=1,c_pars%step_period

    row = 3
    col = ads
    lst = lat%lst(row, col)
    lat%occupations(row,col) = 1
    lat%ads_list(1)%id  = 1
    lat%ads_list(1)%row = row
    lat%ads_list(1)%col = col
    lat%ads_list(1)%ast = 1
    lat%n_ads(1) = 1

    write(6,'(A)',advance='no') '  Checking nn_new lists for ' // lat_site_names(lst)

    allocate(pos_nn(lat%n_nn(lst,1)+1,2))

    pos_nn(1,:) = [row,col]
    do i=1, lat%n_nn(lst,1)
      call lat%neighbor(1,i,row,col)
      pos_nn(1+i,:) = [row,col]
    end do

    do i=1, lat%n_nn(lst,1)

      ! Put another adsorbate into a new site
      lat%occupations(pos_nn(i+1,1),pos_nn(i+1,2)) = i+1
      lat%ads_list(i+1)%id  = 1
      lat%ads_list(i+1)%row = pos_nn(i+1,1)
      lat%ads_list(i+1)%col = pos_nn(i+1,2)
      lat%ads_list(i+1)%ast = 1
      lat%n_ads(1) = lat%n_ads(1) + 1

      lst_new = lat%lst(pos_nn(i+1,1),pos_nn(i+1,2))
      allocate(nn_mask(lat%n_nn(lst_new,1)))
      nn_mask = .true.

      do j=1,lat%n_nn(lst_new,1)

        call lat%neighbor(i+1,j,row_new,col_new)
        do k=1,lat%n_nn(lst,1)+1
          if (pos_nn(k,1)==row_new .and. pos_nn(k,2)==col_new) then
            nn_mask(j) = .false.
          end if
        end do

      end do



      do k=1,lat%n_nn_new(lst,i)
        if (.not.nn_mask( lat%nn_new(lst,i,k) )) checked = .false.
      end do

      if (.not.checked) then

        write(6,*)
        write(6,'(A,2i3)') '   adsorbate position: ', pos_nn(1,:)
        write(6,*) '  nn list:'
        do j=2, size(pos_nn,1)
          write(6,'(3i4)') j-1, pos_nn(j,:)
        end do

        write(6,'(A,i2,A)') '   nn ',i, ' lst ' // lat_site_names(lst_new)

        do j=1,lat%n_nn(lst_new,1)
          call lat%neighbor(i+1,j,row_new,col_new)
          write(6,'(3i4,l3)') j,row_new,col_new, nn_mask(j)
        end do

        write(6,'(A)',advance='no') '  mask   list: '
        do k=1,lat%n_nn(lst_new,1)
          if (nn_mask(k)) write(6,'(i4)',advance='no') k
        end do
        write(6,*)
        write(6,'(A,10i4)') '  nn_new list: ',(lat%nn_new(lst,i,k), k=1,lat%n_nn_new(lst,i))

      end if

      deallocate(nn_mask)

    end do

    if (checked) then
      write(6,'(A)') ' correct.'
    else
      checked = .true.
      write(6,*)
    end if

    deallocate(pos_nn)

  end do

end subroutine

end module
