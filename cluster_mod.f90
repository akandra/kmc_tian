module cluster_mod

  use mc_lat_class
!  use energy_parameters_class

  implicit none

contains


integer function lfind(x, labels)

  integer, intent(inout) :: x
  integer, dimension(:), intent(inout) :: labels

  integer :: y, z

  y = x
    do while (labels(y) /= y)
        y = labels(y)
    end do

!    do while (labels(x) /= x)
!        z = labels(x)
!        labels(x) = y
!        x = z
!    end do
!
    lfind = y

end function

subroutine lunion(x, y, labels)

  integer, intent(inout) :: x, y
  integer, dimension(:), intent(inout) :: labels

    if (x > y) then
        labels( lfind(x, labels) ) = lfind(y, labels)
    else if ( x < y) then
        labels( lfind(y, labels) ) = lfind(x, labels)
    endif

end subroutine

subroutine hoshen_kopelman(lat, species, cluster_label, largest_label)

type(mc_lat), intent(in) :: lat
integer, intent(in) :: species
integer, dimension(:,:), intent(out) :: cluster_label
integer, intent(out) :: largest_label

integer :: i, j, m, n, nnn2, row, col, row_new,col_new
integer :: itemp, ip, jp, i_ads, i_ads_nn

integer, dimension(lat%n_ads(species)) :: labels
integer, dimension(lat%n_nn(1)/2) :: scanned_nn_occs, row_nn, col_nn

nnn2 = lat%n_nn(1)/2
scanned_nn_occs = 0
cluster_label   = 0
do i=1,lat%n_ads(species)
    labels(i) = i
end do

largest_label = 0
call lat%print_ocs

do row=1,lat%n_rows
do col=1,lat%n_cols

  i_ads = lat%occupations(row,col)

  if ( i_ads > 0 .and. lat%ads_list(i_ads)%id == species) then

    do m=1,nnn2

      ! Take previously scan neighbours
      call lat%hop(i_ads,nnn2+m,row_nn(m),col_nn(m))
      i_ads_nn = lat%occupations(row_nn(m), col_nn(m))
      if ( i_ads_nn > 0 .and. lat%ads_list(i_ads_nn)%id == species ) then
          scanned_nn_occs(m) = 1
      else
          scanned_nn_occs(m) = 0
      end if

      ! Switch off periodic boundary conditions
      ! Warning: works only for hexagonal lattice,
      !          maybe, because of Corona virus
      if (row==1) scanned_nn_occs(2:3) = 0
      if (col==1) scanned_nn_occs(1) = 0
      if (col==lat%n_cols) scanned_nn_occs(3) = 0

    end do

    select case (sum(scanned_nn_occs))

      case (0)
        !print*
        !print*,"Case 0"
        largest_label = largest_label + 1
        cluster_label(row,col) = largest_label

      case (1)
        !print*
        !print*,"Case 1"
        do m=1,nnn2
          if (scanned_nn_occs(m) == 1) then
            cluster_label(row,col) = &
                    lfind( cluster_label(row_nn(m),col_nn(m)), labels)
            end if
        end do

      case (2)
        !print*
        !print*,"Case 2"
        do m=1,nnn2-1
          itemp = cluster_label(row_nn(m),col_nn(m))
        do n=m+1,nnn2
          if (scanned_nn_occs(m) == 1 .and. scanned_nn_occs(n) == 1) then
            call lunion(itemp,cluster_label(row_nn(n),col_nn(n)),labels)
            cluster_label(row,col) = lfind(itemp, labels)
          end if
        end do
        end do

      case (3)
        !print*
        !print*,"Case 3"
        itemp = cluster_label(row_nn(1),col_nn(1))
        call lunion(itemp, cluster_label(row_nn(2),col_nn(2)),labels)
        call lunion(itemp, cluster_label(row_nn(3),col_nn(3)),labels)
        cluster_label(row,col) = lfind(itemp, labels)

      case default
        stop 'Hoshen-Kopelman: too many neighbors'

    end select

  end if

end do
end do

print*
call lat%print_ocs
print*
write(*,'(8i4)') transpose(cluster_label)
print*
do i=1,largest_label
    print('(i3,a5,i3)'), i, ' --> ' , labels(i)
end do

! Apply PBC
do row=1,lat%n_rows
  ip = cluster_label(row,1)
  if (ip > 0) then
    do m=1,nnn2
      call lat%hop(lat%occupations(row,1),nnn2+m,row_new,col_new)
      jp = cluster_label(row_new,col_new)
      if (jp > 0) call lunion(ip,jp,labels)
    end do
  end if
end do

do col=2,lat%n_cols-1
  ip = cluster_label(1,col)
  if (ip > 0) then
    do m=1,nnn2
      call lat%hop(lat%occupations(1,col),nnn2+m,row_new,col_new)
      jp = cluster_label(row_new,col_new)
      if (jp > 0) call lunion(ip,jp,labels)
    end do
  end if
end do

print*
do i=1,largest_label
    print('(i3,a5,i3)'), i, ' --> ' , labels(i)
end do

! Going down to roots
do row=1,lat%n_rows
do col=1,lat%n_cols
  if (cluster_label(row,col) > 0) &
     cluster_label(row,col) = lfind(cluster_label(row,col), labels)
end do
end do

print*
write(*,'(8i4)') transpose(cluster_label)

end subroutine


!subroutine count_cluster_sizes(cluster_sizes, cluster_label, ads_list, nads, nlat)
!
!integer, dimension(nads), intent(out) :: cluster_sizes
!integer, intent(in) :: nads, nlat
!integer, dimension(nads,2), intent(in) :: ads_list
!integer, dimension(nlat,nlat), intent(in) :: cluster_label
!
!
!integer :: i, ic
!
!    cluster_sizes   = 0
!
!   do i=1, nads
!        ic = cluster_label(ads_list(i,1),ads_list(i,2))
!        cluster_sizes(ic) = cluster_sizes(ic) + 1
!    end do
!
!end subroutine
!

end module cluster_mod
