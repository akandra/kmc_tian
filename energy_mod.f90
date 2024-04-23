module energy_mod

  use constants
  use mc_lat_class
  use energy_parameters_class

  implicit none

  private
  public  :: energy
  public  :: total_energy

contains

real(dp) function energy(ref, lat, e_pars)

  integer, intent(in)      :: ref
  type(mc_lat), intent(in) :: lat
  type(energy_parameters ), intent(in) :: e_pars

  integer :: i, i_shell, i_n
  integer :: ref_row, ref_col, ref_ast, ref_id, ref_lst
  integer :: row, col, id
  real(8) :: energy_acc

  ref_row = lat%ads_list(ref)%row
  ref_col = lat%ads_list(ref)%col
  ref_ast = lat%ads_list(ref)%ast
  ref_id  = lat%ads_list(ref)%id
  ref_lst = lat%lst(ref_row,ref_col)

  energy_acc = e_pars%ads_energy(ref_id, ref_lst, ref_ast)

  if (e_pars%is_interaction) then

    do i_shell=1,n_shells
    do i=1, lat%n_nn(ref_lst,i_shell) ! Sum up the int. energy over all nearest neighbors

      call lat%neighbor(ref, i, row, col, i_shell)
      i_n = lat%occupations(row,col)

      if (i_n > 0) then

        id = lat%ads_list(i_n)%id
        if (e_pars%int_energy_skip(ref_id, ref_lst, id, lat%lst(row,col), i_shell)) cycle

        if (e_pars%int_energy_law_id(ref_id, ref_lst, id, lat%lst(row,col)) /= linear_id) then
          print'(A,i3,A)', "module energy.f90: function energy: interaction law with id = ", &
                  e_pars%int_energy_law_id(ref_id, ref_lst, id, lat%lst(row,col)), ' is not yet implemented.'
          stop

        end if

        energy_acc = energy_acc + e_pars%int_energy_pars(ref_id, ref_lst, id, lat%lst(row,col),i_shell)

      end if

    end do
    end do

  end if

  energy = energy_acc

end function

real(dp) function total_energy(lat, e_pars)

  type(mc_lat), intent(in) :: lat
  type(energy_parameters ), intent(in) :: e_pars

  integer :: i, i_shell, i_n, ref
  integer :: ref_row, ref_col, ref_ast, ref_lst, ref_id
  integer :: row, col, id
  real(8) :: energy_acc

  energy_acc = 0
  do ref=1,lat%n_ads_tot()

    ref_row = lat%ads_list(ref)%row
    ref_col = lat%ads_list(ref)%col
    ref_ast = lat%ads_list(ref)%ast
    ref_id  = lat%ads_list(ref)%id
    ref_lst = lat%lst(ref_row,ref_col)

    energy_acc = energy_acc + e_pars%ads_energy(ref_id, ref_lst, ref_ast)

    if (e_pars%is_interaction) then

      do i_shell=1,n_shells
      do i=1, lat%n_nn(ref_lst,i_shell)

        call lat%neighbor(ref, i, row, col, i_shell)
        i_n = lat%occupations(row,col)

        if (i_n > 0) then

          if (i_n < ref) cycle ! Avoid double counting
          id = lat%ads_list(i_n)%id
          if (e_pars%int_energy_skip(ref_id, ref_lst, id, lat%lst(row,col),i_shell)) cycle

          if (e_pars%int_energy_law_id(ref_id, ref_lst, id, lat%lst(row, col)) /= linear_id) then
            print'(A,i3,A)', "module energy.f90: function energy: interaction law with id = ", &
                    e_pars%int_energy_law_id(ref_id, ref_lst, id, lat%lst(row,col)), ' is not yet implemented.'
            stop

          end if
          energy_acc = energy_acc + e_pars%int_energy_pars(ref_id, ref_lst, id, lat%lst(row,col),i_shell)

        end if

      end do
      end do

    end if

  end do

  total_energy = energy_acc

end function

!
!real(8) function total_energy(nlat, nads, nnn, occupations, site_type, &
!                              ads_list, nn_list, ads_energy, int_energy)
!
!integer, intent(in) :: nlat, nads, nnn
!integer, dimension(nlat,nlat), intent(in) :: occupations
!integer(1), dimension(nlat,nlat), intent(in) :: site_type
!integer, dimension(nads,2), intent(in) :: ads_list
!integer, dimension(nnn,2), intent(in) :: nn_list
!real(8), dimension(3), intent(in) :: ads_energy
!real(8), dimension(3,3), intent(in) :: int_energy
!
!integer :: i, ic, jc, iads, jads
!real(8) :: energy_acc
!
!    energy_acc = 0.0d0
!    do inx=1,nads
!
!        iads = ads_list(inx,1)
!        jads = ads_list(inx,2)
!
!        energy_acc = energy_acc + ads_energy(site_type(iads,jads))
!
!        do i=1, nnn/2 ! count interaction with neighbors once
!
!            ic = modulo(iads+nn_list(i,1)-1,nlat) + 1
!            jc = modulo(jads+nn_list(i,2)-1,nlat) + 1
!            if (occupations(ic,jc) > 0) &
!                energy_acc = energy_acc +&
!                        int_energy(site_type(iads,jads),site_type(ic,jc))
!
!        end do
!
!    end do
!
!    total_energy = energy_acc
!
!end function

end module energy_mod
