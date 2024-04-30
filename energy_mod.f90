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
  integer :: row, col, id, lst, ast
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
      lst = lat%lst(row,col)

      if (i_n > 0) then

        id  = lat%ads_list(i_n)%id
        ast = lat%ads_list(i_n)%ast
        if (e_pars%int_energy_skip(ref_id, ref_lst, ref_ast, id, lst, ast, i_shell)) cycle

        if (e_pars%int_energy_law_id(ref_id, ref_lst, ref_ast, id, lst, ast) /= linear_id) then
          print'(A,i3,A)', "module energy.f90: function energy: interaction law with id = ", &
                  e_pars%int_energy_law_id(ref_id, ref_lst, ref_ast, id, lst, ast), ' is not yet implemented.'
          stop

        end if

        energy_acc = energy_acc + e_pars%int_energy_pars(ref_id, ref_lst, ref_ast, id, lst, ast, i_shell)

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
  integer :: row, col, id, lst, ast
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
        lst = lat%lst(row,col)

        if (i_n > 0) then

          if (i_n < ref) cycle ! Avoid double counting
          id  = lat%ads_list(i_n)%id
          ast = lat%ads_list(i_n)%ast
          if (e_pars%int_energy_skip(ref_id, ref_lst, ref_ast, id, lst, ast,i_shell)) cycle

          if (e_pars%int_energy_law_id(ref_id, ref_lst, ref_ast, id, lst, ast) /= linear_id) then
            print'(A,i3,A)', "module energy.f90: function energy: interaction law with id = ", &
                    e_pars%int_energy_law_id(ref_id, ref_lst, ref_ast, id, lst, ast), ' is not yet implemented.'
            stop

          end if
          energy_acc = energy_acc + e_pars%int_energy_pars(ref_id, ref_lst, ref_ast, id, lst, ast, i_shell)

        end if

      end do
      end do

    end if

  end do

  total_energy = energy_acc

end function


end module energy_mod
