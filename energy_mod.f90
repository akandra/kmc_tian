module energy_mod

  use constants
  use mc_lat_class
  use energy_parameters_class
  use control_parameters_class

  implicit none

  private
  public  :: energy
  public  :: total_energy

contains

real(dp) function energy(ref, lat, c_pars, e_pars)

  integer, intent(in)      :: ref
  type(mc_lat), intent(in) :: lat
  type(control_parameters), intent(in) :: c_pars
  type(energy_parameters ), intent(in) :: e_pars

  integer :: i, species, n
  integer :: ref_row, ref_col, ref_ast, ref_species, ref_lst
  integer :: row, col, lst, ast, i_ast
  integer, dimension(2) :: direction
  real(8) :: energy_acc

  ref_row      = lat%ads_list(ref)%row
  ref_col      = lat%ads_list(ref)%col
  ref_ast      = lat%ads_list(ref)%ast
  ref_species  = lat%ads_list(ref)%id
  ref_lst      = lat%lst(ref_row,ref_col)

  energy_acc = e_pars%ads_energy(ref_species, ref_lst, ref_ast)

  if (e_pars%is_interaction) then

    do species=1,c_pars%n_species
    do lst=1,n_max_lat_site_types
    do i_ast=1,size(lat%avail_ads_sites(species,lst)%list)
      
      ast = this%avail_ads_sites(species,lst)%list(i_ast)
      do n =1,e_pars%n_interactions(ref_species, ref_lst, ref_ast, species, lst, ast)

        direction = e_pars%neighbors(ref_species, ref_lst, ref_ast, &
                                        species,     lst,     ast, n, :)
        call lat%neighbor2(ref, direction, row, col)

        if (lat%occupations(row,col) > 0) then

          energy_acc = energy_acc + &
            e_pars%int_energy_pars(ref_species, ref_lst, ref_ast,&
                                        species,     lst,     ast,  n)

        end if

      end do
    end do
    end do
    end do

  end if

  energy = energy_acc

end function

real(dp) function total_energy(lat, c_pars, e_pars)

  type(mc_lat), intent(in) :: lat
  type(control_parameters), intent(in) :: c_pars
  type(energy_parameters ), intent(in) :: e_pars

  integer :: i, species, n, ref
  integer :: ref_row, ref_col, ref_ast, ref_species, ref_lst
  integer :: row, col, lst, ast, i_ast
  integer, dimension(2) :: direction

  real(8) :: energy_acc

  energy_acc = 0
  do ref=1,lat%n_ads_tot()

    ref_row      = lat%ads_list(ref)%row
    ref_col      = lat%ads_list(ref)%col
    ref_ast      = lat%ads_list(ref)%ast
    ref_species  = lat%ads_list(ref)%id
    ref_lst      = lat%lst(ref_row,ref_col)

    energy_acc = energy_acc + e_pars%ads_energy(ref_species, ref_lst, ref_ast)

    if (e_pars%is_interaction) then

      do species=1,c_pars%n_species
      do lst=1,n_max_lat_site_types
      do i_ast=1,size(lat%avail_ads_sites(species,lst)%list)
        ast = this%avail_ads_sites(species,lst)%list(i_ast)
        do n =1,e_pars%n_interactions(ref_species, ref_lst, ref_ast, species, lst, ast)

          ast = this%avail_ads_sites(species,lst)%list(i_ast)
          direction = e_pars%neighbors(ref_species, ref_lst, ref_ast, &
                                          species,     lst,     ast, n, :)
          call lat%neighbor2(ref, direction, row, col)

          if (lat%occupations(row,col) > 0) then

            if (lat%occupations(row,col) < ref) cycle ! Avoid double counting
            energy_acc = energy_acc + &
              e_pars%int_energy_pars(ref_species, ref_lst, ref_ast,&
                                        species,     lst,     ast,  n)

          end if

        end do
      end do
      end do
      end do

    end if

  end do

  total_energy = energy_acc

end function


end module energy_mod
