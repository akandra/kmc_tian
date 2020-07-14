module reaction_class
  !---------------------------------------------------------------------------
  !  Module for binary dissocation reaction  r --> p1 + p2
  !---------------------------------------------------------------------------
  !  * controlled by the dissociation keyword section of the .reaction file
  !  * p1 is assumed to be at the same lattice site as r with its adsorption
  !    site type (ast) choosen from those specified in the .reaction file
  !  * p2 lattice site is taken as a nearest neighbor of r with its ast
  !    taken from those specified in the .reaction file
  !---------------------------------------------------------------------------
  use constants
  use control_parameters_class
  use mc_lat_class
  use energy_parameters_class
  use energy_mod
  use open_file
  use utilities
  use rates_hopping_class
  use rates_desorption_class
  use rates_dissociation_class

  implicit none

  private
  public    :: reaction_init, reaction_type


  type :: reaction_type

    type(     hopping_type) :: hopping
    type(  desorption_type) :: desorption
    type(dissociation_type) :: dissociation

    real(dp) :: beta                    ! inverse thermodynamic temperature
    integer  :: n_ads_total             ! total number of adsorbates
    real(dp) :: total_rate              ! sum over all relevant rates

  contains
    procedure :: construct
    procedure :: do_reaction

  end type


contains
!-------------------------------------1-----------------------------------------
  function reaction_init(c_pars, lat, e_pars)
!------------------------------------------------------------------------------
    type(reaction_type) reaction_init

    type(control_parameters), intent(in)    :: c_pars
    type(mc_lat)            , intent(in)    :: lat
    type(energy_parameters) , intent(in)    :: e_pars


    ! Create a rate structure
    reaction_init%hopping      = hopping_init(c_pars, lat, e_pars)
    reaction_init%desorption   = desorption_init(c_pars, lat, e_pars)
    reaction_init%dissociation = dissociation_init(c_pars, lat, e_pars)
!    call reaction_init%hopping%print(c_pars)
!    call reaction_init%desorption%print(c_pars)
!    call reaction_init%dissociation%print(c_pars)

    reaction_init%beta = 1.0_dp/(kB*c_pars%temperature)
    reaction_init%n_ads_total = lat%n_ads_tot()

  end function reaction_init

!-----------------------------------------------------------------------------
  subroutine construct(this, lat, e_pars)
!-----------------------------------------------------------------------------
    class(reaction_type),     intent(inout) :: this
    class(mc_lat),            intent(inout) :: lat
    class(energy_parameters), intent(in)    :: e_pars

    integer       :: i
    integer       :: ads, chan, proc
    integer       :: r_id, r_lst, r_ast
    integer       :: p1_id, p1_lst, p1_ast, p2_id, p2_lst, p2_ast
    real(dp)      :: rate
    integer       :: m
    character(8)  :: lst_name, lst_p2_name, ast_r_name, ast_p1_name, ast_p2_name

    ! Hopping
    if (this%hopping%is_defined) then
      do ads = 1, this%n_ads_total
        call this%hopping%construct(ads, lat, e_pars, this%beta)
      end do
    end if

    ! Desorption
    if (this%desorption%is_defined) then
      do ads = 1, this%n_ads_total
        call this%desorption%construct(ads, lat, e_pars, this%beta)
      end do
    end if

    ! Dissociation
    if (this%dissociation%is_defined) then

      do ads = 1, this%n_ads_total
        call this%dissociation%construct(ads, lat, e_pars, this%beta)
      end do

    end if

! Debugging printout desorption
print*
call lat%print_ocs
print *
print '(a)' ,' Dissociation channels from reaction class'
print '(a)' ,' --------------------------------------------------------------------'
print '(a)', '  ads# proc#  m  rate       lst      ast_r   ast_p1   lst_p2   ast_p2'
print '(a)' ,' --------------------------------------------------------------------'

do ads = 1, this%n_ads_total
  do chan= 1, this%dissociation%rate_info(ads)%n_channels

    i = this%dissociation%rate_info(ads)%list(chan)%proc

    r_id  = this%dissociation%channels(i)%r
    r_ast = this%dissociation%channels(i)%r_ast
    r_lst = this%dissociation%channels(i)%r_lst

    p1_id  = this%dissociation%channels(i)%p1
    p1_ast = this%dissociation%channels(i)%p1_ast
    p1_lst = this%dissociation%channels(i)%p1_lst

    p2_id  = this%dissociation%channels(i)%p2
    p2_ast = this%dissociation%channels(i)%p2_ast
    p2_lst = this%dissociation%channels(i)%p2_lst

    rate   = this%dissociation%channels(i)%rate

    m      = this%dissociation%rate_info(ads)%list(chan)%m

    lst_name    = site_names(r_lst)
    lst_p2_name = site_names(p2_lst)
    ast_r_name  = ads_site_names(r_ast)
    ast_p1_name = ads_site_names(p1_ast)
    ast_p2_name = ads_site_names(p2_ast)

    print '(t4,i0, t9,i0, t15,i0, t16,1pe10.2, t29,A7, 2x,A3, t46,A3, t55,A7, 2x,A3)', &
            ads, i, m, rate, lst_name, ast_r_name, ast_p1_name, lst_p2_name, ast_p2_name

  end do
  print*
end do

  end subroutine construct

!-----------------------------------------------------------------------------
  subroutine do_reaction(this, rand, lat, e_pars)
!-----------------------------------------------------------------------------
    class(reaction_type), intent(inout)     :: this
    real(dp),             intent(in)        :: rand
    class(mc_lat),        intent(inout)     :: lat
    class(energy_parameters), intent(in)    :: e_pars

    integer :: i, m, ads, iads, reaction_id, i_change, channel
    integer :: n_nn, n_nn2, m_nn, col_p2, row_p2, proc
    integer :: id_r, id_p1, id_p2, ast_r, ast_p1, ast_p2
    integer :: row, col, row_new, col_new, lst_new, ast_new, id
    real(dp), dimension(n_reaction_types) :: acc_rate
    integer, dimension(2*lat%n_nn(1)) :: change_list

    real(dp) :: u, temp_dp

    n_nn  = lat%n_nn(1)
    n_nn2 = n_nn/2

    ! -------- calculate total rate

    ! initialize the accumulated rate for hopping
    acc_rate(hopping_id) = 0.0_dp
    ! rate for hopping reactions
    if (this%hopping%is_defined) then
      do ads=1,this%n_ads_total
      do m=1,n_nn
        acc_rate(hopping_id) = acc_rate(hopping_id) + sum(this%hopping%rates(ads,m)%list)
      end do
      end do
    end if
    ! accumulate rate for desorption (= rate for desorption reactions + hopping reactions)
    acc_rate(desorption_id) = acc_rate(hopping_id)
    if (this%desorption%is_defined) then
      do ads=1,this%n_ads_total
        acc_rate(desorption_id) = acc_rate(desorption_id) + this%desorption%rates(ads)
      end do
    end if

    acc_rate(dissociation_id) = acc_rate(desorption_id)
    if (this%dissociation%is_defined) then
      do ads=1,this%n_ads_total
      do channel = 1, this%dissociation%rate_info(ads)%n_channels
        acc_rate(dissociation_id) = acc_rate(dissociation_id) &
                                  + this%dissociation%rate_info(ads)%list(channel)%rate
      end do
      end do
    end if

    ! Save the total rate value
    this%total_rate = acc_rate(n_reaction_types)
    ! do nothing if there is nothing to do
    if (this%total_rate == 0.0_dp) return

    ! ------- Select the process

    ! random rate value to select a process
    u = rand*this%total_rate
    ! determine the type of reaction
    do reaction_id=1,n_reaction_types
      if (u < acc_rate(reaction_id)) exit
    end do

    select case (reaction_id)

      case(hopping_id)
        ! determine hopping channel (adsorbate, direction, ads. site of available ones)
        !                           (      ads,      m_nn,                        iads)
        temp_dp = 0.0_dp
        extloop: do ads=1,this%n_ads_total
          do m_nn=1,n_nn
          do iads=1,size(this%hopping%rates(ads,m_nn)%list) ! Warning: check timing of size calculation
            temp_dp = temp_dp + this%hopping%rates(ads,m_nn)%list(iads)
            if (u < temp_dp) exit extloop
          end do
          end do
        end do extloop

        ! Update hopping rate constants
        ! Particle (ads) is going to hop in direction (m_nn) to ads. site (iads)

        ! create a list of adsorbates affected by hop

        change_list = 0
        i_change = 1
        ! Put the hopping particle iads into the list
        change_list(i_change) = ads

        ! scan over old neighbors
        do m=1,n_nn
          ! position of neighbor m
          call lat%neighbor(ads,m,row,col)
          if (lat%occupations(row,col) > 0) then
            i_change = i_change + 1
            change_list(i_change) = lat%occupations(row,col)
          end if
        end do

        ! a new position of particle (ads) after a hop to a neighbor (m_nn)
        call lat%neighbor(ads,m_nn,row_new,col_new)
        id = lat%ads_list(ads)%id
        lst_new = lat%lst(row_new,col_new)
        ast_new = lat%avail_ads_sites(id,lst_new)%list(iads)

        ! Make a hop:
        ! Delete an adsorbate from its old position
        lat%occupations(lat%ads_list(ads)%row, lat%ads_list(ads)%col) = 0
        ! Put the adsorbate in a new position
        lat%ads_list(ads)%row  = row_new
        lat%ads_list(ads)%col  = col_new
        lat%ads_list(ads)%ast  = ast_new
        ! Update adsorbate position
        lat%occupations(row_new, col_new) = ads

        ! scan over additional new neighbors
        do m=1,n_nn2
          ! position of neighbor nn_new(m_nn,m)
          call lat%neighbor( ads, this%hopping%nn_new(m_nn,m), row, col)
          if (lat%occupations(row,col) > 0) then
            i_change = i_change + 1
            change_list(i_change) = lat%occupations(row,col)
          end if
        end do

        ! Update rate array for the affected adsorbates

        do i=1,i_change
          call      this%hopping%construct(change_list(i), lat, e_pars, this%beta)
          call   this%desorption%construct(change_list(i), lat, e_pars, this%beta)
          call this%dissociation%construct(change_list(i), lat, e_pars, this%beta)
        end do

      case(desorption_id)
        ! determine desorption channel (adsorbate)
        !                              (      ads)
        temp_dp = acc_rate(hopping_id)
        do ads=1,this%n_ads_total
          temp_dp = temp_dp + this%desorption%rates(ads)
          if (u < temp_dp) exit
        end do

        ! Update desorption rate constants
        ! Particle (ads) is going to desorb to nowhere

        ! create a list of adsorbates affected by desorption
        change_list = 0
        i_change = 0

        ! scan over neighbors
        do m=1,n_nn
          ! position of neighbor m
          call lat%neighbor(ads,m,row,col)
          if (lat%occupations(row,col) > 0) then
            i_change = i_change + 1
            change_list(i_change) = lat%occupations(row,col)
          end if
        end do

        ! Do desorption:
        ! Delete an adsorbate from the lattice
        lat%occupations(lat%ads_list(ads)%row, lat%ads_list(ads)%col) = 0
        ! Update the number of adsorbates in the lat structure
        lat%n_ads(lat%ads_list(ads)%id) = lat%n_ads(lat%ads_list(ads)%id) - 1
        ! Rearrange adsorbates except when the last adsorbate desorbs
        if (ads < this%n_ads_total) then
          ! Put the last adsorbate in place of ads
          lat%ads_list(ads) = lat%ads_list(this%n_ads_total)
          ! Update the adsorbate number in the lattice
          lat%occupations(lat%ads_list(ads)%row, lat%ads_list(ads)%col) = ads
          ! Update the rates arrays
          this%hopping%rates(ads,:)        = this%hopping%rates(this%n_ads_total,:)
          this%desorption%rates(ads)       = this%desorption%rates(this%n_ads_total)
          this%dissociation%rate_info(ads) = this%dissociation%rate_info(this%n_ads_total)
        end if

        ! Update rate array for the affected adsorbates
        do i=1,i_change
          ! account for the tightening the ads. list
          if (change_list(i) == this%n_ads_total) change_list(i) = ads
          call      this%hopping%construct(change_list(i), lat, e_pars, this%beta)
          call   this%desorption%construct(change_list(i), lat, e_pars, this%beta)
          call this%dissociation%construct(change_list(i), lat, e_pars, this%beta)
        end do

        ! Clean up rate arrays for the last adsorbate
        do m=1,n_nn
          this%hopping%rates(this%n_ads_total,m)%list = 0.0_dp
        end do
        this%desorption%rates(this%n_ads_total) = 0.0_dp
        this%dissociation%rate_info(this%n_ads_total)%list = rate_info(0,0,0.0_dp)
        this%dissociation%rate_info(this%n_ads_total)%n_channels = 0

        ! Update the the number of adsorbates in this
        this%n_ads_total = this%n_ads_total - 1

      case(dissociation_id)

print*
print*, 'Before dissociation'
call lat%print_ocs
call lat%print_ads

        ! determine dissociation channel (ads, channel)
        temp_dp = acc_rate(desorption_id)
        extloop2: do ads=1,this%n_ads_total
        do channel=1,this%dissociation%rate_info(ads)%n_channels
          temp_dp = temp_dp + this%dissociation%rate_info(ads)%list(channel)%rate
          if (u < temp_dp) exit extloop2
        end do
        end do extloop2

        ! create a list of adsorbates affected by dissociation

        change_list = 0
        i_change = 1
        ! Put the reactant's (which to become p1) ads  into the list
        change_list(i_change) = ads

        ! scan over neighbors of reactant
        do m=1,n_nn
          ! position of neighbor m
          call lat%neighbor(ads,m,row,col)
          if (lat%occupations(row,col) > 0) then
            i_change = i_change + 1
            change_list(i_change) = lat%occupations(row,col)
          end if
        end do

        ! Dissociate:
        ! Find process number
        proc = this%dissociation%rate_info(ads)%list(channel)%proc
        ! Get reactant and products info
        id_r   = lat%ads_list(ads)%id
        id_p1  = this%dissociation%channels(proc)%p1
        ast_p1 = this%dissociation%channels(proc)%p1_ast
        id_p2  = this%dissociation%channels(proc)%p2
        ast_p2 = this%dissociation%channels(proc)%p2_ast
        ! Replace reactant with product1 preserving its number
        lat%ads_list(ads)%id  = id_p1
        lat%ads_list(ads)%ast = ast_p1
        lat%n_ads(id_r)  = lat%n_ads(id_r)  - 1
        lat%n_ads(id_p1) = lat%n_ads(id_p1) + 1
        ! direction for dissociation
        m_nn = this%dissociation%rate_info(ads)%list(channel)%m
        ! position for product 2
        call lat%neighbor(ads, m_nn ,row_p2,col_p2)
        ! Add product 2
        ! increase the total number of adsorbates
        this%n_ads_total = this%n_ads_total  + 1
        ! update occupations
        lat%occupations(row_p2,col_p2) = this%n_ads_total
        ! update adsorbate list
        lat%ads_list(this%n_ads_total)%id  =  id_p2
        lat%ads_list(this%n_ads_total)%row = row_p2
        lat%ads_list(this%n_ads_total)%col = col_p2
        lat%ads_list(this%n_ads_total)%ast = ast_p2
        lat%n_ads(id_p2) = lat%n_ads(id_p2) + 1

        ! scan over neighbors of product 2
        do m=1,n_nn2
          ! position of neighbor with direction nn_new(m_nn,m)
          call lat%neighbor( this%n_ads_total, this%dissociation%nn_new(m_nn,m), row, col)
          if (lat%occupations(row,col) > 0) then
            i_change = i_change + 1
            change_list(i_change) = lat%occupations(row,col)
          end if
        end do
        ! add product 2 to change_list
        i_change = i_change + 1
        change_list(i_change) = lat%occupations(row_p2,col_p2)

print*, 'After dissociation'
call lat%print_ocs
print*, 'Change list:'
print*, 'Number of adsorbates which feel dissociations:', i_change
print*, change_list
call lat%print_ads
pause

        ! Reset the rate info for ads (product 1)
        do m=1,n_nn
          this%hopping%rates(ads,m)%list = 0.0_dp
        end do
        this%desorption%rates(ads) = 0.0_dp
        this%dissociation%rate_info(ads)%list = rate_info(0,0,0.0_dp)
        this%dissociation%rate_info(ads)%n_channels = 0
        ! Update rate array for the affected adsorbates
        do i=1,i_change
          call      this%hopping%construct(change_list(i), lat, e_pars, this%beta)
          call   this%desorption%construct(change_list(i), lat, e_pars, this%beta)
          call this%dissociation%construct(change_list(i), lat, e_pars, this%beta)
        end do

      case default
        print*
        print*, "reaction id is ", reaction_id
        print*, "number of reactions is ", n_reaction_types
        stop 'do_reaction (of reaction_class): must never occur!'

    end select

  end subroutine do_reaction

end module reaction_class
