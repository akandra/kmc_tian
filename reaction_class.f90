module reaction_class
  !---------------------------------------------------------------------------
  !  Module for dealing with reactiona
  !---------------------------------------------------------------------------
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
  use rates_association_class
  use rates_bimolecular_class

  implicit none

  private
  public    :: reaction_init, reaction_type


  type :: reaction_type

    type(     hopping_type) :: hopping
    type(  desorption_type) :: desorption
    type(dissociation_type) :: dissociation
    type( association_type) :: association
    type( bimolecular_type) :: bimolecular

    real(dp) :: beta                                    ! inverse thermodynamic temperature
    integer  :: n_ads_total                             ! total number of adsorbates
    real(dp), dimension(n_reaction_types) :: acc_rate   ! accumulated rate for reaction types
    ! reaction counters
    !                  reaction
    !                  . species
    !                  . .
    integer, dimension(:,:), allocatable :: counter

  contains
    procedure :: accumulate_rates
    procedure :: construct
    procedure :: construct_1
    procedure :: cleanup_rates
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
    reaction_init%hopping      =      hopping_init(c_pars, lat, e_pars)
    reaction_init%desorption   =   desorption_init(c_pars, lat, e_pars)
    reaction_init%dissociation = dissociation_init(c_pars, lat, e_pars)
    reaction_init%association  =  association_init(c_pars, lat, e_pars)
    reaction_init%bimolecular  =  bimolecular_init(c_pars, lat, e_pars)

!    call reaction_init%association%print(c_pars)

    reaction_init%beta = 1.0_dp/(kB*c_pars%temperature)
    reaction_init%n_ads_total = lat%n_ads_tot()

    allocate(reaction_init%counter(n_reaction_types,c_pars%n_species))
    reaction_init%counter = 0

  end function reaction_init

!-----------------------------------------------------------------------------
  subroutine accumulate_rates(this, lat)
!-----------------------------------------------------------------------------
    class(reaction_type),     intent(inout) :: this
    class(mc_lat),            intent(inout) :: lat

!dja    integer :: i, m, ads, channel, n_nn
    integer :: m, ads, channel, n_nn

    ! initialize the accumulated rate for hopping
    this%acc_rate(hopping_id) = 0.0_dp
    ! rate for hopping reactions
    if (this%hopping%is_defined) then
      do ads=1,this%n_ads_total
        n_nn = lat%n_nn( lat%lst( lat%ads_list(ads)%row,lat%ads_list(ads)%col ), 1)
        do m=1,n_nn+1 ! '+1' accounts for intra-site hops
          this%acc_rate(hopping_id) = this%acc_rate(hopping_id) + sum(this%hopping%rates(ads,m)%list)
        end do
      end do
    end if

    ! accumulate rate for desorption (= rate for desorption + hopping)
    this%acc_rate(desorption_id) = this%acc_rate(hopping_id)
    if (this%desorption%is_defined) then
      do ads=1,this%n_ads_total
        this%acc_rate(desorption_id) = this%acc_rate(desorption_id) + this%desorption%rates(ads)
      end do
    end if

    ! accumulate rate for dissociation (= rate for dissociation + desorption + hopping)
    this%acc_rate(dissociation_id) = this%acc_rate(desorption_id)
    if (this%dissociation%is_defined) then
      do ads=1,this%n_ads_total
      do channel = 1, this%dissociation%rate_info(ads)%n_channels
        this%acc_rate(dissociation_id) = this%acc_rate(dissociation_id) &
                                       + this%dissociation%rate_info(ads)%list(channel)%rate
      end do
      end do
    end if

    ! accumulate rate for association (= rate for association + dissociation + desorption + hopping)
    this%acc_rate(association_id) = this%acc_rate(dissociation_id)
    if (this%association%is_defined) then
      do ads=1,this%n_ads_total
      do channel = 1, this%association%rate_info(ads)%n_channels
        this%acc_rate(association_id) = this%acc_rate(association_id) &
                                      + this%association%rate_info(ads)%list(channel)%rate
      end do
      end do
    end if

    ! accumulate rate for bimolecular reaction (= rate for bimolecular + association + dissociation + desorption + hopping)
    this%acc_rate(bimolecular_id) = this%acc_rate(association_id)
    if (this%bimolecular%is_defined) then
      do ads=1,this%n_ads_total
      do channel = 1, this%bimolecular%rate_info(ads)%n_channels
        this%acc_rate(bimolecular_id) = this%acc_rate(bimolecular_id) &
                                      + this%bimolecular%rate_info(ads)%list(channel)%rate
      end do
      end do
    end if

!    ! Debug printing of accumulated rates
!    print*
!    print*,'accumulate_rates debugging output:'
!    print'(3A,1pe12.2)','Total rate for ',reaction_names(1), ' is ', this%acc_rate(1)
!    do i=2, n_reaction_types
!      print'(3A,1pe12.2)','Total rate for ',reaction_names(i), ' is ', this%acc_rate(i) - this%acc_rate(i-1)
!    end do

  end subroutine accumulate_rates


!-----------------------------------------------------------------------------
  subroutine construct(this, lat, e_pars)
!-----------------------------------------------------------------------------
    class(reaction_type),     intent(inout) :: this
    class(mc_lat),            intent(inout) :: lat
    class(energy_parameters), intent(in)    :: e_pars

    integer       :: ads
    ! Debugging variables
!dja    integer :: chan, i, m
!dja    integer :: r_id, r_ast, r_lst
!dja    integer :: p1_id, p1_ast, p1_lst
!dja    integer :: p2_id, p2_ast, p2_lst
!dja    integer :: r1_id, r1_ast, r1_lst
!dja    integer :: r2_id, r2_ast, r2_lst
!dja    real(dp) :: rate
!dja    character(len=max_string_length) :: lst_name, lst_p1_name, lst_p2_name, lst_r1_name, lst_r2_name
!dja    character(len=max_string_length) :: ast_r_name, ast_r1_name, ast_r2_name, ast_p1_name, ast_p2_name

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
!dja    call this%dissociation%construct(ads, lat, e_pars, this%beta)
        call this%dissociation%construct(ads, lat)
      end do

    end if

    ! Association
    if (this%association%is_defined) then

      do ads = 1, this%n_ads_total
!dja    call this%association%construct(ads, lat, e_pars, this%beta)
        call this%association%construct(ads, lat)
      end do

    end if

    ! Bimolecular
    if (this%bimolecular%is_defined) then

      do ads = 1, this%n_ads_total
!dja    call this%bimolecular%construct(ads, lat, e_pars, this%beta)
        call this%bimolecular%construct(ads, lat)
      end do

    end if

    call this%accumulate_rates(lat)

! Debugging printout dissociation
!print*
!print*, 'construct debugging printout: '
!call lat%print_ocs
!print *
!print '(a)' ,' Dissociation channels from reaction class'
!print '(a)' ,' --------------------------------------------------------------------'
!print '(a)', '  ads# proc#  m  rate       lst      ast_r   ast_p1   lst_p2   ast_p2'
!print '(a)' ,' --------------------------------------------------------------------'
!
!do ads = 1, this%n_ads_total
!  do chan= 1, this%dissociation%rate_info(ads)%n_channels
!
!    i = this%dissociation%rate_info(ads)%list(chan)%proc
!
!    r_id  = this%dissociation%channels(i)%r
!    r_ast = this%dissociation%channels(i)%r_ast
!    r_lst = this%dissociation%channels(i)%r_lst
!
!    p1_id  = this%dissociation%channels(i)%p1
!    p1_ast = this%dissociation%channels(i)%p1_ast
!    p1_lst = this%dissociation%channels(i)%p1_lst
!
!    p2_id  = this%dissociation%channels(i)%p2
!    p2_ast = this%dissociation%channels(i)%p2_ast
!    p2_lst = this%dissociation%channels(i)%p2_lst
!
!    rate   = this%dissociation%channels(i)%rate
!
!    m      = this%dissociation%rate_info(ads)%list(chan)%m
!
!    lst_name    = lat_site_names(r_lst)
!    lst_p2_name = lat_site_names(p2_lst)
!    ast_r_name  = ads_site_names(r_ast)
!    ast_p1_name = ads_site_names(p1_ast)
!    ast_p2_name = ads_site_names(p2_ast)
!
!    print '(t4,i0, t9,i0, t15,i0, t16,1pe10.2, t29,A7, 2x,A3, t46,A3, t55,A7, 2x,A3)', &
!            ads, i, m, rate, lst_name, ast_r_name, ast_p1_name, lst_p2_name, ast_p2_name
!
!  end do
!  if ( this%dissociation%rate_info(ads)%n_channels>0 ) print*
!end do
!
!! Debugging printout association
!print*
!print '(a)' ,' Association channels from reaction class'
!print '(a)'    ,' --------------------------------------------------------------------------------------'
!print '(a)'    ,'  ads# proc#       m   rate       lst_r1   ast_r1    lst_r2   ast_r2   lst_p1   ast_p1'
!print '(a)'    ,' --------------------------------------------------------------------------------------'
!
!do ads = 1, this%n_ads_total
!  do chan= 1, this%association%rate_info(ads)%n_channels
!
!    i = this%association%rate_info(ads)%list(chan)%proc
!
!    r1_id  = this%association%channels(i)%r1
!    r1_ast = this%association%channels(i)%r1_ast
!    r1_lst = this%association%channels(i)%r1_lst
!
!    r2_id  = this%association%channels(i)%r2
!    r2_ast = this%association%channels(i)%r2_ast
!    r2_lst = this%association%channels(i)%r2_lst
!
!    p1_id  = this%association%channels(i)%p1
!    p1_ast = this%association%channels(i)%p1_ast
!    p1_lst = this%association%channels(i)%p1_lst
!
!    rate   = this%association%channels(i)%rate
!
!    m      = this%association%rate_info(ads)%list(chan)%m
!
!    lst_r1_name = lat_site_names(r1_lst)
!    lst_r2_name = lat_site_names(r2_lst)
!    lst_p1_name = lat_site_names(p1_lst)
!    ast_r1_name = ads_site_names(r1_ast)
!    ast_r2_name = ads_site_names(r2_ast)
!    ast_p1_name = ads_site_names(p1_ast)
!
!  print '(t4,i0, t9,i0, t15,i0, t20,1pe10.2, t35,A7, 2x,A3, t54,A7, 2x, A3, t72,A7, 2x,A3)', &
!        ads, i, m, rate, lst_r1_name, ast_r1_name, lst_r2_name, ast_r2_name, lst_p1_name, ast_p1_name
!
!  end do
!  if ( this%dissociation%rate_info(ads)%n_channels>0 ) print*
!end do
!! Debugging printout bimolecular
!print*
!print '(a)' ,' Bimolecular channels from reaction class'
!print '(a)'    ,' -------------------------------------------------------------------------------------------------------'
!print '(a)'    ,'  ads# proc#       m   rate       lst_r1   ast_r1    lst_r2   ast_r2   lst_p1   ast_p1   lst_p2   ast_p2'
!print '(a)'    ,' -------------------------------------------------------------------------------------------------------'
!
!do ads = 1, this%n_ads_total
!  do chan= 1, this%bimolecular%rate_info(ads)%n_channels
!
!    i = this%bimolecular%rate_info(ads)%list(chan)%proc
!
!    r1_id  = this%bimolecular%channels(i)%r1
!    r1_ast = this%bimolecular%channels(i)%r1_ast
!    r1_lst = this%bimolecular%channels(i)%r1_lst
!
!    r2_id  = this%bimolecular%channels(i)%r2
!    r2_ast = this%bimolecular%channels(i)%r2_ast
!    r2_lst = this%bimolecular%channels(i)%r2_lst
!
!    p1_id  = this%bimolecular%channels(i)%p1
!    p1_ast = this%bimolecular%channels(i)%p1_ast
!    p1_lst = this%bimolecular%channels(i)%p1_lst
!
!    p2_id  = this%bimolecular%channels(i)%p2
!    p2_ast = this%bimolecular%channels(i)%p2_ast
!    p2_lst = this%bimolecular%channels(i)%p2_lst
!
!    rate   = this%bimolecular%channels(i)%rate
!
!    m      = this%bimolecular%rate_info(ads)%list(chan)%m
!
!    lst_r1_name = lat_site_names(r1_lst)
!    lst_r2_name = lat_site_names(r2_lst)
!    lst_p1_name = lat_site_names(p1_lst)
!    lst_p2_name = lat_site_names(p2_lst)
!    ast_r1_name = ads_site_names(r1_ast)
!    ast_r2_name = ads_site_names(r2_ast)
!    ast_p1_name = ads_site_names(p1_ast)
!    ast_p2_name = ads_site_names(p2_ast)
!
!  print '(t4,i0, t9,i0, t15,i0, t20,1pe10.2, t35,A7, 2x,A3, t54,A7, 2x, A3, t72,A7, 2x,A3, t90,A7, 2x,A3)', &
!        ads, i, m, rate, lst_r1_name, ast_r1_name, lst_r2_name, ast_r2_name, lst_p1_name, ast_p1_name, lst_p2_name, ast_p2_name
!
!  end do
!  if ( this%association%rate_info(ads)%n_channels>0 ) print*
!end do

  end subroutine construct

!-----------------------------------------------------------------------------
  subroutine construct_1(this, ads, lat, e_pars)
!-----------------------------------------------------------------------------
    class(reaction_type),     intent(inout) :: this
    integer,                  intent(in)    :: ads
    class(mc_lat),            intent(inout) :: lat
    class(energy_parameters), intent(in)    :: e_pars

    ! Hopping
    if (this%hopping%is_defined)      call this%hopping%construct(ads, lat, e_pars, this%beta)
    ! Desorption
    if (this%desorption%is_defined)   call this%desorption%construct(ads, lat, e_pars, this%beta)
    ! Dissociation
!djaif (this%dissociation%is_defined) call this%dissociation%construct(ads, lat, e_pars, this%beta)
    if (this%dissociation%is_defined) call this%dissociation%construct(ads, lat)
    ! Association
!djaif (this%association%is_defined)  call this%association%construct(ads, lat, e_pars, this%beta)
    if (this%association%is_defined)  call this%association%construct(ads, lat)
    ! Bimolecular
!djaif (this%bimolecular%is_defined)  call this%bimolecular%construct(ads, lat, e_pars, this%beta)
    if (this%bimolecular%is_defined)  call this%bimolecular%construct(ads, lat)

  end subroutine construct_1

!-----------------------------------------------------------------------------
  subroutine cleanup_rates(this, ads, lat)
!-----------------------------------------------------------------------------
    class(reaction_type),     intent(inout) :: this
    integer,                  intent(in)    :: ads
    class(mc_lat),            intent(inout) :: lat

    integer :: m

    do m=1,lat%n_max_nn+1 ! '+1' accounts for intra-site hops
      this%hopping%rates(ads,m)%list = 0.0_dp
    end do
    this%desorption%rates(ads) = 0.0_dp
    this%dissociation%rate_info(ads)%list = rate_info_dissociation(0,0,0.0_dp)
    this%dissociation%rate_info(ads)%n_channels = 0
    this%association%rate_info(ads)%list  = rate_info_association(0,0,0.0_dp)
    this%association%rate_info(ads)%n_channels = 0
    this%bimolecular%rate_info(ads)%list  = rate_info_bimolecular(0,0,0.0_dp)
    this%bimolecular%rate_info(ads)%n_channels = 0

  end subroutine cleanup_rates

!-----------------------------------------------------------------------------
  subroutine do_reaction(this, rand, lat, e_pars)
!-----------------------------------------------------------------------------
    class(reaction_type),     intent(inout)     :: this
    real(dp),                 intent(in)        :: rand
    class(mc_lat),            intent(inout)     :: lat
    class(energy_parameters), intent(in)    :: e_pars

    integer :: i, m, ads, iads, reaction_id, channel
!dja    integer :: n_nn, n_nn2, m_nn, col_p2, row_p2, proc, ads_r2
    integer :: n_nn, m_nn, col_p2, row_p2, proc, ads_r2
    integer :: id_r, id_r1, id_r2, id_p1, id_p2
    integer :: ast_p1, ast_p2
    integer :: row, col, row_new, col_new, lst_new, ast_new
    integer :: lst, lst_r2
    logical, dimension(:), allocatable :: rate_update_q

    real(dp) :: u, temp_dp

    ! do nothing if there is nothing to do
    if (this%acc_rate(n_reaction_types) == 0.0_dp) return

    allocate(rate_update_q(lat%n_rows*lat%n_cols))

    ! ------- Select the process

    ! random rate value to select a process
    u = rand*this%acc_rate(n_reaction_types)
    ! determine the type of reaction
    do reaction_id=1,n_reaction_types
      if (u < this%acc_rate(reaction_id)) exit
      print*,reaction_id, this%acc_rate(reaction_id), u
    end do

! Debug printing
!debug(1) = (this%counter(reaction_id,1) == 9 .or. this%counter(reaction_id,1) == 8)
!debug(2) = (this%counter(reaction_id,1) == 8)
!debug(1) = .true.
if (debug(1)) then
  print*
  print*,'do reaction debugging printout:'
  print *, 'Before ', trim(reaction_names(reaction_id)),' #',this%counter(reaction_id,1)+1, ':'
  call lat%print_ocs
  call lat%print_ads
end if
    select case (reaction_id)

      case(hopping_id)
        ! determine hopping channel (adsorbate, direction, ads. site of available ones)
        !                           (      ads,      m_nn,                        iads)
        temp_dp = 0.0_dp
        extloop: do ads=1,this%n_ads_total
          lst = lat%lst(lat%ads_list(ads)%row, lat%ads_list(ads)%col)
          n_nn = lat%n_nn(lst,1)
          id_r  = lat%ads_list(ads)%id

          ! hop to a neighbor cell
          do m_nn=1, n_nn
            ! a new position of particle (ads) after a hop to a neighbor (m_nn)
            call lat%neighbor(ads,m_nn,row_new,col_new)
            lst_new = lat%lst(row_new,col_new)
            do iads=1,size(lat%avail_ads_sites(id_r,lst_new)%list)
              ast_new = lat%avail_ads_sites(id_r,lst_new)%list(iads)
              temp_dp = temp_dp + this%hopping%rates(ads,m_nn)%list(iads)
              if (u < temp_dp) then
                if (debug(1)) then
                  print*
                  print*, 'ads ', ads, 'm_nn ', m_nn, 'iads', iads
                  do i=1, size(this%hopping%rates(ads,:))
                    print*, 'rate ', this%hopping%rates(ads,i)%list(iads)
                  end do
                end if
                exit extloop
              end if
            end do
          end do

          ! hop within the same cell
          do iads=1,size(lat%avail_ads_sites(id_r,lst)%list)
            ast_new = lat%avail_ads_sites(id_r,lst)%list(iads)
            temp_dp = temp_dp + this%hopping%rates(ads,n_nn+1)%list(iads)
            if (u < temp_dp) then
              if (debug(1)) then
                print*
                print*, 'ads ', ads, 'intra-site ', 'iads', iads
                do i=1, size(this%hopping%rates(ads,:))
                  print*, 'rate ', this%hopping%rates(ads,i)%list(iads)
                end do
              end if
              m_nn = n_nn + 1
              exit extloop
            end if
          end do

        end do extloop

        ! Update hopping rate constants
        ! Particle (ads) is going to hop in direction (m_nn) to ads. site (iads)

        ! initialize a logical array indicating adsorbates affected by hop
        rate_update_q = .false.
        ! Put the hopping particle into the list
        rate_update_q(ads) = .true.
        ! scan over neighbors before hop
        call lat%update_neighbors(rate_update_q, ads, lst)

        ! make intra-site hop
        if (m_nn == n_nn + 1) then
          ! Put the adsorbate into a new ast
          lat%ads_list(ads)%ast  = ast_new

        ! make hop to a neighbor
        else
          ! Delete an adsorbate from its old position
          lat%occupations(lat%ads_list(ads)%row, lat%ads_list(ads)%col) = 0
          ! Put the adsorbate in a new position
          lat%ads_list(ads)%row  = row_new
          lat%ads_list(ads)%col  = col_new
          lat%ads_list(ads)%ast  = ast_new
          ! Update adsorbate position
          lat%occupations(row_new, col_new) = ads

          ! scan over neighbors after hop
          call lat%update_neighbors(rate_update_q, ads, lst_new)

        end if

        ! Reset the rate info for ads
        call this%cleanup_rates(ads,lat)

        ! Update rate array for the affected adsorbates
        do i=1,this%n_ads_total
          if (rate_update_q(i)) call this%construct_1(i, lat, e_pars)
        end do

        ! Update reaction counter
        this%counter(reaction_id,id_r) = this%counter(reaction_id,id_r) + 1

      case(desorption_id)
        ! determine desorption channel (adsorbate)
        !                              (      ads)
        temp_dp = this%acc_rate(hopping_id)
        do ads=1,this%n_ads_total
          temp_dp = temp_dp + this%desorption%rates(ads)
          if (u < temp_dp) exit
        end do

        ! Update desorption rate constants
        ! Particle (ads) is going to desorb to nowhere

        lst = lat%lst(lat%ads_list(ads)%row, lat%ads_list(ads)%col)
        ! initialize a logical array indicating adsorbates affected by desorption
        rate_update_q = .false.
        ! scan over neighbors before desorption
        call lat%update_neighbors(rate_update_q, ads, lst)

        ! ---------- Do desorption:
        id_r = lat%ads_list(ads)%id
        ! Delete an adsorbate from the lattice
        lat%occupations(lat%ads_list(ads)%row, lat%ads_list(ads)%col) = 0
        ! Update the number of adsorbates in the lat structure
        lat%n_ads(id_r) = lat%n_ads(id_r) - 1
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
          this%association%rate_info(ads)  = this%association%rate_info(this%n_ads_total)
          this%bimolecular%rate_info(ads)  = this%bimolecular%rate_info(this%n_ads_total)
        end if

        ! account for the tightening the ads. list
        rate_update_q(ads) = rate_update_q(this%n_ads_total)
        ! Update rate array for the affected adsorbates
        do i=1,this%n_ads_total - 1
          if (rate_update_q(i)) call this%construct_1(i, lat, e_pars)
        end do
        ! Reset the rate info for the last adsorbate
        call this%cleanup_rates(this%n_ads_total,lat)
        ! Update the the number of adsorbates in this
        this%n_ads_total = this%n_ads_total - 1

        ! Update reaction counter
        this%counter(reaction_id,id_r) = this%counter(reaction_id,id_r) + 1

      case(dissociation_id)
        ! determine dissociation channel (ads, channel)
        temp_dp = this%acc_rate(desorption_id)
        extloop2: do ads=1,this%n_ads_total
        do channel=1,this%dissociation%rate_info(ads)%n_channels
          temp_dp = temp_dp + this%dissociation%rate_info(ads)%list(channel)%rate
          if (u < temp_dp) exit extloop2
        end do
        end do extloop2

        lst = lat%lst(lat%ads_list(ads)%row, lat%ads_list(ads)%col)
        ! initialize a logical array indicating adsorbates affected by dissociation
        rate_update_q = .false.
        ! Put the dissociating particle into the list
        rate_update_q(ads) = .true.
        ! scan over neighbors before dissociation
        call lat%update_neighbors(rate_update_q, ads, lst)

        ! Do dissociation:
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

        ! scan over neighbors after dissociation
        call lat%update_neighbors(rate_update_q, this%n_ads_total, lat%lst(row_p2,col_p2))

        ! Reset the rate info for ads (product 1)
        call this%cleanup_rates(ads, lat)
        ! Update rate array for the affected adsorbates
        do i=1,this%n_ads_total
          if (rate_update_q(i)) call this%construct_1(i, lat, e_pars)
        end do

        ! Update reaction counter
        this%counter(reaction_id,id_r) = this%counter(reaction_id,id_r) + 1

      case(association_id)
        ! determine association channel (ads, channel)
        temp_dp = this%acc_rate(dissociation_id)
        extloop3: do ads=1,this%n_ads_total
        do channel=1,this%association%rate_info(ads)%n_channels
          temp_dp = temp_dp + this%association%rate_info(ads)%list(channel)%rate
          if (u < temp_dp) exit extloop3
        end do
        end do extloop3

        ! Get the position and number of reactant 2
        m_nn   = this%association%rate_info(ads)%list(channel)%m
        call lat%neighbor(ads,m,row,col)
        ads_r2 = lat%occupations(row,col)
        ! Get reactants' lst
        lst    = lat%lst(lat%ads_list(ads)%row,    lat%ads_list(ads)%col)
        lst_r2 = lat%lst(lat%ads_list(ads_r2)%row, lat%ads_list(ads_r2)%col)

        ! initialize a logical array indicating adsorbates affected by dissociation
        rate_update_q = .false.
        ! account for reactants
        rate_update_q(ads)    = .true.
        rate_update_q(ads_r2) = .true.
        ! scan over reactants' neighbors
        call lat%update_neighbors(rate_update_q,    ads,    lst)
        call lat%update_neighbors(rate_update_q, ads_r2, lst_r2)

        ! --------------Do association:
        ! Find process number
        proc = this%association%rate_info(ads)%list(channel)%proc
        ! Get reactants and products info
        id_r1  = lat%ads_list(ads   )%id
        id_r2  = lat%ads_list(ads_r2)%id
        id_p1  = this%association%channels(proc)%p1
        ast_p1 = this%association%channels(proc)%p1_ast
        ! Replace reactant 1 with product 1 preserving reactant's number
        lat%ads_list(ads)%id  = id_p1
        lat%ads_list(ads)%ast = ast_p1
        lat%n_ads(id_r1) = lat%n_ads(id_r1) - 1
        lat%n_ads(id_p1) = lat%n_ads(id_p1) + 1
        ! Delete reactant 2
        lat%occupations( lat%ads_list(ads_r2)%row, lat%ads_list(ads_r2)%col ) = 0
        lat%n_ads(id_r2) = lat%n_ads(id_r2) - 1
        ! Rearrange adsorbates except when reactant 2 is the last adsorbate
        if (ads_r2 < this%n_ads_total) then
          ! Put the last adsorbate in place of ads_r2
          lat%ads_list(ads_r2) = lat%ads_list(this%n_ads_total)
          ! Update the last adsorbate number in the lattice
          lat%occupations(lat%ads_list(ads_r2)%row, lat%ads_list(ads_r2)%col) = ads_r2
          ! Update the rates arrays
          this%hopping%rates(ads_r2,:)        = this%hopping%rates(this%n_ads_total,:)
          this%desorption%rates(ads_r2)       = this%desorption%rates(this%n_ads_total)
          this%dissociation%rate_info(ads_r2) = this%dissociation%rate_info(this%n_ads_total)
          this%association%rate_info(ads_r2)  = this%association%rate_info(this%n_ads_total)
        end if

        ! Reset the rate info for ads (product 1) and for the last adsorbate
        call this%cleanup_rates(ads,lat)
        call this%cleanup_rates(this%n_ads_total,lat)

        ! account for the tightening the ads. list
        rate_update_q(ads_r2) = rate_update_q(this%n_ads_total)
        ! Update rate array for the affected adsorbates
        do i=1,this%n_ads_total - 1
          if (rate_update_q(i)) call this%construct_1(i, lat, e_pars)
        end do
        ! Update the the number of adsorbates
        this%n_ads_total = this%n_ads_total - 1

        ! Update reaction counter
        this%counter(reaction_id,id_r1) = this%counter(reaction_id,id_r1) + 1
        this%counter(reaction_id,id_r2) = this%counter(reaction_id,id_r2) + 1

      case(bimolecular_id)
        ! determine bimolecular reaction channel (ads, channel)
        temp_dp = this%acc_rate(association_id)
        extloop4: do ads=1,this%n_ads_total
        do channel=1,this%bimolecular%rate_info(ads)%n_channels
          temp_dp = temp_dp + this%bimolecular%rate_info(ads)%list(channel)%rate
          if (u < temp_dp) exit extloop4
        end do
        end do extloop4

        ! Get the position and number of reactant 2
        m_nn   = this%bimolecular%rate_info(ads)%list(channel)%m
        call lat%neighbor(ads,m,row,col)
        ads_r2 = lat%occupations(row,col)
        ! Get reactants' lst
        lst    = lat%lst(lat%ads_list(ads)%row,    lat%ads_list(ads)%col)
        lst_r2 = lat%lst(lat%ads_list(ads_r2)%row, lat%ads_list(ads_r2)%col)

        ! initialize a logical array indicating adsorbates affected by bimolecular reaction
        rate_update_q = .false.
        ! account for reactants
        rate_update_q(ads)    = .true.
        rate_update_q(ads_r2) = .true.
        ! scan over reactants' neighbors
        call lat%update_neighbors(rate_update_q,    ads,    lst)
        call lat%update_neighbors(rate_update_q, ads_r2, lst_r2)


        ! --------------Do bimolecular reaction:
        ! Find process number
        proc = this%bimolecular%rate_info(ads)%list(channel)%proc
        ! Get reactants' and products' info
        id_r1  = lat%ads_list(ads   )%id
        id_r2  = lat%ads_list(ads_r2)%id
        id_p1  = this%bimolecular%channels(proc)%p1
        ast_p1 = this%bimolecular%channels(proc)%p1_ast
        id_p2  = this%bimolecular%channels(proc)%p2
        ast_p2 = this%bimolecular%channels(proc)%p2_ast
        ! Replace reactant 1 with product 1 preserving reactant's number
        lat%ads_list(ads)%id  = id_p1
        lat%ads_list(ads)%ast = ast_p1
        lat%n_ads(id_r1) = lat%n_ads(id_r1) - 1
        lat%n_ads(id_p1) = lat%n_ads(id_p1) + 1
        ! Replace reactant 2 with product 2 preserving reactant's number
        lat%ads_list(ads_r2)%id  = id_p2
        lat%ads_list(ads_r2)%ast = ast_p2
        lat%n_ads(id_r2) = lat%n_ads(id_r2) - 1
        lat%n_ads(id_p2) = lat%n_ads(id_p2) + 1

        ! Reset the rate info for ads (product 1) and for ads_r2 (product 2)
        call this%cleanup_rates(ads,lat)
        call this%cleanup_rates(ads_r2,lat)
        ! Update rate array for the affected adsorbates
        do i=1,this%n_ads_total
          if (rate_update_q(i)) call this%construct_1(i, lat, e_pars)
        end do

        ! Update reaction counter
        this%counter(reaction_id,id_r1) = this%counter(reaction_id,id_r1) + 1
        this%counter(reaction_id,id_r2) = this%counter(reaction_id,id_r2) + 1

      case default
        print*
        print*, "reaction id is ", reaction_id
        print*, "number of reactions is ", n_reaction_types
        stop 'do_reaction (of reaction_class): must never occur!'

    end select

    ! Update the accumulated rates
    call this%accumulate_rates(lat)

    deallocate(rate_update_q)

! Debug printing
if (debug(1)) then
  print*
  print*, 'After ', trim(reaction_names(reaction_id)),' #',this%counter(reaction_id,1), ':'
  call lat%print_ocs
  call lat%print_ads
  print*, 'reactant 1:', ads
!  !if (reaction_id==bimolecular_id) pause
end if


  end subroutine do_reaction

end module reaction_class
