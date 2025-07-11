module rates_desorption_class

  use constants
  use control_parameters_class
  use mc_lat_class
  use energy_parameters_class
  use energy_mod
  use open_file
  use utilities
  use rate_constant_laws

  implicit none

  private
  public    :: desorption_init, desorption_type

  type :: desorption_type

    logical :: is_defined = .false.
    ! Desorption Rates
    !                          n_adsorbate              -> which particle
    !                          .
    real(dp), dimension(:), allocatable :: rates

    !                         (   n_species                       -> which species
    !                             .  n_site_type                  -> where from
    !                             .  .  n_adsorption_sites )
    !                             .  .  .
    real(dp),           dimension(:, :, :), allocatable :: process
    type(int_law_pars), dimension(:, :, :), allocatable :: rate_corr_pars

  contains
    procedure :: construct
    procedure :: print

  end type


contains
!------------------------------------------------------------------------------
  function desorption_init(c_pars, lat, e_pars)
!------------------------------------------------------------------------------
    type(desorption_type) desorption_init

    type(control_parameters), intent(in)    :: c_pars
    type(mc_lat)            , intent(in)    :: lat
    type(energy_parameters) , intent(in)    :: e_pars

    integer :: ios, nwords, line_number, i, j, i1, i2

    integer :: species, st1, ast1
    logical :: e_defined1, r_defined, undefined_rate, undefined_energy

    character(len=max_string_length)                :: buffer
    character(len=max_string_length)                :: words(100)
    character(len=len(trim(c_pars%rate_file_name))) :: file_name

    character(len=10)     :: current_species_name
    integer               :: current_species_id
    integer               :: current_law_id
    integer               :: rct_law_id, rct_law_id_glob
    integer               :: rcic_law_id, rcic_law_id_glob

    integer               :: parse_state
    integer, parameter    :: parse_state_ignore     = -1
    integer, parameter    :: parse_state_default    =  0
    integer, parameter    :: parse_state_desorption =  desorption_id

    logical :: rct_law_defined  = .false.
    logical :: rcic_law_defined = .false.

    real(dp), dimension(n_max_rct_pars )::  rct_pars,  rct_pars_glob
    real(dp), dimension(n_max_rcic_pars):: rcic_pars, rcic_pars_glob

    integer,  parameter   :: default_int = 0
    real(dp), parameter   :: default_rate  = -1.0_dp


    ! Allocate rates array
    allocate( desorption_init%rates(lat%n_rows*lat%n_cols) )
    allocate( desorption_init%process( c_pars%n_species,&
                                 n_max_lat_site_types, n_max_ads_sites) )

    desorption_init%process = default_rate
    desorption_init%rates    = 0.0_dp

    allocate(desorption_init%rate_corr_pars( c_pars%n_species,&
                                 n_max_lat_site_types, n_max_ads_sites ))

    desorption_init%rate_corr_pars = int_law_pars(default_int, default_rate)


    !  read rate definitions from the input file
    file_name = c_pars%rate_file_name
    call open_for_read(inp_unit, file_name )

    ios = 0
    parse_state = parse_state_default
    line_number = 0
    undefined_energy = .false.

    !---------------------------------------------------------------------------
    ! loop over all lines of rate input file
    !---------------------------------------------------------------------------

    do while (ios == 0)

      read(inp_unit, '(A)', iostat=ios) buffer

      line_number = line_number + 1
      ! ios = 0: valid record read
      ! ios < 0: end of record condition encountered or end of file condition detected
      ! ios > 0: an error is detected
      !if (ios < 0) buffer = section_end  ! treat end of file as the section end

        ! Split an input string
        words = ''
        call split_string(buffer, words, nwords)

      ! skip comments
      if (nwords == 0) cycle

      select case (parse_state)

        case(parse_state_default)
          ! in parse state default:
          !    word 'desorption' to mark beginning of a desorption section
          !    ignore anything else until desorption section begins

          if (words(1) == reaction_names(desorption_id)) then

            desorption_init%is_defined = .true.
            parse_state = parse_state_desorption
            ! reset necessary to allow multiple desorption sections
            rct_law_defined  = .false.
            rcic_law_defined = .false.
            rct_law_id_glob  = 0
            rcic_law_id_glob = 0

            if (nwords == 2) then
              read(words(2),'(A)') current_species_name
              current_species_id = get_index(current_species_name, c_pars%ads_names)
              if (current_species_id == 0) call error_message(file_name, line_number, buffer, &
                                                              "unknown species in desorption section definition")
            else
              call error_message(file_name, line_number, buffer, &
                         "desorption key must have 1 parameter -- species")
            end if !nwords == 2

          end if  ! words(1)

        case(parse_state_desorption)
          ! process:
          !    temperature law records,
          !    interaction law records,
          !    'from' records,
          !    section end


          if (words(1) == section_end) then
            parse_state = parse_state_default

          elseif (words(1) =='temperature_law') then
            rct_law_id_glob = get_index(words(2), rct_law_names)
            if (rct_law_id_glob == 0) then
              call error_message(file_name, line_number, buffer,&
                                 "invalid temperature law statement")
            else
              rct_law_defined = .true.
              select case (rct_law_id_glob)
                case (Arrhenius_id)
                  if (nwords/=4) call error_message(file_name, line_number, buffer,&
                                                    "Arrhenius must have 2 parameters")
                  read(words(3),*) rct_pars_glob(1)
                  read(words(4),*) rct_pars_glob(2)
                case (extArrhenius_id)
                  if (nwords/=5) call error_message(file_name, line_number, buffer,&
                                                    "extArrhenius must have 3 parameters")
                  read(words(3),*) rct_pars_glob(1)
                  read(words(4),*) rct_pars_glob(2)
                  read(words(5),*) rct_pars_glob(3)
              end select
            endif

          elseif (words(1) =='interaction_law') then
            rcic_law_id_glob = get_index(words(2), rcic_law_names)
            if (rcic_law_id_glob == 0) then
              call error_message(file_name, line_number, buffer,&
                                 "invalid interaction law statement")
            else
              rcic_law_defined = .true.
              select case (rcic_law_id_glob)
                case (rcic_linear_id)
                  if (nwords/=4) call error_message(file_name, line_number, buffer,&
                                                    "linear interaction law must have 2 parameters")
                  read(words(3),*) rcic_pars_glob(1)
                  read(words(4),*) rcic_pars_glob(2)
              end select
            endif

          else
            ! ---------------------------------------------------------
            ! check if we have a valid from rate record
            ! ---------------------------------------------------------
            if (nwords < 2) &
              call error_message(file_name, line_number, buffer, &
                                 "invalid 'from' rate record in the desorption section")
            i1 = get_index(words(1),lat_site_names)
            i2 = get_index(words(2),ads_site_names)

            ! check for invalid site name (invalid lst or ast in either from or to site)
            if ( i1==0 .or. i2==0) then
              print *, 'parse state: ', parse_state
              call error_message(file_name, line_number, buffer, &
                                 "wrong site name in the desorption section")
            end if

            ! ---------------------------------------------------------
            ! check for duplicate entry
            ! ---------------------------------------------------------
            if (desorption_init%process(current_species_id,i1,i2 ) /= default_rate)&
            call error_message(file_name, line_number, buffer, "duplicated entry")

            ! check energy is defined for initial lst and ast
            if( e_pars%ads_energy(current_species_id, i1, i2) == e_pars%undefined_energy ) then
              call error_message(file_name, line_number, buffer, &
                               "rate defined for a site with undefined adsorption energy", &
                               stop=.false., warning=.false.)

              undefined_energy = .true.
            end if


            ! ---------------------------------------------------------
            ! we have a valid rate record. Process it
            ! ---------------------------------------------------------
            ! if record only has 'from' site information then set the law ids and pars to global default values
            if (nwords == 2) then
              if (rct_law_defined .and. rcic_law_defined) then
                rct_law_id  = rct_law_id_glob
                rct_pars    = rct_pars_glob
                rcic_law_id = rcic_law_id_glob
                rcic_pars   = rcic_pars_glob
              elseif ( .not. rct_law_defined .and. .not. rcic_law_defined ) then
                  call error_message(file_name, line_number, buffer, &
                                 "temperature and interaction laws are not defined")
              elseif (.not. rct_law_defined) then
                  call error_message(file_name, line_number, buffer, &
                                 "temperature law is not defined")
              else
                  call error_message(file_name, line_number, buffer, &
                                 "interaction law is not defined")
              end if

            else
              ! set rct_law_id and rcic_law_id to 0 to indicate no law on this line found yet
              rct_law_id  = 0
              rcic_law_id = 0

              do i=3,nwords
                ! check  for rct law on this line
                if (get_index(words(i), rct_law_names) /= 0) then
                  rct_law_id = get_index(words(i), rct_law_names)
                  select case (rct_law_id)
                    case (Arrhenius_id)
                      do j=1,2
                        if ( .not. read_num(words(i+j),rct_pars(j)) )&
                          call error_message(file_name, line_number, buffer,&
                                                  "Arrhenius must have 2 numerical parameters")
                      end do
                    case (extArrhenius_id)
                      do j=1,3
                        if ( .not. read_num(words(i+j),rct_pars(j)) )&
                          call error_message(file_name, line_number, buffer,&
                                                  "extArrhenius must have 3 numerical parameters")
                      end do
                    case default
                      call error_message(file_name, line_number, buffer, "This should not happen! Check the code!")
                  end select

                end if

                ! check if we have an rcic law on this line
                if (get_index(words(i), rcic_law_names) /= 0) then
                  rcic_law_id = get_index(words(i), rcic_law_names)
                  select case (rcic_law_id)
                    case (rcic_linear_id)
                      do j=1,2
                        if ( .not. read_num(words(i+j),rcic_pars(j)) )&
                          call error_message(file_name, line_number, buffer,&
                                                  "linear interaction must have 2 numerical parameters")
                      end do
                    case default
                      call error_message(file_name, line_number, buffer, "This should not happen! Check the code!")
                  end select
                end if

              end do ! i=3,nwords

              ! if rct_law  or rcic_law are not defined on this line, set to global defaults
              if (rct_law_id == 0) then
                rct_law_id = rct_law_id_glob
                rct_pars   = rct_pars_glob
              end if
              if (rcic_law_id == 0) then
                rcic_law_id = rcic_law_id_glob
                rcic_pars   = rcic_pars_glob
              end if

              ! check for valid temperature and interaction laws
              if (rct_law_id == 0) &
                call error_message(file_name, line_number, buffer, &
                                  "no temperature law is specified")
              if (rcic_law_id == 0) &
                call error_message(file_name, line_number, buffer, &
                                   "no interaction law is specified")

            end if ! (nwords==2)

            ! Check if rct and rcic laws are properly set
            if (rct_law_id == 0) &
              call error_message(file_name, line_number, buffer, "invalid temperature law")
            if (rcic_law_id == 0) &
              call error_message(file_name, line_number, buffer, "invalid interaction law")

          ! Set rate constants and rcic for desorption channels
            desorption_init%process(current_species_id,i1,i2 ) = &
                        rct_law(rct_law_id, c_pars%temperature, rct_pars)

            desorption_init%rate_corr_pars(current_species_id,i1,i2)%id   = rcic_law_id
            desorption_init%rate_corr_pars(current_species_id,i1,i2)%pars = rcic_pars

          endif ! words(1) == section_end

      end select

    end do ! while ios=0

    close(inp_unit)

    if (parse_state /= parse_state_default) then
      print *, 'parse state: ', parse_state
      call error_message(file_name, line_number, buffer, &
           'desorption: error, incomplete desorption section')
    endif

    if (.not. undefined_energy) then

      if (desorption_init%is_defined) then
        write(6, '(A/)') ' desorption: passed check that energies are defined for all rates'
      else
        write(6, '(A)') ' no desorption reactions'
      end if

    end if

    ! ---------------------------------------------------------------------------------------------
    ! Check the input consistency
    ! ---------------------------------------------------------------------------------------------
    !
    ! Check if rates are defined for all values of site_types, ads_sites
    !   for which adsorption energies are defined
    ! Adsorpton energies (n_species x n_site_type x n_adsorption_sites)
    ! real(dp), dimension(:,:,:), allocatable :: ads_energy
    !
    ! Note ads_energies could be allocated of n_site_types rather than n_max_lat_site_types
    !

    if (desorption_init%is_defined) then

      undefined_rate = .false.
      do species   = 1, c_pars%n_species
      do st1       = 1, n_max_lat_site_types
      do ast1      = 1, n_max_ads_sites

        e_defined1 = e_pars%ads_energy(species, st1, ast1) /= e_pars%undefined_energy
        r_defined  = desorption_init%process (species, st1, ast1) /= default_rate
        ! Replace default rates with zeros to escape negative rates
        if (.not. r_defined) desorption_init%process (species, st1, ast1) = 0.0_dp

        if ( e_defined1 .and. (.not. r_defined)) then
          if (.not. undefined_rate) then
            undefined_rate = .true.
            print '(2A)',  ' Desorption: warning missing rate definitions in the file ', file_name
            print '(A)',  ' Missing definitions:'
            print '(6x, A)', 'ads  lat_site    ads_site'
          end if

          print '(6x, a5, A10, 2x, a3, 6x, a10, 2x, a3, 6x, L1)' ,            &
                  c_pars%ads_names(species),             &
                  lat_site_names(st1), ads_site_names(ast1)
        end if

      end do
      end do
      end do

      if(undefined_rate) then
        print '(A)', ' Add rate definitions if they should be present'
        print *
      else
        print '(A)', ' Desorption: passed rates consistency check'
        print *
      end if

      if(undefined_energy) then
        print '(A)', ' Desorption: error - undefined energy -- aborting execution'
        stop 997
      end if

    end if

  end function desorption_init

!-----------------------------------------------------------------------------
  subroutine construct(this, ads, lat, e_pars, beta)
!-----------------------------------------------------------------------------
    class(desorption_type), intent(inout) :: this
    integer, intent(in) :: ads
    class(mc_lat), intent(inout) :: lat
    class(energy_parameters), intent(in) :: e_pars
    real(dp), intent(in) :: beta

    integer :: id
    integer :: row, col, lst, ast
    real(dp):: int_energy, int_energy_ts, delta_eps

    row = lat%ads_list(ads)%row
    col = lat%ads_list(ads)%col
    lst = lat%lst(row,col)
    ast = lat%ads_list(ads)%ast
    id  = lat%ads_list(ads)%id

    ! Calculate interaction correction
    int_energy    = energy(ads,lat,e_pars) - e_pars%ads_energy(id, lst, ast)
    int_energy_ts = rcic_law(this%rate_corr_pars(id, lst, ast), int_energy, 0.0_dp)

    ! Barrier correction due to the perturbation
    delta_eps = int_energy_ts - int_energy


    this%rates(ads) = this%process(id, lst, ast)*exp( -beta*delta_eps)
!   print*
!   print*, 'id ',id, ' lst ',lst,' ads. site ', ast
!   print*, 'rate ',this%process(id, lst, ast)

!   print '(A,i4)', 'ads ',ads
!   print '(A,e18.4)',' Energy = ', energy(ads,lat,e_pars)
!   print *,' rate = ', this%rates(ads)
!   write(*,*) 'pause'
!   read(*,*)


  end subroutine construct

!------------------------------------------------------------------------------
  subroutine print(this, c_pars)
!------------------------------------------------------------------------------
    class(desorption_type), intent(in) :: this

    class(control_parameters), intent(in) :: c_pars

    integer :: i, i1, i2

    print*, 'Desorption Rates:'
    do i=1,c_pars%n_species
      print '(/A)','---------------------------'
      print '( A,A)','species: ', c_pars%ads_names(i)
      print '(A)', '---------------------------'
      do i1=1,n_max_lat_site_types
      do i2=1,n_max_ads_sites
        if (this%process(i,i1,i2)< 0.0_dp) then
          cycle
        else
          write(*,'(A,A,1pe11.2)') &
              lat_site_names(i1), ads_site_names(i2), this%process(i,i1,i2)
        end if
      end do
      end do
    end do
    print*

  end subroutine print

end module rates_desorption_class
