module dynamics
  use pbc
  use su3facts
  use parameters
  use hybridMC
  use lua
  use observables
  use GF
  implicit none
  integer, parameter, private :: dp = 8
contains

  subroutine initialize()
    use arrays
    integer :: i

    call create_gellmann_matrices()
    call create_gamma_matrices()
    
    call set_pbc(L)
    
    
    allocate(beta(nbeta))
    allocate(U(4,Lt,Lx,Ly,Lz))

    allocate(a_plqv(N_measurements))
    allocate(a_acc_rate(N_measurements))
    if(nbeta == 1) then
       beta(1) = betai
    else
       beta = [(betai + i*(betaf-betai)/(nbeta-1), i = 0, nbeta - 1)]
    end if
  end subroutine initialize
  
  subroutine sweeps(U,beta,acceptance_rate)
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp), intent(in) :: beta
    real(dp), intent(out), optional :: acceptance_rate

    select case(trim(algorithm))
    case("hmc")
       if(present(acceptance_rate)) then
          call hmc(U,beta,acceptance_rate)
       else
          call hmc(U,beta)
       end if
    case("metropolis")
       call sweeps_metropolis(U,beta)
    case("heatbath")
       call sweeps_heatbath(U,beta)
    end select
       
  end subroutine sweeps

  subroutine thermalization(U,beta,acceptance_rate)
    use starts
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp), intent(out), optional :: acceptance_rate
    real(dp) :: beta
    integer :: it
    
    do it = 1, N_thermalization
       if(present(acceptance_rate))then
          call sweeps(U,beta,acceptance_rate)
       else
          call sweeps(U,beta)
       end if
    end do
    print*, "Thermalization done!"

  end subroutine thermalization

  subroutine measurements(U,beta)
    use starts
    use save
    use arrays, only : a_plqv, a_acc_rate
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp) :: acceptance_rate
    real(dp) :: beta
    integer :: im, iskip
    
    do im = 1, N_measurements
       do iskip = 1, n_skip
          call sweeps(U,beta,acceptance_rate)
       end do
       a_plqv(im) = plaquette_value(U)
       a_acc_rate(im) = acceptance_rate
       if(saveconf) call save_configuration(U,beta)
    end do

  end subroutine measurements

  subroutine simulation(U,beta)
    use starts
    use statistics
    use arrays, only : a_plqv, a_acc_rate
    use save
   
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp) :: beta(:)
    integer :: ib, inunit
    character(:), allocatable :: filename
    
    open(newunit=inunit, file= "data/plaquette_value_"//trim(algorithm)//".dat")
    write(inunit,nml = lattice)

    select case(trim(start))
    case("hot")
       call hot_start(U)
    case("cold")
       call cold_start(U)
    case default 
       stop "Not a valid start. Choose 'hot' or 'cold'."
    end select

    if(GFON) then
       !call thermalization(U,beta(1))
       filename = "data/configurations/Lt="//int2str(Lt)// &
            "/Lx="//int2str(Lx)//"/Ly="//int2str(Ly)//"/Lz="//int2str(Lz)// &
            "/beta="//real2str(beta(1),1,4)//"/U_1.bin"
       call read_configuration(U,filename)
       call wilson_flow_rk3(U)
       return
    end if
    
    do ib = 1, size(beta)
       call thermalization(U,beta(ib))
       call measurements(U,beta(ib))
       print*, beta(ib), sum(a_plqv)/real(size(a_plqv)),sum(a_acc_rate)/real(size(a_acc_rate))
       write(inunit,*) beta(ib), avr(a_plqv), jackknife2(a_plqv), avr(a_acc_rate), jackknife2(a_acc_rate)
       flush(inunit)
    end do

  end subroutine simulation
     
end module dynamics
