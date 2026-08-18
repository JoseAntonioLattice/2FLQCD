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

    if(nbeta == 1) then
       beta(1) = betai
    else
       beta = [(betai +i*(betaf-betai)/(nbeta-1), i = 0, nbeta - 1)]
    end if
  end subroutine initialize
  
  subroutine sweeps(U,beta)
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp), intent(in) :: beta
    !character(*), intent(in) :: algorithm

    select case(trim(algorithm))
    case("hmc")
       call hmc(U,beta)
    case("metropolis")
       call sweeps_metropolis(U,beta)
    case("heatbath")
       call sweeps_heatbath(U,beta)
    end select
       
  end subroutine sweeps

  subroutine thermalization(U,beta)
    use starts
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp) :: beta
    integer :: it
    
    do it = 1, N_thermalization
       call sweeps(U,beta)
    end do
    print*, "Thermalization done!"

  end subroutine thermalization

  subroutine measurements(U,beta)
    use starts
    use arrays, only : a_plqv
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp) :: beta
    integer :: im, iskip
    
    do im = 1, N_measurements
       do iskip = 1, n_skip
          call sweeps(U,beta)
       end do
       a_plqv(im) = plaquette_value(U)
       !print*, a_plqv(im)
    end do

  end subroutine measurements


  subroutine simulation(U,beta)
    use starts
    use statistics
    use arrays, only : a_plqv
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp) :: beta(:)
    integer :: ib, inunit

    open(newunit=inunit, file= "data/plaquette_value_"//trim(algorithm)//".dat")

    select case(trim(start))
    case("hot")
       call hot_start(U)
    case("cold")
       call cold_start(U)
    case default 
       stop "Not a valid algorithm"
    end select

    if(GFON) then
       call thermalization(U,beta(1))
       call wilson_flow_rk3(U)
       return
    end if
    
    do ib = 1, size(beta)
       call thermalization(U,beta(ib))
       call measurements(U,beta(ib))
       print*, beta(ib), sum(a_plqv)/real(size(a_plqv))
       write(inunit,*) beta(ib), avr(a_plqv), jackknife2(a_plqv)
       flush(inunit)
    end do

  end subroutine simulation
   
  
end module dynamics
