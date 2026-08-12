module dynamics
  use pbc
  use su3facts
  use parameters
  use hybridMC
  use lua
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
    
    beta = [(betai +i*(betaf-betai)/(nbeta-1), i = 0, nbeta - 1)]
    
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
    use arrays, only : a_plqv
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp) :: beta(:)
    integer :: ib, inunit

    open(newunit=inunit, file= "data/plaquette_value.dat")
    
    call hot_start(U)
    do ib = 1, size(beta)
       call thermalization(U,beta(ib))
       call measurements(U,beta(ib))
       print*, beta(ib), sum(a_plqv)/real(size(a_plqv))
       write(inunit,*) beta(ib), sum(a_plqv)/real(size(a_plqv))
       flush(inunit)
    end do

  end subroutine simulation
   
  function plaquette_value(U) result(plqv)
    use gauge, only : plaquette
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp) :: plqv
    integer :: t,x,y,z, mu, nu

    plqv = 0.0_dp
    do t = 1, Lt
       do x = 1, Lx
          do y = 1, Ly
             do z = 1, Lz
                do mu = 1, 3
                   do nu = mu + 1, 4
                      plqv = plqv + real(tr(plaquette(U,[t,x,y,z],mu,nu)))
                   end do
                end do
             end do
          end do
       end do
    end do
    
    plqv = plqv/(3*6.0_dp*product(L))
             
  end function plaquette_value
  
end module dynamics
