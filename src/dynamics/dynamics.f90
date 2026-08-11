module dynamics
  use pbc
  use su3facts
  use parameters
  use hybridMC
  implicit none

contains

  subroutine initialize()
    use arrays
    
    call set_pbc(L)

    allocate(U(4,Lt,Lx,Ly,Lz)) 
    
  end subroutine initialize
  
  subroutine sweeps(U,beta)
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp) :: beta

    call hmc(U,beta)

  end subroutine sweeps

  subroutine thermalization(U,beta)
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp) :: beta
    integer :: it
    
    do it = 1, N_thermalization
       call sweeps(U,beta)
    end do

  end subroutine thermalization

  subroutine measurements(U,beta)
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp) :: beta
    integer :: im, iskip
    
    do im = 1, N_measurements
       do iskip = 1, n_skip
          call sweeps(U,beta)
       end do
       
    end do

  end subroutine measurements


  subroutine simulation(U,beta)
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp) :: beta(:)
    integer :: ib

    do ib = 1, size(beta)
       call measurements(U,beta(ib))
    end do


  end subroutine simulation
    
  
end module dynamics
