module hybridMC
  use parameters, only : L, Lt, Lx, Ly, Lz, epsilon, N
  use random
  use dirac
  implicit none
  
  integer, parameter :: dp = 8
  real(dp), parameter :: pi = acos(-1.0_dp) 
  complex(dp) :: i = (0.0_dp,1.0_dp)

contains
  
  subroutine leapfrog(U)
    complex(dp), intent(in) :: U(4,Lt,Lx,Ly,Lz)
    complex(dp) :: Unew(4,Lt,Lx,Ly,Lz) 
    integer :: t
    real(dp), dimension(4,Lx,Ly,Lz) :: P, Pnew, r1, r2

    Unew = U
    !Create the field chi with Gaussian distribution
    call random_number(r1)
    call random_number(r2)
    call elemental_rgauss(chi,r1,r2)

    !Create the pseudo-fermion phi
    phi = D(chi)

    !Create the momenta with Gaussian distribution
    call random_number(r1)
    call random_number(r2)
    call elemental_rgauss(P,r1,r2)
    
    Pnew = P - 0.5*epsilon * F(Unew,phi)
    
    do k = 1, N - 1
       Unew = exp( i*epsilon * Pnew)
       Pnew = Pnew - epsilon * F(Unew,phi)
    end do
    Unew = exp(i*epsilon*Pnew)*Unew
    Pnew = Pnew - 0.5*epsilon*F(Unew,phi)
  end subroutine leapfrog

  subroutine hmc()

    call leapfrog()
  end subroutine hmc


  

  
end module hybridMC
