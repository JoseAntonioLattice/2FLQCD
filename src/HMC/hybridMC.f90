module hybridMC
  use parameters, only : L, Lt, Lx, Ly, Lz, epsilon, N
  use random
  !use dirac
  use su3facts
  implicit none
  
  integer, parameter :: dp = 8
  real(dp), parameter :: pi = acos(-1.0_dp) 
  complex(dp) :: i = (0.0_dp,1.0_dp)

contains
  
  subroutine leapfrog(U,beta)
    type(su3), intent(in) :: U(4,Lt,Lx,Ly,Lz)
    real(dp), intent(in) :: beta
    type(su3) :: Unew(4,Lt,Lx,Ly,Lz) 
    integer :: t
    type(su3alg), dimension(4,Lt,Lx,Ly,Lz) :: P, Pnew
    real(dp), dimension(4,Lt,Lx,Ly,Lz) :: r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12, &
         r13,r14,r15,r16
    real(dp), dimension(4,Lt,Lx,Ly,Lz) :: g1,g2,g3,g4,g5,g6,g7,g8
    integer :: k 
    Unew = U
    !Create the field chi with Gaussian distribution
    !call random_number(r1)
    !call random_number(r2)
    !call elemental_rgauss(chi,r1,r2)

    !Create the pseudo-fermion phi
    !phi = D(chi)

    !Create the momenta with Gaussian distribution
    call random_number(r1)
    call random_number(r2)
    call random_number(r3)
    call random_number(r4)
    call random_number(r5)
    call random_number(r6)
    call random_number(r7)
    call random_number(r8)
    call random_number(r9)
    call random_number(r10)
    call random_number(r11)
    call random_number(r12)
    call random_number(r13)
    call random_number(r14)
    call random_number(r15)
    call random_number(r16)
    
    call elemental_rgauss(g1,r1,r2)
    call elemental_rgauss(g2,r3,r4)
    call elemental_rgauss(g3,r5,r6)
    call elemental_rgauss(g4,r7,r8)
    call elemental_rgauss(g5,r9,r10)
    call elemental_rgauss(g6,r11,r12)
    call elemental_rgauss(g7,r13,r14)
    call elemental_rgauss(g8,r15,r16)

    call P%init_su3alg(g1,g2,g3,g4,g5,g6,g7,g8)
    
    Pnew = P !- 0.5*epsilon * F(Unew,beta)
    
    do k = 1, N - 1
       Unew = my_exp(  Pnew)
       Pnew = Pnew - epsilon * F(Unew,beta)
    end do
    Unew = my_exp(Pnew)*Unew
    Pnew = Pnew - 0.5*epsilon*F(Unew,beta)
  end subroutine leapfrog

  !subroutine hmc()
  !  call leapfrog()
  !end subroutine hmc


  function F(U,beta)!,chi)
    type(su3alg), dimension(4,Lt,LX,Ly,Lz) :: F
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(in) :: U
    real(dp), intent(in) :: beta
    !complex(dp), dimension(Lt,LX,Ly,Lz), intent(in) :: chi
    type(su3alg) :: Zeta
    integer :: t,x,y,z,mu
    
    do t = 1, Lt
       do x = 1, Lx
          do y = 1, Ly
             do z = 1, Lz
                do mu = 1, 4
                   !F(mu,t,x,y,z) = -0.25_dp*beta*Z(U,[t,x,y,z],mu)
                end do
             end do
          end do
       end do
    end do
    
    
  end function F

  
end module hybridMC
