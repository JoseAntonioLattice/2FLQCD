module hybridMC
  use parameters, only : L, Lt, Lx, Ly, Lz, epsilon, N
  use random
  use dirac
  use su3facts
  use CG
  use gauge,  Zeta => Z
  implicit none
  
  integer, parameter, private :: dp = 8
  real(dp), parameter, private :: pi = acos(-1.0_dp) 
  complex(dp), parameter, private :: i = (0.0_dp,1.0_dp)

contains
  
  subroutine leapfrog(U,beta)
    type(su3), intent(inout) :: U(4,Lt,Lx,Ly,Lz)
    real(dp), intent(in) :: beta
    type(su3) :: Unew(4,Lt,Lx,Ly,Lz) 
    integer :: t
    type(su3alg), dimension(4,Lt,Lx,Ly,Lz) :: P, Pnew
    real(dp), dimension(4,Lt,Lx,Ly,Lz) :: r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12, &
         r13,r14,r15,r16
    real(dp), dimension(4,Lt,Lx,Ly,Lz) :: g1,g2,g3,g4,g5,g6,g7,g8
    complex(dp), dimension(4,3,Lt,Lx,Ly,Lz) :: chi , phi, psi ,psiold
    real(dp), dimension(4,3,Lt,Lx,Ly,Lz) :: s1,s2,s3,s4
    integer :: k
    real(dp) :: DS, r

    !Create the field chi with Gaussian distribution
    call random_number(s1)
    call random_number(s2)
    call random_number(s3)
    call random_number(s4)
    chi = (0.0_dp,0.0_dp)
    call elemental_rgauss(chi%re,s1,s2)
    call elemental_rgauss(chi%im,s3,s4)
    
    !Create the pseudo-fermion phi
    phi = D(chi,U)

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
    
    Unew  = U
    psi = conjugate_gradient(phi,Unew)
    chi = Ddagger(psi,Unew)
    
    Pnew = P - 0.5*epsilon * F(U,psi,chi,beta)
    psiold = psi
    do k = 1, N - 1
       Unew = exp( i*epsilon*Pnew)*Unew
       psi = conjugate_gradient(phi,Unew)
       chi = Ddagger(psi,Unew)
       Pnew = Pnew - epsilon * F(Unew,psi,chi,beta)
    
    end do
    Unew = exp(i*epsilon*Pnew)*Unew
       
    psi = conjugate_gradient(phi,Unew)
    chi = Ddagger(psi,Unew)
    Pnew = Pnew - 0.5*epsilon*F(Unew,psi,chi,beta)

    DS = sum(tr(P*P - Pnew*Pnew)) - DeltaS(U,Unew,beta) + sum(conjg(phi)*psiold) &
         -sum(conjg(psi)*phi)

    call random_number(r)
    if( r <= min(1.0_dp,exp(DS))) U = Unew
    
  end subroutine leapfrog

  subroutine hmc(U,beta)
    type(su3), intent(inout) :: U(4,Lt,Lx,Ly,Lz)
    real(dp), intent(in) :: beta
    call leapfrog(U,beta)
  end subroutine hmc


  function F(U,psi,chi,beta)
    type(su3alg), dimension(4,Lt,Lx,Ly,Lz) :: F
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(in) :: U
    !complex(dp), dimension(4,3,Lt,Lx,Ly,Lz), intent(in) :: phi
    complex(dp), dimension(4,3,Lt,Lx,Ly,Lz), intent(in) :: chi, psi
    real(dp), intent(in) :: beta
    integer :: t,x,y,z,mu,a,b, alpha,bet, xp(4)
    complex(dp), dimension(3,3) :: res1, res2
    type(su3alg) :: W, ZetaU
    type(matrix3x3) :: WTA
    type(su3) :: dagU
    


    do t = 1, Lt
       do x = 1, Lx
          do y = 1, Ly
             do z = 1, Lz
                do mu = 1, 4
                   xp = ip([t,x,y,z],mu)
                   dagU = dagger(U(mu,t,x,y,z))
                   res1 = (0.0_dp,0.0_dp)
                   res2 = (0.0_dp,0.0_dp)
                   do alpha = 1, 4
                      do bet = 1, 4
                         do a = 1, 3
                            do b = 1, 3
                               res1(a,b) = res1(a,b) + (delta_4x4(alpha,bet)-gamma(mu)%mat(alpha,bet))* &
                                    chi(bet,a,xp(1),xp(2),xp(3),xp(4))*conjg(psi(alpha,b,t,x,y,z))
                               res2(a,b) = res2(a,b) + (delta_4x4(alpha,bet)+gamma(mu)%mat(alpha,bet))* &
                                    chi(bet,a,t,x,y,z)*conjg(psi(alpha,b,xp(1),xp(2),xp(3),xp(4)))
                            end do
                         end do
                      end do
                   end do

                   W%mat = sgnp(mu,t)*matmul(U(mu,t,x,y,z)%mat,res1) - sgnm(mu,xp(1))*matmul(res2,dagU%mat)
                   WTA = i*TA(W)
                   ZetaU = beta*i/6.0_dp*Zeta(U,[t,x,y,z],mu)
                   f(mu,t,x,y,z)%mat = ZetaU%mat + WTA%mat  
                end do
             end do
          end do
       end do
    end do
    
    
  end function F

  function DeltaS(U,Unew,beta)
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(in) :: U, Unew 
    real(dp), intent(in) :: beta
    integer :: t,x,y,z,mu, nu
    real(dp) :: DeltaS
    
    DeltaS = 0.0_dp
    do t = 1, Lt
       do x = 1, Lx
          do y = 1, Ly
             do z = 1, Lz
                do mu = 1, 3
                   do nu = mu+1, 4
                      DeltaS = DeltaS + real(tr(plaquette(U,[t,x,y,z],mu,nu) - plaquette(Unew,[t,x,y,z],mu,nu)))
                   end do
                end do
             end do
          end do
       end do
    end do

    DeltaS = beta/3.0_dp * DeltaS
    
    
  end function DeltaS

  
end module hybridMC
