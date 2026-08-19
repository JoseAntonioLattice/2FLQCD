module GF
  use parameters, only : Lt, Lx, Ly, Lz, epsilon, N
  use su3facts
  use observables
  use gauge, Zeta => Z
  implicit none
  integer, parameter, private :: dp = 8
contains

  subroutine wilson_flow_euler(U)
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    type(su3), dimension(4,Lt,Lx,Ly,Lz) :: V
    type(su3) :: expB
    type(su3alg) :: B
    integer :: x(4), mu
    integer :: i, x1, x2, x3, x4
    real(dp) :: S

    open(unit = 666, file = "data/WF.dat")
    print '(A6,2X,A14,2X,A18,2X,A18)', "# step", "t", "action"
    S = energy_density(U)
    write(666, '(I6,2X,F14.6,2X,F18.10,2X,F18.10)') 0, 0.0_dp, S, energy_density_clover(U)

    do i = 1, n
       do x1 = 1, Lt
          do x2 = 1, Lx
             do x3 = 1, Ly
                do x4 = 1, Lz
                   x = [x1,x2,x3,x4]
                   do mu = 1, 4
                      B = Zeta(U,x,mu)
                      B = B*epsilon
                      expB = exp(B)
                      V(mu,x1,x2,x3,x4) = expB * U(mu,x1,x2,x3,x4)
                   end do
                end do
             end do
          end do
       end do
             
       U = V
       S = energy_density(U)
              
       write(666, '(I6,2X,F14.6,2X,F18.10,2X,F18.10)') i, i*epsilon, S, energy_density_clover(U)
       flush(666)
    end do
    close(666)
  end subroutine wilson_flow_euler

  
  subroutine wilson_flow_rk3(U)
    type(su3), intent(inout) :: U(4,Lt,Lx,Ly,Lz)
    type(su3), dimension(4,Lt,Lx,Ly,Lz) :: W1, W2, W3
    type(su3alg), dimension(4,Lt,Lx,Ly,Lz) :: Z0, Z1, Z2, B
    integer :: x,y,t,z, mu, it
    real(dp) :: S
    

    open(unit = 666, file = "data/WF.dat")
    print '(A6,2X,A14,2X,A18,2X,A18)', "# step", "t", "action"
    S = energy_density(U)
    write(666, '(I6,2X,F14.6,2X,F18.10,2X,F18.10,2X)') 0, 0.0_dp, S, energy_density_clover(U)

    wilson_time: do it = 1, N
          
       do t = 1, Lt
          do x = 1, Lx
             do y = 1, Ly
                do z = 1, Lz
                   do mu = 1, 4
                      Z0(mu,t,x,y,z) = epsilon*Zeta(U,[t,x,y,z],mu)
                   end do
                end do
             end do
          end do
       end do
       
       W1 = exp(0.25_dp*Z0)*U
       
       do t = 1, Lt
          do x = 1, Lx
             do y = 1, Ly
                do z = 1, Lz
                   do mu = 1, 4
                      Z1(mu,t,x,y,z) = epsilon*Zeta(W1,[t,x,y,z],mu)
                   end do
                end do
             end do
          end do
       end do
       
       
       B = 8.0_dp/9.0_dp*Z1-17.0_dp/36.0_dp*Z0
       W2 = exp(B)*W1
       
       do t = 1, Lt
          do x = 1, Lx
             do y = 1, Ly
                do z = 1, Lz
                   do mu = 1, 4
                      Z2(mu,t,x,y,z) = epsilon*Zeta(W2,[t,x,y,z],mu)
                   end do
                end do
             end do
          end do
       end do
       
       U = exp(0.75_dp*Z2 - B)*W2
       S = energy_density(U)
       write(666, '(I6,2X,F14.6,2X,F18.10,2X,F18.10)') it, it*epsilon, S, energy_density_clover(U)
       flush(666)
    end do wilson_time
    close(666)
  end subroutine wilson_flow_rk3


end module GF
