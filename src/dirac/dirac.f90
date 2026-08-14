module dirac
  use parameters, only: m0, Lx, Ly,Lz, Lt
  use pbc
  use su3facts, gamma => dirac_matrix
  
  implicit none
  integer, parameter, private :: dp = 8
contains



  function D(chi,U)
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(in) :: U
    complex(dp), dimension(4,3,Lt,Lx,Ly,Lz), intent(in) :: chi
    complex(dp), dimension(4,3,Lt,Lx,Ly,Lz) :: D
    integer :: t, x, y, z, mu, alpha, beta, xp(4), xm(4)

    complex(dp) :: res(3), dag(3,3)
    
    do t = 1, Lt
       do x = 1, Lx
          do y = 1, Ly
             do z = 1, Lz
                do alpha = 1, 4
                   res = (0.0_dp,0.0_dp)
                   do mu = 1, 4
                      xp = ip([t,x,y,z],mu)
                      xm = im([t,x,y,z],mu)
                      dag = U(mu,xm(1),xm(2),xm(3),xm(4))%mat
                      dag = conjg(transpose(dag))
                      do beta = 1, 4
                         
                         res = res  - 0.5_dp*( &
                              (delta_4x4(alpha,beta) - gamma(mu)%mat(alpha,beta)) * &
                              matmul(U(mu,t,x,y,z)%mat,chi(beta,:,xp(1),xp(2),xp(3),xp(4)) * sgnp(mu,t) ) + &
                              (delta_4x4(alpha,beta) + gamma(mu)%mat(alpha,beta)) * &
                              matmul(dag,chi(beta,:,xm(1),xm(2),xm(3),xm(4)) * sgnm(mu,t)) )
                      end do
                   end do
                   D(alpha,:,t,x,y,z) = res + (m0 + 4.0_dp)*chi(alpha,:,t,x,y,z)
                end do
             end do
          end do
       end do
    end do
                      
  end function D

  function Ddagger(chi,U)
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(in) :: U
    complex(dp), dimension(4,3,Lt,Lx,Ly,Lz), intent(in) :: chi
    integer :: t, x, y, z, mu, alpha, beta, xp(4), xm(4)
    complex(dp) :: Ddagger(4,3,Lt,Lx,Ly,Lz), res(3), dag(3,3)
    
    do t = 1, Lt
       do x = 1, Lx
          do y = 1, Ly
             do z = 1, Lz
                do alpha = 1, 4
                   res = (0.0_dp,0.0_dp)
                   do mu = 1, 4
                      xp = ip([t,x,y,z],mu)
                      xm = im([t,x,y,z],mu)
                      dag = U(mu,xm(1),xm(2),xm(3),xm(4))%mat
                      dag = conjg(transpose(dag))
                      do beta = 1, 4
                         res = res  - 0.5_dp*( &
                              (delta_4x4(alpha,beta) - gamma(mu)%mat(alpha,beta)) * &
                              matmul(dag, chi(beta,:,xp(1),xp(2),xp(3),xp(4)) * sgnp(mu,t)) + &
                              (delta_4x4(alpha,beta) + gamma(mu)%mat(alpha,beta)) * &
                              matmul(U(mu,t,x,y,z)%mat, chi(beta,:,xm(1),xm(2),xm(3),xm(4)) * sgnm(mu,t)) )
                      end do
                   end do
                   Ddagger(alpha,:,t,x,y,z) = res + (m0 + 4.0_dp)*chi(alpha,:,t,x,y,z)
                end do
             end do
          end do
       end do
    end do
                      
  end function Ddagger
  
end module dirac
