module CG
  use parameters, only : Lt, Lx, Ly, Lz, m0
  use dirac, Diracmat => D
  use su3facts
  implicit none
  integer, parameter, private :: dp = 8
  integer, parameter, private :: max_iter = 1000
  real(dp), parameter, private :: tol = 1.0e-6_dp
contains

  function conjugate_gradient(phi,U) result(x) 

    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(in) :: U
    complex(dp), dimension(4,3,Lt,Lx,Ly,Lz), intent(in) :: phi
    complex(dp), dimension(4,3,Lt,Lx,Ly,Lz) :: x
    complex(dp), dimension(4,3,Lt,Lx,Ly,Lz) :: r,Ap, p
    real(dp) :: alpha, pAp
    real(dp)    :: phi_norm2, beta, rr, rr_new
    integer :: k

    
    !Ap = Ddagger(Diracmat(phi,U),U)
    !x = phi
        
    !r = x - Ap
    x = (0.0_dp,0.0_dp)
    r = phi
    p = r

    rr = real(sum(r*conjg(r)))
    phi_norm2 = real(sum(phi*conjg(phi)))


    do k = 1, max_iter
       
       Ap = Diracmat(Ddagger(p,U),U)
       pAp = sum(p*conjg(Ap))
       
       alpha = rr/pAp
       x = x + alpha*p
       r = r - alpha*Ap
       
       rr_new = real(sum(r*conjg(r)))
       
       if( rr_new < tol**2*phi_norm2 ) return
       
       beta = rr_new/rr
       p = r + beta*p
       
       rr = rr_new
    end do
    print*, "Max iterations", max_iter," reached, No convergence. Residual: ", sqrt(rr_new/phi_norm2)
    
  end function conjugate_gradient

  
end module CG
