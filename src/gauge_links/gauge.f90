module gauge
  use pbc, only : ip, im
  implicit none
  integer, parameter, private :: dp = 8
contains

  function dagger(U)
    
  end function dagger
  
  function plaquette(U,x,mu,nu)
    complex(dp), dimension(4,Lt,Lx,Ly,Lz), intent(in) :: U
    integer, intent(in) :: x(4), mu, nu
    integer, dimension(4) :: x2, x3

    x2 = ip(x,mu)
    x3 = ip(x,nu)
    
    plaquette = U(mu,x(1),x(2),x(3),x(4)) * &
                U(nu,x2(1),x2(2),x2(3),x2(4)) * &
         dagger(U(mu,x3(1),x3(2),x3(3),x3(4)) * &
         dagger(U(mu,x(1),x(2),x(3),x(4))
  end function plaquette

end module gauge
