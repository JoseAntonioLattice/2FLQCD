module gauge
  use su3facts 
  use pbc, only : ip, im
  implicit none
  integer, parameter, private :: dp = 8
contains

  
  function plaquette(U,x,mu,nu)
    type(su3) :: plaquette
    type(su3), dimension(:,:,:,:,:), intent(in) :: U
    integer, intent(in) :: x(4), mu, nu
    integer, dimension(4) :: x2, x3

    x2 = ip(x,mu)
    x3 = ip(x,nu)
    
    plaquette = U(mu,x(1),x(2),x(3),x(4)) * &
                U(nu,x2(1),x2(2),x2(3),x2(4)) * &
         dagger(U(mu,x3(1),x3(2),x3(3),x3(4))) * &
         dagger(U(nu,x(1),x(2),x(3),x(4)))
  end function plaquette

  function staples(U,x,mu)
    type(matrix3x3) :: staples
    type(su3), dimension(:,:,:,:,:), intent(in) :: U
    integer, intent(in) :: x(4), mu
    integer, dimension(4) :: x2, x3,x4,x5
    integer :: nu
    type(su3) :: sum
    
    call sum%zero
    x3 = ip(x,mu)
    do nu = 1, 4
       if( nu == mu ) cycle
       x2 = ip(x, nu)
       x4 = im(x, nu)
       x5 = im(x3,nu)
       sum = sum + &
            U(nu,x(1),x(2),x(3),x(4))* &
            U(mu,x2(1),x2(2),x2(3),x2(4))* &
     dagger(U(nu,x3(1),x3(2),x3(3),x3(4))) + &
     dagger(U(nu,x4(1),x4(2),x4(3),x4(4)))* &
            U(mu,x4(1),x4(2),x4(3),x4(4))* &
            U(nu,x5(1),x5(2),x5(3),x5(4))
    end do
    staples%mat = sum%mat
  end function staples

  function Z(U,x,mu)
    type(su3alg) :: Z
    type(su3), dimension(:,:,:,:,:), intent(in) :: U
    integer, intent(in) :: x(4), mu
    type(matrix3x3) :: Z1
    Z1 = -TA(U(mu,x(1),x(2),x(3),x(4)) * dagger(staples(U,x,mu)) )
    Z%mat = Z1%mat
  end function Z
  
    
end module gauge
