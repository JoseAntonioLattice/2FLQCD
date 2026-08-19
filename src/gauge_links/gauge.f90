module gauge
  use su3facts 
  use pbc, only : ip, im
  use parameters, only : Lt,lx,ly,lz
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

  function plaquette_oriented(U,x,mu,nu)
    type(su3) :: plaquette_oriented
    type(su3), dimension(:,:,:,:,:), intent(in) :: U
    type(su3) :: U1, U2, U3, U4
    integer, intent(in) :: x(4), mu, nu
    integer, dimension(4) :: xnew

    xnew = x
   
    U1 = oriented_link(U,xnew,mu)
    U2 = oriented_link(U,xnew,nu)
    U3 = oriented_link(U,xnew,-mu)
    U4 = oriented_link(U,xnew,-nu)

    plaquette_oriented = U1*U2*U3*U4
  end function plaquette_oriented

  function oriented_link(U,x,mu)
    type(su3), intent(in) :: U(:,:,:,:,:)
    integer, intent(inout) :: x(4)
    integer, intent(in) :: mu
    type(su3) :: oriented_link

    if( mu < 0 )then
       x = im(x,-mu)
       oriented_link = dagger(U(-mu,x(1),x(2),x(3),x(4)))
    else
       oriented_link = U(mu,x(1),x(2),x(3),x(4))
       x = ip(x,mu)
    end if
    
  end function oriented_link
  
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

  function clover(U,x,mu,nu)
    type(su3), dimension(:,:,:,:,:), intent(in) :: U
    integer, intent(in) :: x(4), mu, nu
    type(su3) :: tmp
    type(matrix3x3) :: clover
    integer, dimension(4) :: ipx_mu, ipx_nu,imx_mu, imx_nu, &
         x_im_mu_ip_nu, x_im_mu_im_nu,  x_ip_mu_im_nu


    ipx_mu = ip(x,mu)
    ipx_nu = ip(x,nu)
    imx_mu = im(x,mu)
    imx_nu = im(x,nu)

    x_im_mu_ip_nu = im(ipx_nu,mu)
    x_im_mu_im_nu = im(imx_nu,mu)
    x_ip_mu_im_nu = ip(imx_nu,mu)

    
    tmp = U(mu,x(1),x(2),x(3),x(4)) * &
          U(nu,ipx_mu(1),ipx_mu(2),ipx_mu(3),ipx_mu(4)) * &
          dagger(U(mu,ipx_nu(1),ipx_nu(2),ipx_nu(3),ipx_nu(4))) *&
          dagger(U(nu,x(1),x(2),x(3),x(4))) + &
          U(nu,x(1),x(2),x(3),x(4)) * &
          dagger(U(mu,x_im_mu_ip_nu(1),x_im_mu_ip_nu(2),x_im_mu_ip_nu(3),x_im_mu_ip_nu(4))) * &
          dagger(U(nu,imx_mu(1),imx_mu(2),imx_mu(3),imx_mu(4))) &
          * U(mu,imx_mu(1),imx_mu(2),imx_mu(3),imx_mu(4)) + &
          dagger(U(mu,imx_mu(1),imx_mu(2),imx_mu(3),imx_mu(4))) &
          * dagger(U(nu,x_im_mu_im_nu(1),x_im_mu_im_nu(2),x_im_mu_im_nu(3),x_im_mu_im_nu(4))) * &
          U(mu,x_im_mu_im_nu(1),x_im_mu_im_nu(2),x_im_mu_im_nu(3),x_im_mu_im_nu(4)) &
          * U(nu,imx_nu(1),imx_nu(2),imx_nu(3),imx_nu(4)) + &
          dagger(U(nu,imx_nu(1),imx_nu(2),imx_nu(3),imx_nu(4)))&
          * U(mu,imx_nu(1),imx_nu(2),imx_nu(3),imx_nu(4)) * &
          U(nu,x_ip_mu_im_nu(1),x_ip_mu_im_nu(2),x_ip_mu_im_nu(3),x_ip_mu_im_nu(4)) &
          * dagger(U(mu,x(1),x(2),x(3),x(4)))

    clover%mat = tmp%mat
    
  end function clover

  function clover2(U,x,mu,nu)
    type(su3), dimension(:,:,:,:,:), intent(in) :: U
    integer, intent(in) :: x(4), mu, nu
    type(matrix3x3) :: clover2
    type(su3) :: tmp
    
    tmp = plaquette_oriented(U,x, mu, nu) + plaquette_oriented(U,x, nu,-mu) + &
          plaquette_oriented(U,x,-mu,-nu) + plaquette_oriented(U,x,-nu, mu) 

    clover2%mat = tmp%mat
    
  end function clover2
    
end module gauge
