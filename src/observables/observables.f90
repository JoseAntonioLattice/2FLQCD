module observables
  
  use su3facts
  use parameters, only : L, Lt, Lx, Ly, Lz
  use gauge, only : plaquette
  implicit none

  integer, parameter, private :: dp = 8
  
contains

  function plaquette_value(U) result(plqv)
    use gauge, only : plaquette
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp) :: plqv
    integer :: t,x,y,z, mu, nu

    plqv = 0.0_dp
    do t = 1, Lt
       do x = 1, Lx
          do y = 1, Ly
             do z = 1, Lz
                do mu = 1, 3
                   do nu = mu + 1, 4
                      plqv = plqv + real(tr(plaquette(U,[t,x,y,z],mu,nu)))
                   end do
                end do
             end do
          end do
       end do
    end do
    
    plqv = plqv/(3*6.0_dp*product(L))
             
  end function plaquette_value


  function action(U,beta) result(S)
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp),intent(in) :: beta
    real(dp) :: S
    integer :: t,x,y,z, mu, nu

    S = 0.0_dp
    do t = 1, Lt
       do x = 1, Lx
          do y = 1, Ly
             do z = 1, Lz
                do mu = 1, 3
                   do nu = mu + 1, 4
                      S = S + real(tr(plaquette(U,[t,x,y,z],mu,nu)))
                   end do
                end do
             end do
          end do
       end do
    end do
    
    S  = beta/3.0_dp*(18.0_dp*product(L) - S)
             
  end function action

  function energy_density(U) result(E)
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp) :: E
    integer :: t,x,y,z, mu, nu

    E = 0.0_dp
    do t = 1, Lt
       do x = 1, Lx
          do y = 1, Ly
             do z = 1, Lz
                do mu = 1, 3
                   do nu = mu + 1, 4
                      E = E + real(tr(plaquette(U,[t,x,y,z],mu,nu)))
                   end do
                end do
             end do
          end do
       end do
    end do

    E = 2*(18.0_dp - E/product(L))
  end function energy_density
  
end module observables

