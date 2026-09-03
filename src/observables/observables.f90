module observables
  
  use su3facts
  use parameters, only : L, Lt, Lx, Ly, Lz
  use gauge
  implicit none

  integer, parameter, private :: dp = 8
  real(dp), parameter, private :: pi = acos(-1.0_dp)
  
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


  !!Computation of energy density with clover definition
  ! E_{mu,nu} = 1/4 G_{mu,nu}^a G_{mu,nu}^a = -1/2 Tr(G_{mu,nu} G_{mu,nu}) = -1/128 tr((Q_{mu,nu} - dagger(Q_{mu,nu}))^2)
  ! G_{mu,nu} = (Q_{mu,nu} - dagger(Q_{mu,nu}))/8
  ! Q_{mu,nu}: clover
  ! E = \sum_{mu,nu} E_{mu,nu} = 2 \sum_{mu<nu} E_{mu,nu} = -1/64 \sum_{mu<nu} tr((Q_{mu,nu} - dagger(Q_{mu,nu}))^2)
  function energy_density_clover(U) result(E)
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp) :: E
    integer :: t,x,y,z, mu, nu
    type(matrix3x3) :: Q, tmp
    
    E = 0.0_dp
    do t = 1, Lt
       do x = 1, Lx
          do y = 1, Ly
             do z = 1, Lz
                do mu = 1, 3
                   do nu = mu+1, 4
                      Q = clover(U,[t,x,y,z],mu,nu)
                      tmp = Q - dagger(Q)
                      !tmp%mat = tmp%mat -tr(tmp)/3.0_dp*delta_3x3
                      E = E + real(tr(tmp*tmp))
                   end do
                end do
             end do
          end do
       end do
    end do

    E = -E/(64.0_dp*product(L))
  end function energy_density_clover

  
  function energy_density_clover2(U) result(E)
    type(su3), dimension(4,Lt,Lx,Ly,Lz), intent(inout) :: U
    real(dp) :: E
    integer :: t,x,y,z, mu, nu
    type(matrix3x3) :: Q, tmp

    E = 0.0_dp
    do t = 1, Lt
       do x = 1, Lx
          do y = 1, Ly
             do z = 1, Lz
                do mu = 1, 3
                   do nu = mu+1, 4
                      Q = clover2(U,[t,x,y,z],mu,nu)
                      tmp = Q - dagger(Q)
                      E = E + tr(tmp*tmp)
                   end do
                end do
             end do
          end do
       end do
    end do

    E = -E/(64.0_dp*product(L))
  end function energy_density_clover2

  function unormalized_topological_charge_density_clover(U) result(q)
    type(su3), dimension(4,Lt,Lx,Ly,Lz) :: U
    integer :: t, x, y, z
    type(matrix3x3) :: Q34, Q24, Q23
    real(dp) :: q(Lt,Lx,Ly,Lz)
    
    do t = 1, Lt
       do x = 1, Lx
          do y = 1, Ly
             do z = 1, Lz
                Q34 = clover(U,[t,x,y,z],3,4)
                Q34 = Q34 - dagger(Q34)

                Q24 = clover(U,[t,x,y,z],2,4)
                Q24 = Q24 - dagger(Q24)

                Q23 = clover(U,[t,x,y,z],2,3)
                Q23 = Q23 - dagger(Q23)
        
                q(t,x,y,z) = real(tr(clover(U,[t,x,y,z],1,2)*Q34) &
                                - tr(clover(U,[t,x,y,z],1,3)*Q24) &
                                + tr(clover(U,[t,x,y,z],1,4)*Q23))
             end do
          end do
       end do
    end do
        
  end function unormalized_topological_charge_density_clover

  function topological_charge_clover(U) result(q)
    type(su3), dimension(4,Lt,Lx,Ly,Lz) :: U
    integer :: t, x, y, z
    type(matrix3x3) :: Q34, Q24, Q23
    real(dp) :: q

    q = -1/(128*pi**2)*sum(unormalized_topological_charge_density_clover(U))
    
  end function topological_charge_clover
  
end module observables

