module su3facts
  implicit none
  integer, parameter, private :: dp = 8
  complex(dp), parameter, private :: i = (0.0_dp,1.0_dp)
  complex(dp), dimension(5,4,4) :: gamma
  real(dp), dimension(4,4) :: delta
contains


  subroutine create_kronecker_delta()
    delta = 0.0_dp
    delta(1,1) = 1.0_dp
    delta(2,2) = 1.0_dp
    delta(3,3) = 1.0_dp
    delta(4,4) = 1.0_dp
  end subroutine create_kronecker_delta
  
  ! Create gamma matrices in Euclidean space-time
  subroutine create_gamma_matrices()

    gamma = (0.0_dp,0.0_dp)
    
               
    !           [0 0  0 -i]
    ! gamma_1 = [0 0 -i  0]
    !           [0 i  0  0]
    !           [i 0  0  0]
    gamma(1,4,1) = i
    gamma(1,3,2) = i
    gamma(1,2,3) = -i
    gamma(1,1,4) = -i
    
    !           [ 0 0 0 -1]
    ! gamma_2 = [ 0 0 1  0]
    !           [ 0 1 0  0]
    !           [-1 0 0  0]
    gamma(2,4,1) = -1.0_dp
    gamma(2,3,2) = 1.0_dp
    gamma(2,2,3) = 1.0_dp
    gamma(2,1,4) = -1.0_dp
    
    !           [0  0 -i 0]
    ! gamma_3 = [0  0  0 i]
    !           [i  0  0 0]
    !           [0 -i  0 0]
    gamma(3,3,1) = i
    gamma(3,4,2) = -i
    gamma(3,1,3) = -i
    gamma(3,2,4) = i

    !           [0 0 1 0]
    ! gamma_4 = [0 0 0 1]
    !           [1 0 0 0]
    !           [0 1 0 0]
    gamma(4,3,1) = 1.0_dp
    gamma(4,4,2) = 1.0_dp
    gamma(4,1,3) = 1.0_dp
    gamma(4,2,4) = 1.0_dp

    
    
  end subroutine create_gamma_matrices
  
end module su3facts
