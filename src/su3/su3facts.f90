module su3facts
  implicit none
  integer, parameter, private :: dp = 8
  complex(dp), parameter, private :: i = (0.0_dp,1.0_dp)
 
  type matrix4x4
     complex(dp), dimension(4,4) :: mat
  end type matrix4x4

  type :: matrix3x3
     complex(dp) :: mat(3,3)
  end type matrix3x3

  type, extends(matrix3x3) :: su3

  end type su3

  type(matrix4x4), dimension(5) :: gamma
  real(dp), dimension(4,4) :: delta
  complex(dp), dimension(3,3) :: delta_3x3
  
  interface operator(+)
     module procedure :: matrix4x4_sum
     module procedure :: matrix3x3_sum
  end interface operator(+)
  
  interface operator(-)
     module procedure :: matrix4x4_substraction
     module procedure :: matrix4x4_substraction2
     module procedure :: matrix3x3_substraction
     module procedure :: matrix3x3_substraction2
  end interface operator(-)
  
  interface operator(*)
     module procedure :: matrix4x4_product
     module procedure :: matrix4x4_product_scalar_real_left
     module procedure :: matrix4x4_product_scalar_real_right
     module procedure :: matrix4x4_product_scalar_complex_left
     module procedure :: matrix4x4_product_scalar_complex_right
     module procedure :: matrix3x3_product
     module procedure :: matrix3x3_product_scalar_real_left
     module procedure :: matrix3x3_product_scalar_real_right
     module procedure :: matrix3x3_product_scalar_complex_left
     module procedure :: matrix3x3_product_scalar_complex_right
  end interface operator(*)
  
contains

  function matrix4x4_sum(A,B)
    type(matrix4x4), intent(in) :: A, B
    type(matrix4x4) :: matrix4x4_sum
    matrix4x4_sum%mat = A%mat+ B%mat
  end function matrix4x4_sum

  function matrix4x4_substraction(A,B)
    type(matrix4x4), intent(in) :: A, B
    type(matrix4x4) :: matrix4x4_substraction
    matrix4x4_substraction%mat = A%mat- B%mat
  end function matrix4x4_substraction
  
  function matrix4x4_substraction2(A)
    type(matrix4x4), intent(in) :: A
    type(matrix4x4) :: matrix4x4_substraction2
    matrix4x4_substraction2%mat = -A%mat
  end function matrix4x4_substraction2
  
  function matrix4x4_product(A,B)
    type(matrix4x4), intent(in) :: A, B
    type(matrix4x4) :: matrix4x4_product
    matrix4x4_product%mat = matmul(A%mat,B%mat)
  end function matrix4x4_product
  
  function matrix4x4_product_scalar_real_left(a,B)
    type(matrix4x4), intent(in) :: B
    real(dp), intent(in) :: a
    type(matrix4x4) :: matrix4x4_product_scalar_real_left
    matrix4x4_product_scalar_real_left%mat = a*B%mat
  end function matrix4x4_product_scalar_real_left

  function matrix4x4_product_scalar_real_right(B,a)
    type(matrix4x4), intent(in) :: B
    real(dp), intent(in) :: a
    type(matrix4x4) :: matrix4x4_product_scalar_real_right
    matrix4x4_product_scalar_real_right%mat = B%mat*a
  end function matrix4x4_product_scalar_real_right

  function matrix4x4_product_scalar_complex_left(a,B)
    type(matrix4x4), intent(in) :: B
    complex(dp), intent(in) :: a
    type(matrix4x4) :: matrix4x4_product_scalar_complex_left
    matrix4x4_product_scalar_complex_left%mat = a*B%mat
  end function matrix4x4_product_scalar_complex_left

  function matrix4x4_product_scalar_complex_right(B,a)
    type(matrix4x4), intent(in) :: B
    complex(dp), intent(in) :: a
    type(matrix4x4) :: matrix4x4_product_scalar_complex_right
    matrix4x4_product_scalar_complex_right%mat = B%mat*a
  end function matrix4x4_product_scalar_complex_right
  

  function matrix3x3_sum(A,B)
    type(matrix3x3), intent(in) :: A, B
    type(matrix3x3) :: matrix3x3_sum
    matrix3x3_sum%mat = A%mat+ B%mat
  end function matrix3x3_sum

  function matrix3x3_substraction(A,B)
    type(matrix3x3), intent(in) :: A, B
    type(matrix3x3) :: matrix3x3_substraction
    matrix3x3_substraction%mat = A%mat- B%mat
  end function matrix3x3_substraction
  
  function matrix3x3_substraction2(A)
    type(matrix3x3), intent(in) :: A
    type(matrix3x3) :: matrix3x3_substraction2
    matrix3x3_substraction2%mat = -A%mat
  end function matrix3x3_substraction2
  
  function matrix3x3_product(A,B)
    type(matrix3x3), intent(in) :: A, B
    type(matrix3x3) :: matrix3x3_product
    matrix3x3_product%mat = matmul(A%mat,B%mat)
  end function matrix3x3_product
  
  function matrix3x3_product_scalar_real_left(a,B)
    type(matrix3x3), intent(in) :: B
    real(dp), intent(in) :: a
    type(matrix3x3) :: matrix3x3_product_scalar_real_left
    matrix3x3_product_scalar_real_left%mat = a*B%mat
  end function matrix3x3_product_scalar_real_left

  function matrix3x3_product_scalar_real_right(B,a)
    type(matrix3x3), intent(in) :: B
    real(dp), intent(in) :: a
    type(matrix3x3) :: matrix3x3_product_scalar_real_right
    matrix3x3_product_scalar_real_right%mat = B%mat*a
  end function matrix3x3_product_scalar_real_right

  function matrix3x3_product_scalar_complex_left(a,B)
    type(matrix3x3), intent(in) :: B
    complex(dp), intent(in) :: a
    type(matrix3x3) :: matrix3x3_product_scalar_complex_left
    matrix3x3_product_scalar_complex_left%mat = a*B%mat
  end function matrix3x3_product_scalar_complex_left

  function matrix3x3_product_scalar_complex_right(B,a)
    type(matrix3x3), intent(in) :: B
    complex(dp), intent(in) :: a
    type(matrix3x3) :: matrix3x3_product_scalar_complex_right
    matrix3x3_product_scalar_complex_right%mat = B%mat*a
  end function matrix3x3_product_scalar_complex_right

  function normalize(U)
    type(matrix3x3) :: normalize
    type(matrix3x3) :: U
    normalize%mat = 1.0_dp
  end function normalize

  function dagger(U)
    type(matrix3x3) :: dagger
    type(matrix3x3) :: U
    dagger%mat = transpose(conjg(U%mat))
  end function dagger
  
  
  subroutine create_kronecker_delta()
    delta = 0.0_dp
    delta(1,1) = 1.0_dp
    delta(2,2) = 1.0_dp
    delta(3,3) = 1.0_dp
    delta(4,4) = 1.0_dp

    delta_3x3 = 0.0_dp
    delta_3x3(1,1) = 1.0_dp
    delta_3x3(2,2) = 1.0_dp
    delta_3x3(3,3) = 1.0_dp
    
  end subroutine create_kronecker_delta
  
  ! Create gamma matrices in Euclidean space-time
  subroutine create_gamma_matrices()
    integer :: mu, m, n
    do mu = 1, 5
       do n = 1, 4
          do m = 1, 4
             gamma(mu)%mat(n,m) = (0.0_dp,0.0_dp)
          end do
       end do
    end do
        
    !           [0 0  0 -i]
    ! gamma_1 = [0 0 -i  0]
    !           [0 i  0  0]
    !           [i 0  0  0]
    gamma(1)%mat(4,1) = i
    gamma(1)%mat(3,2) = i
    gamma(1)%mat(2,3) = -i
    gamma(1)%mat(1,4) = -i
    
    !           [ 0 0 0 -1]
    ! gamma_2 = [ 0 0 1  0]
    !           [ 0 1 0  0]
    !           [-1 0 0  0]
    gamma(2)%mat(4,1) = -1.0_dp
    gamma(2)%mat(3,2) = 1.0_dp
    gamma(2)%mat(2,3) = 1.0_dp
    gamma(2)%mat(1,4) = -1.0_dp
    
    !           [0  0 -i 0]
    ! gamma_3 = [0  0  0 i]
    !           [i  0  0 0]
    !           [0 -i  0 0]
    gamma(3)%mat(3,1) = i
    gamma(3)%mat(4,2) = -i
    gamma(3)%mat(1,3) = -i
    gamma(3)%mat(2,4) = i

    !           [0 0 1 0]
    ! gamma_4 = [0 0 0 1]
    !           [1 0 0 0]
    !           [0 1 0 0]
    gamma(4)%mat(3,1) = 1.0_dp
    gamma(4)%mat(4,2) = 1.0_dp
    gamma(4)%mat(1,3) = 1.0_dp
    gamma(4)%mat(2,4) = 1.0_dp


    !           [1 0  0  0]
    ! gamma_5 = [0 1  0  1]
    !           [0 0 -1  0]
    !           [0 0  0 -1]
    gamma(5)%mat(1,1) = 1.0_dp
    gamma(5)%mat(2,2) = 1.0_dp
    gamma(5)%mat(3,3) = -1.0_dp
    gamma(5)%mat(4,4) = -1.0_dp
  end subroutine create_gamma_matrices
  
end module su3facts
