module su3facts
  implicit none
  integer, parameter, private :: dp = 8
  complex(dp), parameter, private :: i = (0.0_dp,1.0_dp)
 
  type matrix4x4
     complex(dp), dimension(4,4) :: mat
  end type matrix4x4

  type :: matrix3x3
     complex(dp) :: mat(3,3)
   contains
     procedure :: init
     procedure :: zero
     procedure :: TA
  end type matrix3x3

  type, extends(matrix3x3) :: su3
   contains
     procedure :: normalize
     procedure :: tr => tr_su3
     procedure :: dagger => dagger_su3
     procedure :: init_su3
  end type su3

  type, extends(matrix3x3) :: su3alg
     real(dp) :: r(8)
   contains
     procedure :: init_su3alg
     procedure :: tr => tr_su3alg
     procedure :: dagger => dagger_su3alg
  end type su3alg

  type(su3alg) :: gellmann_matrix(8)
  type(matrix4x4) :: dirac_matrix(5)
  complex(dp), dimension(3,3) :: delta_3x3 = &
  reshape([(1.0_dp,0.0_dp),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),&
           (0.0_dp,0.0_dp),(1.0_dp,0.0_dp),(0.0_dp,0.0_dp),&
           (0.0_dp,0.0_dp),(0.0_dp,0.0_dp),(1.0_dp,0.0_dp)],[3,3])
                                           

  interface exp
     module procedure :: exp_su3alg
  end interface exp
  
  interface dagger
     module procedure :: dagger_mat3x3
     module procedure :: dagger_su3
     module procedure :: dagger_su3alg
  end interface dagger  
  
  interface tr
     module procedure :: tr_mat3x3
     module procedure :: tr_su3
     module procedure :: tr_su3alg
  end interface tr
  
  interface operator(+)
     module procedure :: matrix4x4_sum
     module procedure :: matrix3x3_sum
     module procedure :: su3_sum
     module procedure :: su3alg_sum
  end interface operator(+)
  
  interface operator(-)
     module procedure :: matrix4x4_substraction
     module procedure :: matrix4x4_substraction2
     module procedure :: matrix3x3_substraction
     module procedure :: matrix3x3_substraction2
     module procedure :: su3_substraction
     module procedure :: su3_substraction2
     module procedure :: su3alg_substraction
     module procedure :: su3alg_substraction2
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
     module procedure :: su3_product
     module procedure :: su3_product_scalar_real_left
     module procedure :: su3_product_scalar_real_right
     module procedure :: su3_product_scalar_complex_left
     module procedure :: su3_product_scalar_complex_right
     module procedure :: su3alg_product
     module procedure :: su3alg_product_scalar_real_left
     module procedure :: su3alg_product_scalar_real_right
     module procedure :: su3alg_product_scalar_complex_left
     module procedure :: su3alg_product_scalar_complex_right
     module procedure :: product_hybrid
  end interface operator(*)
  
contains

  elemental pure function matrix4x4_sum(A,B)
    type(matrix4x4), intent(in) :: A, B
    type(matrix4x4) :: matrix4x4_sum
    matrix4x4_sum%mat = A%mat+ B%mat
  end function matrix4x4_sum

  elemental pure function matrix4x4_substraction(A,B)
    type(matrix4x4), intent(in) :: A, B
    type(matrix4x4) :: matrix4x4_substraction
    matrix4x4_substraction%mat = A%mat- B%mat
  end function matrix4x4_substraction
  
  elemental pure function matrix4x4_substraction2(A)
    type(matrix4x4), intent(in) :: A
    type(matrix4x4) :: matrix4x4_substraction2
    matrix4x4_substraction2%mat = -A%mat
  end function matrix4x4_substraction2
  
  elemental pure function matrix4x4_product(A,B)
    type(matrix4x4), intent(in) :: A, B
    type(matrix4x4) :: matrix4x4_product
    matrix4x4_product%mat = matmul(A%mat,B%mat)
  end function matrix4x4_product
  
  elemental pure function matrix4x4_product_scalar_real_left(a,B)
    type(matrix4x4), intent(in) :: B
    real(dp), intent(in) :: a
    type(matrix4x4) :: matrix4x4_product_scalar_real_left
    matrix4x4_product_scalar_real_left%mat = a*B%mat
  end function matrix4x4_product_scalar_real_left

  elemental pure function matrix4x4_product_scalar_real_right(B,a)
    type(matrix4x4), intent(in) :: B
    real(dp), intent(in) :: a
    type(matrix4x4) :: matrix4x4_product_scalar_real_right
    matrix4x4_product_scalar_real_right%mat = B%mat*a
  end function matrix4x4_product_scalar_real_right

  elemental pure function matrix4x4_product_scalar_complex_left(a,B)
    type(matrix4x4), intent(in) :: B
    complex(dp), intent(in) :: a
    type(matrix4x4) :: matrix4x4_product_scalar_complex_left
    matrix4x4_product_scalar_complex_left%mat = a*B%mat
  end function matrix4x4_product_scalar_complex_left

  elemental pure function matrix4x4_product_scalar_complex_right(B,a)
    type(matrix4x4), intent(in) :: B
    complex(dp), intent(in) :: a
    type(matrix4x4) :: matrix4x4_product_scalar_complex_right
    matrix4x4_product_scalar_complex_right%mat = B%mat*a
  end function matrix4x4_product_scalar_complex_right
 
  elemental pure function matrix3x3_sum(A,B)
    type(matrix3x3), intent(in) :: A, B
    type(matrix3x3) :: matrix3x3_sum
    matrix3x3_sum%mat = A%mat+ B%mat
  end function matrix3x3_sum

  elemental pure function matrix3x3_substraction(A,B)
    type(matrix3x3), intent(in) :: A, B
    type(matrix3x3) :: matrix3x3_substraction
    matrix3x3_substraction%mat = A%mat- B%mat
  end function matrix3x3_substraction
  
  elemental pure function matrix3x3_substraction2(A)
    type(matrix3x3), intent(in) :: A
    type(matrix3x3) :: matrix3x3_substraction2
    matrix3x3_substraction2%mat = -A%mat
  end function matrix3x3_substraction2
  
  elemental pure function matrix3x3_product(A,B)
    type(matrix3x3), intent(in) :: A, B
    type(matrix3x3) :: matrix3x3_product
    matrix3x3_product%mat = matmul(A%mat,B%mat)
  end function matrix3x3_product
  
  elemental pure function matrix3x3_product_scalar_real_left(a,B)
    type(matrix3x3), intent(in) :: B
    real(dp), intent(in) :: a
    type(matrix3x3) :: matrix3x3_product_scalar_real_left
    matrix3x3_product_scalar_real_left%mat = a*B%mat
  end function matrix3x3_product_scalar_real_left

  elemental pure function matrix3x3_product_scalar_real_right(B,a)
    type(matrix3x3), intent(in) :: B
    real(dp), intent(in) :: a
    type(matrix3x3) :: matrix3x3_product_scalar_real_right
    matrix3x3_product_scalar_real_right%mat = B%mat*a
  end function matrix3x3_product_scalar_real_right

  elemental pure function matrix3x3_product_scalar_complex_left(a,B)
    type(matrix3x3), intent(in) :: B
    complex(dp), intent(in) :: a
    type(matrix3x3) :: matrix3x3_product_scalar_complex_left
    matrix3x3_product_scalar_complex_left%mat = a*B%mat
  end function matrix3x3_product_scalar_complex_left

  elemental pure function matrix3x3_product_scalar_complex_right(B,a)
    type(matrix3x3), intent(in) :: B
    complex(dp), intent(in) :: a
    type(matrix3x3) :: matrix3x3_product_scalar_complex_right
    matrix3x3_product_scalar_complex_right%mat = B%mat*a
  end function matrix3x3_product_scalar_complex_right

  elemental pure subroutine init(U)
    class(matrix3x3), intent(inout) :: U
    U%mat = delta_3x3
  end subroutine init
  
  pure elemental subroutine zero(U)
    class(matrix3x3), intent(inout) :: U
    U%mat = 0.0_dp
  end subroutine zero
  
  elemental pure function dagger_mat3x3(U)
    type(matrix3x3) :: dagger_mat3x3
    type(matrix3x3), intent(in) :: U
    dagger_mat3x3%mat = transpose(conjg(U%mat))
  end function dagger_mat3x3

  elemental pure function tr_mat3x3(U)
    complex(dp) :: tr_mat3x3
    type(matrix3x3), intent(in) :: U
    tr_mat3x3 = U%mat(1,1) +  U%mat(2,2) +  U%mat(3,3)
  end function tr_mat3x3
  
  elemental pure function su3_sum(A,B)
    type(su3), intent(in) :: A, B
    type(su3) :: su3_sum
    su3_sum%mat = A%mat+ B%mat
  end function su3_sum

  elemental pure function su3_substraction(A,B)
    type(su3), intent(in) :: A, B
    type(su3) :: su3_substraction
    su3_substraction%mat = A%mat- B%mat
  end function su3_substraction
  
  elemental pure function su3_substraction2(A)
    type(su3), intent(in) :: A
    type(su3) :: su3_substraction2
    su3_substraction2%mat = -A%mat
  end function su3_substraction2
  
  elemental pure function su3_product(A,B)
    type(su3), intent(in) :: A, B
    type(su3) :: su3_product
    su3_product%mat = matmul(A%mat,B%mat)
  end function su3_product
  
  elemental pure function su3_product_scalar_real_left(a,B)
    type(su3), intent(in) :: B
    real(dp), intent(in) :: a
    type(su3) :: su3_product_scalar_real_left
    su3_product_scalar_real_left%mat = a*B%mat
  end function su3_product_scalar_real_left

  elemental pure function su3_product_scalar_real_right(B,a)
    type(su3), intent(in) :: B
    real(dp), intent(in) :: a
    type(su3) :: su3_product_scalar_real_right
    su3_product_scalar_real_right%mat = B%mat*a
  end function su3_product_scalar_real_right

  elemental pure function su3_product_scalar_complex_left(a,B)
    type(su3), intent(in) :: B
    complex(dp), intent(in) :: a
    type(su3) :: su3_product_scalar_complex_left
    su3_product_scalar_complex_left%mat = a*B%mat
  end function su3_product_scalar_complex_left

  elemental pure function su3_product_scalar_complex_right(B,a)
    type(su3), intent(in) :: B
    complex(dp), intent(in) :: a
    type(su3) :: su3_product_scalar_complex_right
    su3_product_scalar_complex_right%mat = B%mat*a
  end function su3_product_scalar_complex_right

  elemental function normalize(U)
    type(su3) :: normalize
    class(su3), intent(in) :: U
    normalize%mat = delta_3x3
  end function normalize

  elemental pure function dagger_su3(U)
    type(su3) :: dagger_su3
    class(su3), intent(in) :: U
    dagger_su3%mat = transpose(conjg(U%mat))
  end function dagger_su3

  elemental pure function tr_su3(U)
    complex(dp) :: tr_su3
    class(su3), intent(in) :: U
    tr_su3 = U%mat(1,1) +  U%mat(2,2) +  U%mat(3,3)
  end function tr_su3






 

  
  
  
  
  subroutine create_kronecker_delta(delta,n)
    integer, intent(in) :: n
    complex(dp) :: delta(n,n)
    integer :: k
    
    delta = (0.0_dp,0.0_dp)
    do k = 1, n
       delta(k,k) = (1.0_dp,0.0_dp)
    end do
    
    
  end subroutine create_kronecker_delta

  subroutine create_gellmann_matrices()

    call gellmann_matrix%zero()
    
    !            [0 1 0]
    ! lambda_1 = [1 0 0]
    !            [0 0 0]
    gellmann_matrix(1)%mat(2,1) = 1.0_dp
    gellmann_matrix(1)%mat(1,2) = 1.0_dp

    
    !            [ 0 i 0]
    ! lambda_2 = [-i 0 0]
    !            [ 0 0 0]
    gellmann_matrix(2)%mat(2,1) =-i
    gellmann_matrix(2)%mat(1,2) = i

    !            [1  0 0]
    ! lambda_3 = [0 -1 0]
    !            [0  0 0]
    gellmann_matrix(3)%mat(1,1) = 1.0_dp
    gellmann_matrix(3)%mat(2,2) =-1.0_dp

    !            [0 0 1]
    ! lambda_4 = [0 0 0]
    !            [1 0 0]
    gellmann_matrix(4)%mat(1,3) = 1.0_dp
    gellmann_matrix(4)%mat(3,1) = 1.0_dp

    !            [0 0 -i]
    ! lambda_5 = [0 0  0]
    !            [i 0  0]
    gellmann_matrix(5)%mat(1,3) =-i
    gellmann_matrix(5)%mat(3,1) = i

    !            [0 0 0]
    ! lambda_6 = [0 0 1]
    !            [0 1 0]
    gellmann_matrix(6)%mat(3,2) = 1.0_dp
    gellmann_matrix(6)%mat(2,3) = 1.0_dp

    !            [0 0  0]
    ! lambda_7 = [0 0 -i]
    !            [0 i  0]
    gellmann_matrix(7)%mat(2,3) =-i
    gellmann_matrix(7)%mat(3,2) = i

    !                      [1 0  0]
    ! lambda_8 = 1/sqrt(3) [0 1  0]
    !                      [0 0 -2]
    gellmann_matrix(8)%mat(1,1) =  1.0_dp/sqrt(3.0_dp)
    gellmann_matrix(8)%mat(2,2) =  1.0_dp/sqrt(3.0_dp)
    gellmann_matrix(8)%mat(3,3) = -2.0_dp/sqrt(3.0_dp)
    
  end subroutine create_gellmann_matrices
  
  ! Create dirac_matrix matrices in Euclidean space-time
  subroutine create_gamma_matrices()
    integer :: mu, m, n
    do mu = 1, 5
       do n = 1, 4
          do m = 1, 4
             dirac_matrix(mu)%mat(n,m) = (0.0_dp,0.0_dp)
          end do
       end do
    end do
        
    !           [0 0  0 -i]
    ! gamma_1 = [0 0 -i  0]
    !           [0 i  0  0]
    !           [i 0  0  0]
    dirac_matrix(1)%mat(4,1) = i
    dirac_matrix(1)%mat(3,2) = i
    dirac_matrix(1)%mat(2,3) = -i
    dirac_matrix(1)%mat(1,4) = -i
    
    !           [ 0 0 0 -1]
    ! gamma_2 = [ 0 0 1  0]
    !           [ 0 1 0  0]
    !           [-1 0 0  0]
    dirac_matrix(2)%mat(4,1) = -1.0_dp
    dirac_matrix(2)%mat(3,2) = 1.0_dp
    dirac_matrix(2)%mat(2,3) = 1.0_dp
    dirac_matrix(2)%mat(1,4) = -1.0_dp
    
    !           [0  0 -i 0]
    ! gamma_3 = [0  0  0 i]
    !           [i  0  0 0]
    !           [0 -i  0 0]
    dirac_matrix(3)%mat(3,1) = i
    dirac_matrix(3)%mat(4,2) = -i
    dirac_matrix(3)%mat(1,3) = -i
    dirac_matrix(3)%mat(2,4) = i

    !           [0 0 1 0]
    ! gamma_4 = [0 0 0 1]
    !           [1 0 0 0]
    !           [0 1 0 0]
    dirac_matrix(4)%mat(3,1) = 1.0_dp
    dirac_matrix(4)%mat(4,2) = 1.0_dp
    dirac_matrix(4)%mat(1,3) = 1.0_dp
    dirac_matrix(4)%mat(2,4) = 1.0_dp


    !           [1 0  0  0]
    ! gamma_5 = [0 1  0  1]
    !           [0 0 -1  0]
    !           [0 0  0 -1]
    dirac_matrix(5)%mat(1,1) = 1.0_dp
    dirac_matrix(5)%mat(2,2) = 1.0_dp
    dirac_matrix(5)%mat(3,3) = -1.0_dp
    dirac_matrix(5)%mat(4,4) = -1.0_dp
  end subroutine create_gamma_matrices

  pure elemental function exp_su3alg(X) result(expX)
    ! Computes exp(X) for X in the su(3) Lie algebra using Cayley-Hamilton theorem.
    !
    ! By Cayley-Hamilton, any X in su(3) (traceless, antihermitian) satisfies:
    !   X^3 = -t*X + d*I
    ! where:
    !   t = -1/2 * tr(X^2)
    !   d =  det(X)
    !
    ! Therefore exp(X) = q0*I + q1*X + q2*X^2, with coefficients q0,q1,q2
    ! computed via Horner's method on the Taylor series:
    !   exp(X) = sum_{n=0}^{inf} X^n / n!
    !
    ! Recurrence (running from n=K down to n=0):
    !   q2_new = q1_old
    !   q1_new = q0_old - t * q2_old
    !   q0_new = 1/(i+1)! + d * q2_old
    type(su3alg), intent(in) :: X
    type(su3) :: expX
    type(su3alg) :: B, d3x3
    integer, parameter :: K = 20
    complex(dp) :: q0, q1, q2, q0old, q1old, q2old
    complex(dp) :: d, t, trB
    integer :: l
    
    ! X^2 and its trace (needed for t)
    B  = X*X
    trB = tr(B)
    d3x3%mat = delta_3x3
    ! Cayley-Hamilton coefficients
    ! For traceless X: characteristic poly is lambda^3 + t*lambda - d = 0
    t = -0.5_dp * trB          ! coefficient of lambda (no ii factor needed)
    d = determinant(X%mat)      ! constant term: det(X), no extra ii factor

    !print*, gamma(1.0_dp)
    ! Horner recurrence initialised at step K
    q0old = 1.0_dp / gamma(1.0_dp*(K+1))
   
    q1old = (0.0_dp, 0.0_dp)
    q2old = (0.0_dp, 0.0_dp)

    do l = K-1, 0, -1
       q0 = 1.0_dp/ gamma(1.0_dp*(l+1)) + d * q2old   ! X^3 = -t*X + d*I  => +d (not -ii*d)
       q1 = q0old - t * q2old
       q2 = q1old
       q0old = q0
       q1old = q1
       q2old = q2
    end do

    expX%mat = q0*d3x3%mat + q1*X%mat + q2*B%mat

  end function exp_su3alg
  
  pure function determinant(a) result(det)
    complex(dp), dimension(3,3), intent(in) :: a
    complex(dp) :: det

    det = a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) + a(1,3)*a(2,1)*a(3,2) &
         -a(1,3)*a(2,2)*a(3,1) - a(1,1)*a(2,3)*a(3,2) - a(1,2)*a(2,1)*a(3,3)
   
  end function determinant

  elemental subroutine init_su3(U,r1,r2,r3,r4,r5,r6,r7,r8)
    class(su3), intent(inout) :: U
    real(dp), intent(in) :: r1,r2,r3,r4,r5,r6,r7,r8
    type(su3alg) :: A
    type(su3) :: C
    
    A =     r1*gellmann_matrix(1) + r2*gellmann_matrix(2) + r3*gellmann_matrix(3)
    A = A + r4*gellmann_matrix(4) + r5*gellmann_matrix(5) + r6*gellmann_matrix(6)
    A = A + r7*gellmann_matrix(7) + r8*gellmann_matrix(8)
    C = exp(0.5_dp*i*A)
    U%mat = C%mat
    
  end subroutine init_su3

  elemental subroutine init_su3alg(U,r1,r2,r3,r4,r5,r6,r7,r8)
    class(su3alg), intent(inout) :: U
    real(dp), intent(in) :: r1,r2,r3,r4,r5,r6,r7,r8
    type(su3alg) :: A

    A =     r1*gellmann_matrix(1) + r2*gellmann_matrix(2) + r3*gellmann_matrix(3)
    A = A + r4*gellmann_matrix(4) + r5*gellmann_matrix(5) + r6*gellmann_matrix(6)
    A = A + r7*gellmann_matrix(7) + r8*gellmann_matrix(8)
    U%mat = 0.5_dp*A%mat 
    
  end subroutine init_su3alg
  
  elemental pure function su3alg_sum(A,B)
    type(su3alg), intent(in) :: A, B
    type(su3alg) :: su3alg_sum
    su3alg_sum%mat = A%mat+ B%mat
  end function su3alg_sum

  elemental pure function su3alg_substraction(A,B)
    type(su3alg), intent(in) :: A, B
    type(su3alg) :: su3alg_substraction
    su3alg_substraction%mat = A%mat- B%mat
  end function su3alg_substraction
  
  elemental pure function su3alg_substraction2(A)
    type(su3alg), intent(in) :: A
    type(su3alg) :: su3alg_substraction2
    su3alg_substraction2%mat = -A%mat
  end function su3alg_substraction2
  
  elemental pure function su3alg_product(A,B)
    type(su3alg), intent(in) :: A, B
    type(su3alg) :: su3alg_product
    su3alg_product%mat = matmul(A%mat,B%mat)
  end function su3alg_product
  
  elemental pure function su3alg_product_scalar_real_left(a,B)
    type(su3alg), intent(in) :: B
    real(dp), intent(in) :: a
    type(su3alg) :: su3alg_product_scalar_real_left
    su3alg_product_scalar_real_left%mat = a*B%mat
  end function su3alg_product_scalar_real_left

  elemental pure function su3alg_product_scalar_real_right(B,a)
    type(su3alg), intent(in) :: B
    real(dp), intent(in) :: a
    type(su3alg) :: su3alg_product_scalar_real_right
    su3alg_product_scalar_real_right%mat = B%mat*a
  end function su3alg_product_scalar_real_right

  elemental pure function su3alg_product_scalar_complex_left(a,B)
    type(su3alg), intent(in) :: B
    complex(dp), intent(in) :: a
    type(su3alg) :: su3alg_product_scalar_complex_left
    su3alg_product_scalar_complex_left%mat = a*B%mat
  end function su3alg_product_scalar_complex_left

  elemental pure function su3alg_product_scalar_complex_right(B,a)
    type(su3alg), intent(in) :: B
    complex(dp), intent(in) :: a
    type(su3alg) :: su3alg_product_scalar_complex_right
    su3alg_product_scalar_complex_right%mat = B%mat*a
  end function su3alg_product_scalar_complex_right

  elemental pure function dagger_su3alg(U)
    type(su3alg) :: dagger_su3alg
    class(su3alg), intent(in) :: U
    dagger_su3alg%mat = transpose(conjg(U%mat))
  end function dagger_su3alg

  elemental pure function tr_su3alg(U)
    complex(dp) :: tr_su3alg
    class(su3alg), intent(in) :: U
    tr_su3alg = U%mat(1,1) +  U%mat(2,2) +  U%mat(3,3)
  end function tr_su3alg
  
  function TA(W)
    type(matrix3x3) :: TA
    class(matrix3x3), intent(in) :: W
    integer, dimension(4) :: x2, x3
    type(matrix3x3) :: A
    A = W - dagger(W)
    TA%mat = 0.5_dp*A%mat - (1/6.0_dp)*tr(A)*delta_3x3 
  end function TA

  elemental function product_hybrid(A,B)
    type(matrix3x3) :: product_hybrid
    type(su3), intent(in) :: A
    type(matrix3x3), intent(in) :: B
    product_hybrid%mat = matmul(A%mat,B%mat)
  end function product_hybrid
  
end module su3facts
