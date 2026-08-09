program test
  use random
  use starts
  use pbc
  use su3facts, gamma => dirac_matrix
  implicit none
  integer, parameter :: dp = 8

  
  !call test_mean_zero()
  !call test_variance_one()
  !call test_gamma()
  !call test_gamma5()
  !call test_pbc()

  !call test_sum_4x4matrices
  !call test_substraction_4x4matrices()

  !call test_cold_start()

  !call test_su3_links

  !call test_gellmann_matrices
  !call test_unitarity()

  call test_assert_matrix()
  
contains
  subroutine assert_close(name, measured, expected, tol)
    character(*), intent(in) :: name
    real(dp), intent(in) :: measured, expected, tol

    if( abs(measured - expected) < tol ) then
       print"(2a)", " [PASS] : ", name
    else
       print"(3a,ES12.4,a,ES12.4,a,ES12.4)", " [FAIL] : ", name, " Expected:", expected, " Measured: ", measured, " Tol: ", tol
    end if
  end subroutine assert_close

  subroutine assert_equal_complex(name, measured, expected, tol)
    character(*), intent(in) :: name
    complex(dp), intent(in) :: measured, expected
    real(dp) :: tol

    if( abs(real(measured) - real(expected)) < tol .and. abs(aimag(measured) - aimag(expected)) < tol) then
       print"(2a)", " [PASS] : ", name
    else
       print"(3a,2ES12.4,a,2ES12.4,a,ES12.4)", " [FAIL] : ", name, " Expected:", expected, " Measured: ", measured, " Tol: ", tol
    end if
  end subroutine assert_equal_complex

  
  subroutine assert_Equal_complex_matrix(name, measured, expected,tol)
    character(*), intent(in) :: name
    complex(dp), intent(in) :: measured(:,:), expected(size(measured(:,1)),size(measured(1,:)))
    real(dp), intent(in) :: tol
    integer :: n, m, i, j

    n = size(measured(:,1))
    m = size(measured(1,:))

    do i = 1, n
       do j = 1, m
          call assert_equal_complex(name,measured(i,j), expected(i,j),tol)
       end do
    end do
    
    
  end subroutine assert_Equal_complex_matrix
  
  subroutine assert_equal(name, measured, expected)
    character(*), intent(in) :: name
    integer, intent(in) :: measured, expected

    if( measured == expected ) then
       print"(2a)", " [PASS] : ", name
    else
       print"(3a,i0,a,i0)", " [FAIL] : ", name, " Expected:", expected, " Measured: ", measured
    end if
  end subroutine assert_equal

  subroutine assert_equal_real_matrix(name, measured, expected,tol)
    character(*), intent(in) :: name
    real(dp), intent(in), dimension(:,:) :: measured, expected
    real(dp), intent(in) :: tol
    logical, dimension(size(measured(:,1)),size(measured(1,:))) :: condition
    integer :: i, j, n, m

    n = size(measured(:,1))
    m = size(measured(1,:))
    
    condition = abs(measured - expected) < tol
    if( all(condition) ) then
       print"(2a)", " [PASS] : ", name
    else
       print"(2a)", " [FAIL] : ", name
       do i = 1, n
          do j = 1, m
             if(condition(i,j) .eqv. .false.) &
                  print"(4x,a,i0,a,i0,a,ES12.4,a,ES12.4,a,ES12.4)", &
                  "(",i,",",j,") Expected: ", expected(i,j), " Measured: ", measured(i,j), " tol = ", tol
             
          end do
       end do
    end if
  end subroutine assert_equal_real_matrix

  subroutine test_assert_matrix()
    real(dp), dimension(2,2) :: A, B
    integer :: i, j

    A = 1.0_dp
    B = 1.0_dp
    B(2,2) = 2.0_dp
    B(1,2) = 0.0_dp
    call assert_equal_real_matrix("A = B (real matrix)",A,B,1.0e-12_dp)

  end subroutine test_assert_matrix
  
  subroutine assert_true(name, condition)
    character(*), intent(in) :: name
    logical, intent(in) :: condition
    
    if( condition ) then
       print"(2a)", " [PASS] : ", name
    else
       print"(2a)", " [FAIL] : ", name
    end if
  end subroutine assert_true

  subroutine test_mean_zero()
    real(dp), dimension(10**6) :: r, u1, u2

    call random_number(u1)
    call random_number(u2)
    call elemental_rgauss(r,u1,u2)
    call assert_close("Mean zero",avr(r),0.0_dp,5*std_err(r))
  end subroutine test_mean_zero

  subroutine test_variance_one()
    real(dp), dimension(10**6) :: r, u1, u2

    call random_number(u1)
    call random_number(u2)
    call elemental_rgauss(r,u1,u2)
    call assert_close("Variance one",avr(r),0.0_dp,5*var(r)*sqrt(2.0_dp/size(r)))
  end subroutine test_variance_one
  
  function avr(x)
    real(dp) :: avr
    real(dp) :: x(:)

    avr = sum(x)/size(x)

  end function avr

  function var(x)
    real(dp) :: var
    real(dp) :: x(:)

    var = sum((x - avr(x))**2)/(size(x) - 1)

  end function var

  function std_err(x)
    real(dp) :: std_err
    real(dp) :: x(:)

    std_err = sqrt(var(x)/size(x))
    
  end function std_err

  subroutine test_gamma()
    use su3facts,  gamma => dirac_matrix
    integer :: mu, nu, i, j
    type(matrix4x4) :: product
    logical :: condition
    
    call create_kronecker_delta()
    call create_gamma_matrices()
    
    do mu = 1, 4
       do nu = 1, 4
          product = gamma(mu)*gamma(nu) + gamma(nu)*gamma(mu)
          print*, "mu = ",mu, "nu = ", nu
          do i = 1, 4
             do j = 1, 4
                call assert_close("{gamma_mu,gamma_nu}", real(product%mat(i,j)),2*delta(mu,nu)*delta(i,j), 1.0E-12_dp) 
             end do
          end do
       end do
    end do
    
  end subroutine test_gamma

  subroutine test_gamma5()
    use su3facts,  gamma => dirac_matrix
    integer :: i, j
    type(matrix4x4) :: gamma5

    
    call create_kronecker_delta()
    call create_gamma_matrices()
    gamma5 = gamma(1)*gamma(2)*gamma(3)*gamma(4) 
    do i = 1, 4
       do j = 1, 4
          call assert_close("gamma5 real", real(gamma(5)%mat(i,j)), real(gamma5%mat(i,j)),1.0E-12_dp)
          call assert_close("gamma5 imag",aimag(gamma(5)%mat(i,j)),aimag(gamma5%mat(i,j)),1.0E-12_dp)
       end do
    end do
  end subroutine test_gamma5

  subroutine test_sum_4x4matrices()
    type(matrix4x4) :: A, B, C
    integer :: i, j
    call create_kronecker_delta()

    A%mat = delta
    B%mat = delta
    C = A + B

    do i = 1, 4
       do j = 1, 4
          call assert_close("A+B",real(C%mat(i,j)),2.0_dp*delta(i,j),1.0E-12_dp)
       end do
    end do
        
  end subroutine test_sum_4x4matrices

  subroutine test_substraction_4x4matrices()
    type(matrix4x4) :: A, B, C, D
    integer :: i, j
    call create_kronecker_delta()

    A%mat = delta
    B%mat = delta
    C = A - B
    D%mat = (0.0_dp,0.0_dp)
    call assert_equal_complex_matrix("A-B",C%mat,D%mat,1.0E-12_dp)

    D%mat = -delta
    B = -B
    call assert_equal_complex_matrix("-B",B%mat,D%mat,1.0E-12_dp)
        
  end subroutine test_substraction_4x4matrices

  
  subroutine test_pbc()
    integer :: L(4),p(4),m(4), i
    
    L = [16,15,14,13]
    call init_arrays(L)
    
    do i = 1, 4
       p = ip(L,i)
       m = im([1,1,1,1],i)
       call assert_equal("ip(L)",p(i), 1)
       call assert_equal("im(1)",m(i), L(i))
    end do
    
  end subroutine test_pbc

  subroutine test_cold_start()
    integer, parameter :: n = 4
    type(matrix3x3) :: U(n,n,n,n,n)
    complex(dp) :: D(3,3)
    integer :: i,j,k,l,m

    call create_kronecker_delta()
    call cold_start(U,[n,n,n,n,n])

    do i = 1, n
       do j = 1, n
          do k = 1, n
             do l = 1,  n
                do m = 1, n
                   D = U(i,j,k,l,m)%mat
                   call assert_equal_complex_matrix("cold start",D,delta_3x3,1.0E-12_dp)
                end do
             end do
          end do
       end do
    end do
    
  end subroutine test_cold_start

  subroutine test_su3_links()
    type(matrix3x3) :: A
    type(su3) :: U, V, W
    complex(dp) :: D(3,3), t

    call create_kronecker_delta
    D = delta_3x3
    call U%init()
    call V%init
    W = U*V
    V = U%dagger()
    t = tr(U)
    t = U%tr()
    call assert_equal_complex_matrix("su3 link",U%mat,D,1.0E-12_dp)
    call assert_equal_complex("Tr 1 = 3", t,(3.0_dp,0.0_dp),1.0e-12_dp)
    call assert_equal_complex_matrix("product su3",W%mat,D,1.0e-12_dp)
  end subroutine test_su3_links

  subroutine test_unitarity()
    type(su3) :: U,V
    real(dp) :: r(8)

    call random_number(r)
    call create_gellmann_matrices()
    call create_kronecker_delta()
    call U%init_su3(r(1),r(2),r(3),r(4),r(5),r(6),r(7),r(8))
    
    V = U*U%dagger()

    call assert_equal_complex_matrix("unitarity su3",V%mat,delta_3x3,1.0e-12_dp)
    
    
  end subroutine test_unitarity

  subroutine test_gellmann_matrices
    integer :: i, j
    complex(dp) :: r
    
    call create_gellmann_matrices()
    do i = 1, 8
       do j = 1, 8
          if( i == j )then
             r = (2.0_dp,0.0_dp)
          else
             r = (0.0_dp,0.0_dp)
          end if
          call assert_equal_complex("tr(l_i l_j) = 2d_ij",tr(gellmann_matrix(i)*gellmann_matrix(j)),r,1.0e-12_dp)
       end do
    end do
    
  end subroutine test_gellmann_matrices

  
end program test
