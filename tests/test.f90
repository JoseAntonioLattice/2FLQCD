program test
  use random
  use su3facts
  use pbc
  implicit none
  integer, parameter :: dp = 8

  
  !call test_mean_zero()
  !call test_variance_one()
  !call test_gamma()
  !call test_gamma5()
  !call test_pbc()
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

  
  subroutine assert_equal(name, measured, expected)
    character(*), intent(in) :: name
    integer, intent(in) :: measured, expected

    if( measured == expected ) then
       print"(2a)", " [PASS] : ", name
    else
       print"(3a,i0,a,i0)", " [FAIL] : ", name, " Expected:", expected, " Measured: ", measured
    end if
  end subroutine assert_equal
  
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
    integer :: mu, nu, i, j
    complex(dp), dimension(4,4) :: product
    logical :: condition
    
    call create_kronecker_delta()
    call create_gamma_matrices()
    
    do mu = 1, 4
       do nu = 1, 4
          product = matmul(gamma(mu,:,:),gamma(nu,:,:)) + matmul(gamma(nu,:,:),gamma(mu,:,:))
          print*, "mu = ",mu, "nu = ", nu
          do i = 1, 4
             do j = 1, 4
                call assert_close("{gamma_mu,gamma_nu}", real(product(i,j)),2*delta(mu,nu)*delta(i,j), 1.0E-12_dp) 
             end do
          end do
       end do
    end do
    
  end subroutine test_gamma

  subroutine test_gamma5()
    integer :: i, j
    complex(dp) :: gamma5(4,4)

    
    call create_kronecker_delta()
    call create_gamma_matrices()
    gamma5 = matmul(gamma(1,:,:),gamma(2,:,:))
    gamma5 = matmul(gamma5,gamma(3,:,:))
    gamma5 = matmul(gamma5,gamma(4,:,:))
    do i = 1, 4
       do j = 1, 4
          call assert_close("gamma5 real",real(gamma(5,i,j)),real(gamma5(i,j)),1.0E-12_dp)
          call assert_close("gamma5 imag",aimag(gamma(5,i,j)),aimag(gamma5(i,j)),1.0E-12_dp)
       end do
    end do
  end subroutine test_gamma5

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
  
end program test
