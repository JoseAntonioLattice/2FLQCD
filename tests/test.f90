program test
  use random
  use starts
  use pbc
  use su3facts, gamma => dirac_matrix
  use gauge
  use parameters
  implicit none
  integer, parameter :: dp = 8

  
  call test_mean_zero()
  call test_variance_one()
  call test_gamma()
  call test_gamma5()
  !call test_pbc()

  call test_sum_4x4matrices
  call test_substraction_4x4matrices()

  call test_cold_start()

  call test_su3_links

  call test_gellmann_matrices
  call test_unitarity()

  call test_assert_matrix()
  call test_plaquette()

  call test_TA
  call test_su3alg

  call test_parameters
  
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
    logical, dimension(size(measured(:,1)),size(measured(1,:))) :: condition
    integer :: n, m, i, j

    n = size(measured(:,1))
    m = size(measured(1,:))
    condition = abs(real(measured - expected)) < tol .and. abs(aimag(measured - expected)) < tol
    if( all(condition) ) then
       print"(2a)", " [PASS] : ", name
    else
       print"(2a)", " [FAIL] : ", name
       do i = 1, n
          do j = 1, m
             if(condition(i,j) .eqv. .false.) &
                  print"(4x,a,i0,a,i0,a,2ES12.4,a,2ES12.4,a,ES12.4)", &
                  "(",i,",",j,") Expected: ", expected(i,j), " Measured: ", measured(i,j), " tol = ", tol
             
          end do
       end do
    end if

    
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
    !B(2,2) = 2.0_dp
    !B(1,2) = 0.0_dp
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
    complex(dp) :: delta(4,4)
    character(99) :: str
    
    call create_kronecker_delta(delta,4)
    call create_gamma_matrices()
    
    do mu = 1, 4
       do nu = 1, 4
          product = gamma(mu)*gamma(nu) + gamma(nu)*gamma(mu)          
          write(str,"(a,i0,a,i0,a,2i0)") "{gamma_",mu,",gamma_",nu,"} = delta_",mu,nu
          call assert_equal_complex_matrix(trim(str), product%mat,2*delta(mu,nu)*delta, 1.0E-12_dp) 
       end do
    end do
    
  end subroutine test_gamma

  subroutine test_gamma5()
    integer :: i, j
    type(matrix4x4) :: gamma5

    
    call create_gamma_matrices()
    gamma5 = gamma(1)*gamma(2)*gamma(3)*gamma(4) 
    call assert_equal_complex_matrix("gamma5", gamma(5)%mat,gamma5%mat,1.0E-12_dp)
         
  end subroutine test_gamma5

  subroutine test_sum_4x4matrices()
    type(matrix4x4) :: A, B, C
    integer :: i, j
    complex(dp) :: delta(4,4)
    
    call create_kronecker_delta(delta,4)
   
    A%mat = delta
    B%mat = delta
    C = A + B

    call assert_equal_complex_matrix("A+B",C%mat,2*delta,1.0e-12_dp)
  end subroutine test_sum_4x4matrices

  subroutine test_substraction_4x4matrices()
    type(matrix4x4) :: A, B, C, D
    integer :: i, j
    complex(dp) :: delta(4,4)
    
    call create_kronecker_delta(delta,4)

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
    !deallocate(ip1)!,ip2,ip3,ip4,im1,im2,im3,im4)
  end subroutine test_pbc

  subroutine test_cold_start()
    integer, parameter :: n = 4
    type(matrix3x3) :: U(n,n,n,n,n)
    complex(dp) :: D(3,3)
    integer :: i,j,k,l,m
    logical :: condition(n,n,n,n,n)
    real(dp) :: tol = 1.0e-12_dp
    character(99) :: str

   
    call cold_start(U,[n,n,n,n,n])
    
    print*, "Cold start"
    do i = 1, n
       do j = 1, n
          do k = 1, n
             do l = 1,  n
                do m = 1, n
                   condition(i,j,k,l,m) = all(abs( real(U(i,j,k,l,m)%mat - delta_3x3)) < tol &
                   .and. abs(aimag(U(i,j,k,l,m)%mat - delta_3x3)) < tol)
                   if(.not.condition(i,j,k,l,m)) then
                      D = U(i,j,k,l,m)%mat
                      write(str,"(5(a,i0),a)") "x = (",i,",",j,",",k,",",l,",", m,")"
                      call assert_equal_complex_matrix(trim(str),D,delta_3x3,tol)
                   end if
                end do
             end do
          end do
       end do
    end do
    if(all(condition)) then
       print*, "  [PASS]  : cold start"
    end if
  end subroutine test_cold_start

  subroutine test_su3_links()
    type(matrix3x3) :: A
    type(su3) :: U, V, W
    complex(dp) :: D(3,3), t

    
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
    r(1:5) = 0.0_dp
    call create_gellmann_matrices()
    
    call U%init_su3(r(1),r(2),r(3),r(4),r(5),r(6),r(7),r(8))
    
    V = U*U%dagger()

    print*, "Test su(3) unitarity"
    call assert_equal_complex_matrix("UU^dagger = 1",V%mat,delta_3x3,1.0e-12_dp)
    call assert_equal_complex("det(U) = 1",determinant(U%mat),(1.0_dp,0.0_dp),1.0e-12_dp)
    
  end subroutine test_unitarity

  subroutine test_gellmann_matrices
    integer :: i, j
    complex(dp) :: delta(8,8), trace(8,8)
    logical :: condition(8,8)


    print*, "Test Gell-Mann Matrices"
    call create_gellmann_matrices()
    call create_kronecker_delta(delta,8)
    do i = 1, 8
       do j = 1, 8
         trace(i,j) = tr(gellmann_matrix(i)*gellmann_matrix(j))
       end do
    end do

    call assert_equal_complex_matrix("tr(l_i l_j) = 2d_ij",trace,2*delta,1.0e-12_dp)
    
  end subroutine test_gellmann_matrices

  subroutine test_plaquette()
    integer, parameter :: n = 2
    type(su3), dimension(4,n,n,n,n) :: U, V
    type(su3), dimension(n,n,n,n) :: omega
    real(dp), dimension(n,n,n,n) :: r1,r2,r3,r4,r5,r6,r7,r8
    integer :: i,j,k,l,mu,nu
    logical :: condition(4,4,n,n,n,n)
    type(su3) :: plq
    real(dp) :: tol = 1.0e-12_dp
    character(99) :: str
    call U%init()

    call random_number(r1)
    call random_number(r2)
    call random_number(r3)
    call random_number(r4)
    call random_number(r5)
    call random_number(r6)
    call random_number(r7)
    call random_number(r8)

    call omega%init_su3(r1,r2,r3,r4,r5,r6,r7,r8)

    call init_arrays([n,n,n,n])
    call gauge_transformation(U,omega)

    condition = .true.
    print*, "Plaquette invariance"

    do i = 1, n
       do j = 1, n
          do k = 1, n
             do l = 1,  n
                do mu = 1, 3
                   do nu = mu + 1, 4
                      plq = plaquette(U,[i,j,k,l],mu,nu)
                      condition(mu,nu,i,j,k,l) = all(abs( real(plq%mat - delta_3x3)) < tol) &
                           .and. all(abs(aimag(plq%mat - delta_3x3)) < tol)
                      if(.not.condition(mu,nu,i,j,k,l)) then
                         
                         write(str,"(4(a,i0),a,2i0,a)") " [FAIL] plq inv at x = (",i,",",j,",",k,",",l,")(",mu,nu,")"
                         call assert_equal_complex_matrix(trim(str),plq%mat,delta_3x3,tol)
                      end if
                   end do
                end do
             end do
          end do
       end do
    end do

    
    if(all(condition)) then
       print*, "  [PASS]  : plaquette invariance"
    end if
  end subroutine test_plaquette


  subroutine gauge_transformation(U,omega)
    type(su3), intent(inout) :: U(:,:,:,:,:)
    type(su3), intent(in) :: omega(:,:,:,:)
    integer :: mu, i,j,k,l,m, xp(4)
    
    do i = 1, size(U(1,:,1,1,1))
       do j = 1, size(U(1,1,:,1,1))
          do k = 1, size(U(1,1,1,:,1))
             do l = 1, size(U(1,1,1,1,:))
                do mu = 1, size(U(:,1,1,1,1))
                   xp = ip([i,j,k,l,m], mu)
                   U(mu,i,j,k,l) = omega(i,j,k,l)*U(mu,i,j,k,l)*dagger(omega(xp(1),xp(2),xp(3),xp(4)))
                end do
             end do
          end do
       end do
    end do
    
  end subroutine gauge_transformation

  subroutine test_TA()
    type(matrix3x3) :: W1, TA1, TA2, TA3, DTA2, DTA3
    type(su3) :: W2 
    type(su3alg) :: W3
    complex(dp) :: t1, t2, t3
    real(dp) :: r1,r2,r3,r4,r5,r6,r7,r8

    call random_number(r1)
    call random_number(r2)
    call random_number(r3)
    call random_number(r4)
    call random_number(r5)
    call random_number(r6)
    call random_number(r7)
    call random_number(r8)
    r4 = 0.0_dp
    r5 = 0.0_dp
    r6 = 0.0_dp
    r7 = 0.0_dp
    r8 = 0.0_dp
    call W2%init_su3(r1,r2,r3,r4,r5,r6,r7,r8)

    call random_number(r1)
    call random_number(r2)
    call random_number(r3)
    call random_number(r4)
    call random_number(r5)
    call random_number(r6)
    call random_number(r7)
    call random_number(r8)
    
    call W3%init_su3alg(r1,r2,r3,r4,r5,r6,r7,r8)

    TA2 = TA(W2)
    TA3 = W3%TA()

    DTA2 = dagger(TA2)
    DTA3 = dagger(TA3)
    
    t2 = tr(TA2)
    t3 = tr(TA3)
       
    call assert_equal_complex("tr W2_TA = 0", t2,(0.0_dp,0.0_dp),1.0e-12_dp)
    call assert_equal_complex("tr W3_TA = 0", t3,(0.0_dp,0.0_dp),1.0e-12_dp)

    call assert_equal_complex_matrix("W2^dagger = -W2", DTA2%mat,-TA2%mat,1.0e-12_dp)
    call assert_equal_complex_matrix("W3^dagger = -W3", DTA3%mat,-TA3%mat,1.0e-12_dp)
   
  end subroutine test_TA

  subroutine test_su3alg()
    real(dp) :: r(8)
    type(su3alg) :: P, Pd


    call random_number(r)
    
    call P%init_su3alg(r(1),r(2),r(3),r(4),r(5),r(6),r(7),r(8))
    pd = p%dagger()
    print*, "Test su(3) algebra element"
    call assert_equal_complex("Tr u = 0", tr(p),(0.0_dp,0.0_dp),1.0e-12_dp)
    call assert_equal_complex_matrix("u^dagger = u", pd%mat,p%mat,1.0e-12_dp)

  end subroutine test_su3alg

  subroutine test_parameters

    call read_input()
    
  end subroutine test_parameters
  
  
end program test
