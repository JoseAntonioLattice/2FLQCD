module pbc

  implicit none
  integer, dimension(:), allocatable, private :: ip1, ip2, ip3, ip4, im1, im2, im3, im4
contains


  subroutine set_pbc(L)
    integer, intent(in) :: L(4)
    integer :: i

    allocate(ip1(L(1)),ip2(L(2)),ip3(L(3)),ip4(L(4)))
    allocate(im1(L(1)),im2(L(2)),im3(L(3)),im4(L(4)))

    call set_arrays(ip1,im1,L(1))
    call set_arrays(ip2,im2,L(2))
    call set_arrays(ip3,im3,L(3))
    call set_arrays(ip4,im4,L(4))
    
  end subroutine set_pbc

  subroutine set_arrays(p,m,L)
    integer, intent(in) :: L
    integer, intent(out), dimension(L) :: p, m
    integer :: i

    do i = 1, L
       p(i) = i + 1
       m(i) = i - 1
    end do
    p(L) = 1
    m(1) = L
    
  end subroutine set_arrays
  
  
  function ip(x,mu)
    integer :: ip(4)
    integer, intent(in) :: x(4), mu

    ip = x
    select case(mu)
    case(1)
       ip(1) = ip1(x(1))
    case(2)
       ip(2) = ip2(x(2))
    case(3)
       ip(3) = ip3(x(3))
    case(4)
       ip(4) = ip4(x(4))
    end select
  end function ip

  function im(x,mu)
    integer :: im(4)
    integer, intent(in) :: x(4), mu
    
    im = x
    select case(mu)
    case(1)
       im(1) = im1(x(1))
    case(2)
       im(2) = im2(x(2))
    case(3)
       im(3) = im3(x(3))
    case(4)
       im(4) = im4(x(4))
    end select
  end function im

  
end module pbc
