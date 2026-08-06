module starts
  
  use su3facts
  implicit none

contains

  subroutine cold_start(U,L)
    integer, dimension(5) :: L
    type(matrix3x3), dimension(L(1),L(2),L(3),L(4),L(5)) :: U
    integer :: i,j,k,n, m

    do i = 1, L(1)
       do j = 1, L(2)
          do k = 1, L(3)
             do m = 1, L(4)
                do n = 1, L(5)
                   U(i,j,k,m,n)%mat = delta_3x3
                end do
             end do
          end do
       end do
    end do
    
  end subroutine cold_start
  
end module starts
