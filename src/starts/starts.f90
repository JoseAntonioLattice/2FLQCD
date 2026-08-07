module starts
  
  use su3facts
  implicit none

contains

  subroutine cold_start(U,L)
    integer, dimension(5) :: L
    type(matrix3x3), dimension(L(1),L(2),L(3),L(4),L(5)) :: U
    integer :: i,j,k,n, m

    call u%init()
    
  end subroutine cold_start
  
end module starts
