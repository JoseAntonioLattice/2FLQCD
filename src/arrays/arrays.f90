module arrays
  use su3facts, only: su3, su3alg, matrix3x3
  implicit none
  integer, parameter, private :: dp = 8
  type(su3), allocatable :: U(:,:,:,:,:)
  real(dp), allocatable  :: beta(:)
  real(dp), allocatable  :: a_plqv(:) 
end module arrays
