module random
  implicit none
  integer, parameter, private :: dp = 8
  real(dp), parameter, private :: pi = acos(-1.0_dp)
  real(dp), parameter, private :: twopi = 2*pi

contains
  
  elemental subroutine elemental_rgauss(r,u1,u2)
    real(dp), intent(in) :: u1, u2
    real(dp), intent(out) :: r
    r = sqrt(-2*log(1.0_dp - u1)) * cos(twopi*u2)
  end subroutine elemental_rgauss
  
end module random
