module statistics
    
  implicit none
  integer, parameter, private :: dp = 8

contains
  
  function avr(x)
    real(dp) :: avr, x(:)
    avr = sum(x)/size(x)
  end function avr

  function var(x)
    real(dp) :: var, x(:), xbar
    xbar = avr(x)
    var = 1.0_dp/(size(x) - 1) * sum( (x - xbar)**2 )
  end function var

  function stderr(x)
    real(dp) :: stderr, x(:)
    stderr = sqrt(var(x)/size(x))
  end function stderr

  function jackknife(x,bins)
    real(dp) :: jackknife, x(:)
    integer, intent(in) :: bins
    integer :: MM, NN, i
    real(dp) :: xbar, sum_x
    real(dp) :: x_m(bins)

    NN = size(x)
    MM = NN/bins

    xbar = avr(x)
    sum_x = sum(x)
    x_m = 1.0_dp/(NN-MM) * [(sum_x - sum(x(MM*(i-1)+1:MM*i)),i=1,bins)]

    jackknife = sqrt( real(bins - 1,dp)/bins * sum( (x_m - xbar)**2) )
  end function jackknife

  function jackknife2(x)
    real(dp) :: jackknife2, x(:)
    integer, dimension(4), parameter :: bins = [4,10,20,25]
    real(dp), dimension(size(bins)) :: arr
    integer :: i

    arr = 0.0_dp
    do i = 1, size(bins)
       if(mod(size(x),bins(i)) /= 0) cycle
       arr = jackknife(x,bins(i))
    end do

    jackknife2 = maxval(arr)
    
  end function jackknife2
  
end module statistics
