module format

  use number2string
  use iso_fortran_env, only : dp => real64
  implicit none
  
contains

  function exponent(x)
    integer :: exponent
    real(dp) :: x
    exponent = floor(log10(x))
  end function exponent

  function error(x,n)
    integer :: error
    real(dp) :: x
    integer, intent(in), optional :: n
    integer :: l
    
    if(present(n))then
       l = n
    else
       l = 2
    end if
        
    error = nint(x*10**(-exponent(x)+l-1))
    
  end function error

  function format_error(x,x_err,n)
    character(:), allocatable :: format_error
    real(dp) :: x, x_err
    character(100) :: k, format
    integer, intent(in), optional :: n
    integer :: l

    if(present(n))then
       l = n
    else
       l = 2
    end if
    
    format = '(f20.'//int2str(-exponent(x_err)+l-1)//',"(",i'//int2str(l)//',")" )' 
    write(k,format) x,error(x_err,l)
    format_error = trim(adjustl(k))
    
  end function format_error

  
end module format
