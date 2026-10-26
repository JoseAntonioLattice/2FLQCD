program main

  use number2string
  use statistics
  implicit none
  character(:), allocatable :: filename
  integer, dimension(4) :: L
  integer :: i, j, N, Ntime, inunit, idx
  real(8) :: beta
  real(8), dimension(:,:), allocatable :: E_plq, E_clov, top_char
  real(8), dimension(:), allocatable :: time
  character(512), allocatable :: configs(:)

  
  print*, "Enter lattice array"
  read(*,*) L
  print*, "User typed: L = ", L
  print*, "Enter beta value"
  read(*,*) beta
  print*, "User typed: beta = ", beta
  print*, "Enter number of configurations"
  read(*,*) N
  print*, "User typed: N =", N
  print*, "Enter number of time steps"
  read(*,*) Ntime
  print*, "User typed: Ntime =", Ntime

  allocate(configs(N))
  do i = 1, N
     read(*,'(a)') configs(i)
  end do
  
  do i = 1, N
     print"(a)", trim(configs(i))
  end do

  allocate(time(1:Ntime))
  allocate(E_plq(N,0:Ntime))
  allocate(E_clov,top_char, mold = E_plq)


    
end program main
