program main
  use parameters
  use dynamics
  use arrays
  implicit none

  write(*,*) "2 flavor QCD"
  call read_input()
  call initialize()

  call simulation(U,beta)

end program main
