module parameters
  
  implicit  none
  integer, parameter, private :: dp = 8
  
  integer, dimension(4) :: L
  integer :: Lt,Lx,Ly,Lz
  integer :: N
  real(8) :: epsilon
  integer :: N_measurements, N_thermalization, N_skip
  real(dp):: betai, betaf
  integer :: nbeta
  character(50):: algorithm
  namelist /lattice/ L, N, epsilon, N_measurements, N_thermalization, N_skip, &
       betai, betaf, nbeta, algorithm

contains

  subroutine read_input()
    character(100) :: filename
    integer :: inunit
    
    read(*,"(a)") filename
    print*, 'User typed: ', trim(filename)
    open(newunit = inunit,file = trim(filename), status = 'old')
    read(inunit, nml = lattice)
    write(*,nml = lattice)
    Lt = L(1)
    Lx = L(2)
    Ly = L(3)
    Lz = L(4)
    
  end subroutine read_input
  
end module parameters
