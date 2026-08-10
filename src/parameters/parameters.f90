module parameters

  
  integer, dimension(4) :: L
  integer :: Lt,Lx,Ly,Lz
  integer :: N
  real(8) :: epsilon
  namelist /lattice/ L, N, epsilon 

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
