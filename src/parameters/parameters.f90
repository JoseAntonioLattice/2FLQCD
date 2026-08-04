module parameters

  integer, dimension(4) :: L

  namelist /lattice/ L 

contains

  subroutine read_input()
    character(100) :: filename
    read(*,'(a)') filename
    write(*,*) "User typed: ", trim(filename)
  end subroutine read_input
  
end module parameters
