module save
  use iso_fortran_env, only : dp => real64, i4 => int32
  use number2string
  use check_files_directories 
  use parameters, only : L, Lt, Lx, Ly, Lz, m0
  use su3facts
  implicit none
contains

  subroutine save_configuration(U,beta)
    type(su3), intent(in) :: U(4,Lt,Lx,Ly,Lz)
    real(dp), intent(in) :: beta
    character(:), allocatable :: path,file_name
    character(100), dimension(7) :: directory_array
    integer :: un
        
    directory_array = [character(100):: "data","configurations", &
         "Lt="//trim(int2str(Lt)),"Lx="//trim(int2str(Lx)),"Ly="//trim(int2str(Ly)),"Lz="//trim(int2str(Lz)), &
         "beta="//trim(real2str(beta,1,4))]
       
    call check_directory(directory_array,path)
    call numbered_files(path,"U",".bin",file_name)
    open(newunit = un, file = file_name, access = "sequential", form = "unformatted")
    
    write(un) U
    close(un)
        
  end subroutine save_configuration

  subroutine read_configuration(U,filename)
    type(su3), intent(out) :: U(4,Lt,Lx,Ly,Lz)
    character(*), intent(in) :: filename
    integer :: un
        
    open(newunit = un, file = filename, access = "sequential", form = "unformatted", action = "read")  
    read(un) U
    close(un)
        
  end subroutine read_configuration



end module save

