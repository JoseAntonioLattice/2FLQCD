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

  allocate(time(0:Ntime))
  allocate(E_plq(N,0:Ntime))
  allocate(E_clov,top_char, mold = E_plq)
   
  do i = 1, N
     filename = "data/clusterdata/WF_Lt="//int2str(L(1)) &
          //"_Lx="//int2str(L(2)) &
          //"_Ly="//int2str(L(3)) &
          //"_Lz="//int2str(L(4)) &
          //"_beta="//real2str(beta,1,4) &
          //"_"//int2str(i)//".dat"
     print*, filename
     !call read_file(filename,[idx, time, E_plq, E_clov])
     open(newunit = inunit, file = filename, status = "old")

    
     
     !j = 0
     do j = 0, Ntime
        !read(inunit,*,iostat = stat) d_idx, d_t, d_Eplq, d_Eclov
        !if(stat /= 0) exit
        !j = j + 1
        !idx = [idx, d_idx]
        !E_plq(i,j) = [E_plq(i,j), d_Eplq]
        !E_clov = [E_clov, d_Eclov]

        read(inunit,*) idx, time(j), E_plq(i,j), E_clov(i,j), top_char(i,j)
        
     end do
     close(inunit)
  end do

  open(unit = 69, file = "data/WF.dat", status = "unknown")
  do j = 0, Ntime
     write(69,*) time(j), &
          avr(E_plq(:,j)), stderr(E_plq(:,j)), &
          avr(E_clov(:,j)), stderr(E_clov(:,j)), &
          avr(top_char(:,j)), stderr(top_char(:,j)) 
  end do

  
  
end program main
