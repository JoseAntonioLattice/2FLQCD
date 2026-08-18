module starts
  use parameters, only : Lt, Lx, Ly, Lz
  use su3facts
  implicit none
  integer, parameter, private :: dp = 8
contains

  subroutine cold_start(U)
    type(su3), dimension(4,Lt,Lx,Ly,Lz) :: U

    call u%init()
    
  end subroutine cold_start

  subroutine hot_start(U)    
    type(su3), dimension(4,Lt,Lx,Ly,Lz) :: U
    real(dp), dimension(4,Lt,Lx,Ly,Lz) :: r1,r2,r3,r4,r5,r6,r7,r8

    call random_number(r1)
    call random_number(r2)
    call random_number(r3)
    call random_number(r4)
    call random_number(r5)
    call random_number(r6)
    call random_number(r7)
    call random_number(r8)
    
    call u%init_su3(r1,r2,r3,r4,r5,r6,r7,r8)

  end subroutine hot_start
  
end module starts
