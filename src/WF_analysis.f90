program main

  use number2string
  use statistics
  implicit none
  character(:), allocatable :: filename
  integer, dimension(4) :: L
  integer :: i, j, N, Ntime, inunit, idx, outunit, t, unitslab, ios
  real(8) :: beta
  real(8), dimension(:,:), allocatable :: E_plq, E_clov, top_char, newtopchar
  real(8), dimension(:,:,:), allocatable :: slab_plq,slab_clov, slab_plq_profile,slab_clov_profile
  real(8), dimension(:), allocatable :: time, alpha
  character(512), allocatable :: configs(:), slabs(:)
  
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

  allocate(E_plq(0:Ntime,N))
  allocate(E_clov, top_char,newtopchar, mold = E_plq)
  allocate(configs(N), slabs(N))
  allocate(time(0:Ntime))
  allocate(alpha, mold = time)
  
  allocate(slab_plq(0:ntime,n,L(1)))
  allocate(slab_clov, slab_plq_profile, slab_clov_profile, mold = slab_plq)
  
  do i = 1, N
     read(*,'(a)') configs(i)
  end do
  
  do i = 1, N
     read(*,'(a)') slabs(i)
  end do

  
  do i = 1, N
     print"(*(a))", trim(configs(i))
     open(newunit = inunit, file = trim(configs(i)), status = "old")
     
     do j = 0, Ntime
        read(inunit,*) time(j), E_plq(j,i), E_clov(j,i), top_char(j,i)
     end do
    
  end do

  open(unit=666,file = "topchar.dat", status = "unknown")
  do j = 0, ntime
     call rescale_charge(n,top_char(j,:),alpha(j),newtopchar(j,:))
  end do

  do i =1, n
     do j = 0, ntime
        write(666,*) newtopchar(j,i)
     end do
     write(666,"(2/)")
  end do
  
  print*, newtopchar(0,:)
  do i = 1, N
     open(newunit = unitslab, file = trim(slabs(i)), status = "old")
     print"(a)", trim(slabs(i))
     do j = 0, Ntime
        do t = 1, L(1)
           read(unitslab,*,iostat=ios) idx, slab_clov(j,i,t), slab_clov_profile(j,i,t), slab_plq(j,i,t), slab_plq_profile(j,i,t) 
        end do

        if (ios /= 0) then
           print *, "Error leyendo ", trim(slabs(i)), " en j=", j, " t=", t, " ios=", ios
           stop
        end if
        
     end do
     close(unitslab)
  end do

  open(newunit = outunit, file = "WF.dat" ,status = "unknown")
  do j = 0, ntime
     write(outunit,*) time(j), avr(E_plq(j,:)), stderr(E_plq(j,:)), &
          avr(E_clov(j,:)), stderr(E_clov(j,:)), &
          avr(top_char(j,:)), stderr(top_char(j,:))    
  end do
  
  open(newunit = unitslab, file = "WF_slab.dat" ,status = "unknown")

  do j = 0, ntime
     do t = 1, L(1)
        write(unitslab,*) t, avr(slab_clov(j,:,t))/alpha(j), stderr(slab_clov(j,:,t))/alpha(j), &
             avr(slab_plq(j,:,t))/alpha(j), stderr(slab_plq(j,:,t))/alpha(j),&
             avr(slab_clov_profile(j,:,t))/alpha(j), stderr(slab_clov_profile(j,:,t))/alpha(j)
     end do
     write(unitslab,"(2/)")
  end do


  
  open(newunit = unitslab, file = "WF_slab_new.dat" ,status = "unknown")
  do j = 0, ntime
     do t = 1, L(1)
        write(unitslab,*) t, slab_clov(j,1,t)*alpha(j), &
             slab_plq(j,1,t)*alpha(j), &
        slab_plq_profile(j,1,t)*alpha(j), slab_clov_profile(j,1,t)*alpha(j)
     end do
     write(unitslab,"(2/)")
  end do
contains

  subroutine rescale_charge(nconf, qL, alpha_best, qint)
    implicit none

    integer, intent(in)  :: nconf
    real(8), intent(in)  :: qL(nconf)
    real(8), intent(out) :: alpha_best
    real(8), intent(out) :: qint(nconf)
    real(8) :: alpha, d, f, fbest
    integer :: i, k
    integer, parameter :: nscan = 50000
    real(8), parameter :: amin = 1.0d0, amax = 1.5d0
    
    fbest = huge(1.0d0)
    do k = 0, nscan
       alpha = amin + (amax - amin) * dble(k) / dble(nscan)
       f = 0.0d0
       do i = 1, nconf
          d = alpha*qL(i) - dnint(alpha*qL(i))
          f = f + d*d
       end do
       f = f / dble(nconf)
       if (f < fbest) then
          fbest = f
          alpha_best = alpha
       end if
    end do
    
    do i = 1, nconf
       qint(i) = 1.0d0*nint(alpha_best * qL(i))
    end do
  end subroutine rescale_charge
  
end program main
