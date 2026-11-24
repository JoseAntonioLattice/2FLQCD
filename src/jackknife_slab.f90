!=====================================================================
! Errores jackknife para chi_t (metodo slab y medida directa)
!
! Uso en WF_analysis.f90:
!   use jackknife_slab
!   ...
!   ! prof(i,t) = slab_clov_profile(j,i,t) para TODAS las configuraciones
!   call chi_slab_jk(slab_clov_profile(j,:,:), tmin, tmax, nbin, A, dA)
!   call chi_direct_jk(top_char_clov(j,:), nbin, B, dB)
!   chi = A/V,  error = dA/V      (V = product(L))
!
! nbin = tamano del bloque (en configuraciones). Usa nbin = 1 y luego
! 2, 4, 8, 16: si el error crece con nbin hay autocorrelacion, y debes
! quedarte con el valor donde el error se estabiliza.
!=====================================================================
module jackknife_slab
  implicit none
  integer, parameter, private :: dp = 8
contains

  ! Ajuste de un parametro D(x) = A x(1-x) con pesos FIJOS w(t).
  ! Es lineal en D, asi que se puede reevaluar en cada muestra jackknife.
  pure function fit_A(Dbar, w, tmin, tmax, Lt) result(A)
    real(dp), intent(in) :: Dbar(:), w(:)
    integer,  intent(in) :: tmin, tmax, Lt
    real(dp) :: A, x, f, num, den
    integer  :: t
    num = 0.0_dp; den = 0.0_dp
    do t = tmin, tmax
       x = dble(t)/dble(Lt)
       f = x*(1.0_dp - x)
       num = num + w(t)*f*Dbar(t)
       den = den + w(t)*f*f
    end do
    A = num/den
  end function fit_A

  ! prof(nconf, Lt): perfil del slab por configuracion a un paso de flujo fijo
  ! Regresa A = chi_t V y su error jackknife dA.
  subroutine chi_slab_jk(prof, tmin, tmax, nbin, A, dA)
    real(dp), intent(in)  :: prof(:,:)
    integer,  intent(in)  :: tmin, tmax, nbin
    real(dp), intent(out) :: A, dA
    real(dp), allocatable :: D(:,:), Dsum(:), Djk(:), Ak(:), w(:), var(:)
    integer :: nconf, Lt, nb, n, k, t, i1, i2
    real(dp) :: x, Am

    nconf = size(prof,1); Lt = size(prof,2)
    nb = nconf/nbin            ! numero de bloques
    n  = nb*nbin               ! se descartan las ultimas configs sobrantes
    allocate(D(n,Lt), Dsum(Lt), Djk(Lt), Ak(nb), w(Lt), var(Lt))

    ! D_i(t) = prof_i(t) - x^2 prof_i(Lt)   (resta el sector de cada config)
    do t = 1, Lt
       x = dble(t)/dble(Lt)
       D(:,t) = prof(1:n,t) - x*x*prof(1:n,Lt)
    end do

    ! pesos fijos = 1/varianza de la media (muestra completa)
    Dsum = sum(D, dim=1)
    do t = 1, Lt
       var(t) = sum((D(:,t) - Dsum(t)/n)**2)/(dble(n-1)*dble(n))
       w(t) = 0.0_dp
       if (var(t) > 0.0_dp) w(t) = 1.0_dp/var(t)
    end do

    A = fit_A(Dsum/n, w, tmin, tmax, Lt)

    do k = 1, nb
       i1 = (k-1)*nbin + 1; i2 = k*nbin
       Djk = (Dsum - sum(D(i1:i2,:), dim=1))/dble(n - nbin)
       Ak(k) = fit_A(Djk, w, tmin, tmax, Lt)
    end do
    Am = sum(Ak)/nb
    dA = sqrt(dble(nb-1)/dble(nb)*sum((Ak - Am)**2))
  end subroutine chi_slab_jk

  ! Medida directa: B = <Q^2> - <Q>^2 (= chi_t V), con error jackknife
  subroutine chi_direct_jk(Q, nbin, B, dB)
    real(dp), intent(in)  :: Q(:)
    integer,  intent(in)  :: nbin
    real(dp), intent(out) :: B, dB
    real(dp), allocatable :: Bk(:)
    integer  :: nconf, nb, n, k, i1, i2
    real(dp) :: s1, s2, m1, m2, Bm

    nconf = size(Q); nb = nconf/nbin; n = nb*nbin
    allocate(Bk(nb))
    s1 = sum(Q(1:n)); s2 = sum(Q(1:n)**2)
    B = s2/n - (s1/n)**2
    do k = 1, nb
       i1 = (k-1)*nbin + 1; i2 = k*nbin
       m1 = (s1 - sum(Q(i1:i2)))/dble(n - nbin)
       m2 = (s2 - sum(Q(i1:i2)**2))/dble(n - nbin)
       Bk(k) = m2 - m1**2
    end do
    Bm = sum(Bk)/nb
    dB = sqrt(dble(nb-1)/dble(nb)*sum((Bk - Bm)**2))
  end subroutine chi_direct_jk

end module jackknife_slab
