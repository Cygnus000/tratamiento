program stepanova
    use, intrinsic :: iso_fortran_env
    implicit none

    real(kind=real64), parameter :: x0 = 1.0d0
    real(kind=real64), parameter :: z0 = 0.0d0
    real(kind=real64), parameter :: g = 0.6d0
    real(kind=real64), parameter :: x_max = 7.5d0
    real(kind=real64), parameter :: gama = 1.0d0
    real(kind=real64), parameter :: kI = 0.5d0
    real(kind=real64), parameter :: beta = 0.3d0
    real(kind=real64), parameter :: delta = 0.4d0
    real(kind=real64), parameter :: mu = 0.1d0

    real(kind=real64), parameter :: t0 = 0.0d0
    real(kind=real64), parameter :: tmax = 100.0d0
    integer          , parameter :: N    = 10000
    real(kind=real64), parameter :: dt = (tmax - t0) / dble(N)
    integer          , parameter :: N_equ = 2    ! Numero de ecuaciones

    integer           :: i
    real(kind=real64) :: r(N_equ)
    real(kind=real64) :: t(N), x(N), z(N)
!**********************************************************************
    t = [ ( dt * i, i = 1, N ) ]             ! llenando vector temporal
!**********************************************************************
    r = [ x0, z0 ]                                  ! valores iniciales
!**********************************************************************
    do i = 1, N                                           ! resolviendo
        x(i) = r(1)
        z(i) = r(2)

        r = r + rk4( r, t(i), dt )
!**********************************************************************
        open(1,file='nova.dat')                      ! llenando archivo
          write(1,*) t(i), x(i), z(i)
          print*,    t(i), x(i), z(i)
    end do
    close(1) 
    call system('gnuplot -c nova.p')
!**********************************************************************
contains
!**********************************************************************
    pure function f(r, t) ! Aqui se colocan las ecuaciones a resol
        real(kind=real64), intent(in) :: r(N_equ) ! Valores
        real(kind=real64), intent(in) :: t    ! Paso
        real(kind=real64)             :: f(N_equ)
        real(kind=real64)             :: u,v

        u = r(1)
        v = r(2)

        f(1) = g * u * (1.0 - u/x_max) - gama * u * v
        f(2) = kI * ( u - beta * u**2 ) * v - delta * v + mu
        
    end function f
!**********************************************************************
    pure function rk4(r, t, dt)                         ! Runge-Kutta 4
        real(kind=real64), intent(in) :: r(N_equ) ! Valores
        real(kind=real64), intent(in) :: t    ! Paso
        real(kind=real64), intent(in) :: dt   ! Tamano de paso
        real(kind=real64)             :: rk4(N_equ)
        real(kind=real64)             :: k1(N_equ), k2(N_equ)
        real(kind=real64)             :: k3(N_equ), k4(N_equ)   

        k1 = dt * f( r,            t                )
        k2 = dt * f( r + 0.5d0 * k1, t + 0.5d0 * dt )
        k3 = dt * f( r + 0.5d0 * k2, t + 0.5d0 * dt )
        k4 = dt * f( r + k3        , t + dt         )

        rk4 = ( k1 + ( 2.0d0 * k2 ) + ( 2.0d0 * k3 ) + k4 ) / 6.0d0
    end function rk4
!**********************************************************************
end program stepanova
