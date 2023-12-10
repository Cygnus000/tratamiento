program logkillH
    use, intrinsic :: iso_fortran_env
    implicit none

    real(kind=real64), parameter :: x0 = 30.0d0
    real(kind=real64), parameter :: g = 0.1d0
    real(kind=real64), parameter :: d = 0.04d0

    real(kind=real64), parameter :: t0 = 0.0d0
    real(kind=real64), parameter :: tmax = 100.0d0
    integer          , parameter :: N    = 10000
    real(kind=real64), parameter :: dt = (tmax - t0) / dble(N)
    integer          , parameter :: N_equ = 1    ! Numero de ecuaciones

    integer           :: i
    real(kind=real64) :: r(N_equ)
    real(kind=real64) :: t(N), x(N), farmaco(N)
!**********************************************************************
    t = [ ( dt * i, i = 1, N ) ]             ! llenando vector temporal
!**********************************************************************
    r = [ x0 ]                                      ! valores iniciales
!**********************************************************************
    do i = 1, N                                           ! resolviendo
        farmaco(i) = 30.0d0*(sqrt(t(i))/(sqrt(10.0d0)+sqrt(t(i))))
        x(i) = r(1)

        r = r + rk4( r, t(i), farmaco(i), dt )
!**********************************************************************
        open(1,file='logkill-H.dat')                 ! llenando archivo
          write(1,*) t(i), x(i), farmaco(i)
          print*,    t(i), x(i), farmaco(i)
    end do
    close(1) 
    call system('gnuplot -c logkill-H.p')
!**********************************************************************
contains
!**********************************************************************
    pure function f(r, t, far) ! Aqui se colocan las ecuaciones a resol
        real(kind=real64), intent(in) :: r(N_equ) ! Valores
        real(kind=real64), intent(in) :: t    ! Paso
        real(kind=real64), intent(in) :: far
        real(kind=real64)             :: f(N_equ)
        real(kind=real64)             :: u

        u = r(1)

        f(1) = g*u-d*u*far
        
    end function f
!**********************************************************************
    pure function rk4(r, t, far, dt)                    ! Runge-Kutta 4
        real(kind=real64), intent(in) :: r(N_equ) ! Valores
        real(kind=real64), intent(in) :: t    ! Paso
        real(kind=real64), intent(in) :: far  !funcion cantidad farmaco
        real(kind=real64), intent(in) :: dt   ! Tamano de paso
        real(kind=real64)             :: rk4(N_equ)
        real(kind=real64)             :: k1(N_equ), k2(N_equ)
        real(kind=real64)             :: k3(N_equ), k4(N_equ)   

        k1 = dt * f( r,            t               ,  far )
        k2 = dt * f( r + 0.5d0 * k1, t + 0.5d0 * dt,  far )
        k3 = dt * f( r + 0.5d0 * k2, t + 0.5d0 * dt,  far )
        k4 = dt * f( r + k3        , t + dt        ,  far )

        rk4 = ( k1 + ( 2.0d0 * k2 ) + ( 2.0d0 * k3 ) + k4 ) / 6.0d0
    end function rk4
!**********************************************************************
end program logkillH