program logkillH
    use, intrinsic :: iso_fortran_env, only: qp=>real128
    implicit none

    real(qp), parameter :: x0 = 30.0_qp
    real(qp), parameter :: g = 0.1_qp
    real(qp), parameter :: d = 0.04_qp

    real(qp), parameter :: t0 = 0.0_qp
    real(qp), parameter :: tmax = 100.0_qp
    integer , parameter :: N    = 10000
    real(qp), parameter :: dt = (tmax - t0) / dble(N)
    integer , parameter :: N_equ = 1    ! Numero de ecuaciones

    integer  :: i
    real(qp) :: r(N_equ), t(N), x(N), farmaco(N)
!**********************************************************************
    t = [ ( dt * i, i = 1, N ) ]             ! llenando vector temporal
!**********************************************************************
    r = [ x0 ]                                      ! valores iniciales
!**********************************************************************
    open(1,file='logkill-H.dat')                 ! llenando archivo
    do i = 1, N                                           ! resolviendo
        farmaco(i) = 30.0_qp*(sqrt(t(i))/(sqrt(10.0_qp)+sqrt(t(i))))
        x(i) = r(1)

        r = r + rk4( r, t(i), dt )
        write(1,*) t(i), x(i), farmaco(i)
        print*,    t(i), x(i), farmaco(i)
    end do
    close(1)
    call system('gnuplot -c logkill-H.gplot')
!**********************************************************************
contains
!**********************************************************************
    pure function f(r, t) ! Aqui se colocan las ecuaciones a resol
        real(qp), intent(in) :: r(N_equ) ! Valores
        real(qp), intent(in) :: t    ! Paso
        real(qp)             :: far
        real(qp)             :: f(N_equ)
        real(qp)             :: u
        far = 30.0_qp*(sqrt(t)/(sqrt(10.0_qp)+sqrt(t)))
        u = r(1)

        f(1) = g*u-d*u*far

    end function f
!**********************************************************************
    pure function rk4(r, t, dt)                         ! Runge-Kutta 4
        real(qp), intent(in) :: r(N_equ) ! Valores
        real(qp), intent(in) :: t    ! Paso
        real(qp), intent(in) :: dt   ! Tamano de paso
        real(qp)             :: rk4(N_equ)
        real(qp)             :: k1(N_equ), k2(N_equ)
        real(qp)             :: k3(N_equ), k4(N_equ)

        k1 = dt * f( r              , t               )
        k2 = dt * f( r + 0.5_qp * k1, t + 0.5_qp * dt )
        k3 = dt * f( r + 0.5_qp * k2, t + 0.5_qp * dt )
        k4 = dt * f( r + k3         , t + dt          )

        rk4 = ( k1 + ( 2.0_qp * k2 ) + ( 2.0_qp * k3 ) + k4 ) / 6.0_qp
    end function rk4
!**********************************************************************
end program logkillH
