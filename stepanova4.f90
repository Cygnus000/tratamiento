program stepanova
    use, intrinsic :: iso_fortran_env, only: qp=>real128
    implicit none

    real(qp), parameter :: x0 = 1._qp
    real(qp), parameter :: z0 = 0._qp
    real(qp), parameter :: g = 0.6_qp
    real(qp), parameter :: x_max = 7.5_qp
    real(qp), parameter :: gama = 1.0_qp
    real(qp), parameter :: kI = 0.5_qp
    real(qp), parameter :: beta = 0.3_qp
    real(qp), parameter :: delta = 0.4_qp
    real(qp), parameter :: mu = 0.1_qp

    real(qp), parameter :: t0 = 0.0_qp
    real(qp), parameter :: tmax = 100.0_qp
    integer          , parameter :: N    = 10000
    real(qp), parameter :: dt = (tmax - t0) / dble(N)
    integer          , parameter :: N_equ = 2    ! Numero de ecuaciones

    integer           :: i
    real(qp) :: r(N_equ)
    real(qp) :: t(N), x(N), z(N)
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
        real(qp), intent(in) :: r(N_equ) ! Valores
        real(qp), intent(in) :: t    ! Paso
        real(qp)             :: f(N_equ)
        real(qp)             :: u,v

        u = r(1)
        v = r(2)

        f(1) = g * u * (1.0_qp - u/x_max) - gama * u * v
        f(2) = kI * ( u - beta * u**2 ) * v - delta * v + mu
        
    end function f
!**********************************************************************
    pure function rk4(r, t, dt)                         ! Runge-Kutta 4
        real(qp), intent(in) :: r(N_equ) ! Valores
        real(qp), intent(in) :: t    ! Paso
        real(qp), intent(in) :: dt   ! Tamano de paso
        real(qp)             :: rk4(N_equ)
        real(qp)             :: k1(N_equ), k2(N_equ)
        real(qp)             :: k3(N_equ), k4(N_equ)   

        k1 = dt * f( r,            t                )
        k2 = dt * f( r + 0.5_qp * k1, t + 0.5_qp * dt )
        k3 = dt * f( r + 0.5_qp * k2, t + 0.5_qp * dt )
        k4 = dt * f( r + k3        , t + dt         )

        rk4 = ( k1 + ( 2.0_qp * k2 ) + ( 2.0_qp * k3 ) + k4 ) / 6.0_qp
    end function rk4
!**********************************************************************
end program stepanova
