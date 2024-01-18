program logkillResistencia
    use, intrinsic :: iso_fortran_env, only: qp=>real128
    implicit none

    real(qp), parameter :: x0 = 30._qp
    real(qp), parameter :: g = 0.1_qp
    real(qp), parameter :: d = 0.04_qp
    real(qp), parameter :: lambda = 0.1_qp

    real(qp), parameter :: t0 = 0._qp
    real(qp), parameter :: tmax = 50._qp
    integer , parameter :: N    = 10000
    real(qp), parameter :: dt = (tmax - t0) / dble(N)
    integer , parameter :: N_equ = 1    ! Numero de ecuaciones

    integer           :: i
    real(qp) :: r(N_equ), t(N), x(N), farmaco(N)
!**********************************************************************
    t = [ ( dt * i, i = 1, N ) ]             ! llenando vector temporal
!**********************************************************************
    r = [ x0 ]                                      ! valores iniciales
!**********************************************************************
    do i = 1, N                                           ! resolviendo
        farmaco(i) = 30._qp*(sqrt(t(i))/(sqrt(10._qp)+sqrt(t(i))))
        x(i) = r(1)

        r = r + rk4( r, t(i), farmaco(i), dt )
!**********************************************************************
        open(1,file='resistencia.dat')               ! llenando archivo
          write(1,*) t(i), x(i), farmaco(i)
          print*,    t(i), x(i), farmaco(i)
    end do
    close(1) 
    call system('gnuplot -c resistencia.p')
!**********************************************************************
contains
!**********************************************************************
    pure function f(r, t, far) ! Aqui se colocan las ecuaciones a resol
        real(qp), intent(in) :: r(N_equ) ! Valores
        real(qp), intent(in) :: t    ! Paso
        real(qp), intent(in) :: far
        real(qp)             :: f(N_equ)
        real(qp)             :: u

        u = r(1)

        f(1) = g*u-exp(-lambda*t)*d*u*far
        
    end function f
!**********************************************************************
    pure function rk4(r, t, far, dt)                    ! Runge-Kutta 4
        real(qp), intent(in) :: r(N_equ) ! Valores
        real(qp), intent(in) :: t    ! Paso
        real(qp), intent(in) :: far  !funcion cantidad farmaco
        real(qp), intent(in) :: dt   ! Tamano de paso
        real(qp)             :: rk4(N_equ)
        real(qp)             :: k1(N_equ), k2(N_equ)
        real(qp)             :: k3(N_equ), k4(N_equ)   

        k1 = dt * f( r,            t               ,  far )
        k2 = dt * f( r + 0.5_qp * k1, t + 0.5_qp * dt,  far )
        k3 = dt * f( r + 0.5_qp * k2, t + 0.5_qp * dt,  far )
        k4 = dt * f( r + k3        , t + dt        ,  far )

        rk4 = ( k1 + ( 2._qp * k2 ) + ( 2._qp * k3 ) + k4 ) / 6._qp
    end function rk4
!**********************************************************************
end program logkillResistencia
