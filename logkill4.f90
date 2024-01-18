program logkill
    use, intrinsic :: iso_fortran_env, only: qp=>real128
    implicit none

    real(qp), parameter :: x0 = 30.0_qp
    real(qp), parameter :: g = 0.1_qp
    real(qp), parameter :: d = 0.4_qp

    real(qp), parameter :: t0 = 0.0_qp
    real(qp), parameter :: tmax = 100.0_qp
    integer , parameter :: N  = 10000
    real(qp), parameter :: dt = (tmax - t0) / dble(N)
    integer , parameter :: N_equ = 1    ! Numero de ecuaciones

    integer  :: i
    real(qp) :: r(N_equ), t(N), x(N)
!**********************************************************************
    t = [ ( dt * i, i = 1, N ) ]             ! llenando vector temporal
!**********************************************************************
    r = [ x0 ]                                      ! valores iniciales
!**********************************************************************
    do i = 1, N                                           ! resolviendo
        x(i) = r(1)

        r = r + rk4( r, t(i), dt )
!**********************************************************************
        open(1,file='logkill.dat')                   ! llenando archivo
          write(1,*) t(i), x(i), 1
          print*,    t(i), x(i), 1
    end do
    close(1) 
    call system('gnuplot -c logkill.p')
!**********************************************************************
contains
!**********************************************************************
    pure function f(r, t)   ! Aqui se colocan las ecuaciones a resolver
        real(qp), intent(in) :: r(N_equ) ! Valores
        real(qp), intent(in) :: t    ! Paso
        real(qp)             :: f(N_equ)
        real(qp)             :: u

        u = r(1)

        f(1) = g*u-d*u
        
    end function f
!**********************************************************************
    pure function rk4(r, t, dt)                         ! Runge-Kutta 4
        real(qp), intent(in) :: r(N_equ) ! Valores
        real(qp), intent(in) :: t    ! Paso
        real(qp), intent(in) :: dt   ! Tamano de paso
        real(qp)             :: rk4(N_equ)
        real(qp)             :: k1(N_equ), k2(N_equ)
        real(qp)             :: k3(N_equ), k4(N_equ)   

        k1 = dt * f( r            , t             )
        k2 = dt * f( r + 0.50 * k1, t + 0.50 * dt )
        k3 = dt * f( r + 0.50 * k2, t + 0.50 * dt )
        k4 = dt * f( r + k3       , t + dt        )

        rk4 = ( k1 + ( 2.00 * k2 ) + ( 2.00 * k3 ) + k4 ) / 6.00
    end function rk4
!**********************************************************************
end program logkill
