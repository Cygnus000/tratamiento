PROGRAM angiogenesis
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: qp=>REAL128
    IMPLICIT NONE

    REAL(qp), PARAMETER :: X0 = 17._qp
    REAL(qp), PARAMETER :: Y0 = 30.7_qp
    REAL(qp), PARAMETER :: G = 5.63_qp
    REAL(qp), PARAMETER :: B = 0.829_qp
    REAL(qp), PARAMETER :: ALFA  = 5.5_qp
    REAL(qp), PARAMETER :: DELTA = 0.2_qp
    REAL(qp), PARAMETER :: GAMA  = 0.4_qp

    REAL(qp), PARAMETER :: T0 = 0._qp
    REAL(qp), PARAMETER :: Tmax = 50._qp
    INTEGER , PARAMETER :: N    = 10000
    REAL(qp), PARAMETER :: Dt = (Tmax - T0) / dble(N)
    INTEGER , PARAMETER :: N_equ = 2    ! Numero de ecuaciones

    INTEGER  :: i
    REAL(qp) :: r(N_equ), t(N), x(N), y(N)
!**********************************************************************
    t = [ ( dt * i, i = 1, N ) ]             ! llenando vector temporal
!**********************************************************************
    r = [ X0, Y0 ]                                  ! valores iniciales
!**********************************************************************
    OPEN(1,FILE='angiogenesis.dat')                  ! llenando archivo
    DO i = 1, N                                           ! resolviendo
        x(i) = r(1)
        y(i) = r(2)

        r = r + rk4( r, t(i), dt )
    WRITE(1,*) t(i), x(i), y(i)
    PRINT*,    t(i), x(i), y(i)
    END DO
    CLOSE(1)
    CALL SYSTEM('gnuplot -c angiogenesis.gplot')
    STOP
!**********************************************************************
CONTAINS
!**********************************************************************
    PURE FUNCTION f(r, t)   ! Aqui se colocan las ecuaciones a resolver
        REAL(qp), INTENT(IN) :: r(N_equ) ! Valores
        REAL(qp), INTENT(IN) :: t    ! Paso
        REAL(qp)             :: f(N_equ)
        REAL(qp)             :: u,v

        u = r(1)
        v = r(2)

        f(1) = G*u*LOG(v/u) - ALFA*EXP(-DELTA*t)*u
        f(2) = B*u**(2._qp/3._qp) - GAMA*u

    END FUNCTION f
!**********************************************************************
    PURE FUNCTION rk4(r, t, dt)                    ! Runge-Kutta 4
        REAL(qp), INTENT(IN) :: r(N_equ) ! Valores
        REAL(qp), INTENT(IN) :: t    ! Paso
        REAL(qp), INTENT(IN) :: dt   ! Tamano de paso
        REAL(qp)             :: rk4(N_equ)
        REAL(qp)             :: k1(N_equ), k2(N_equ)
        REAL(qp)             :: k3(N_equ), k4(N_equ)

        k1 = dt * f( r              , t               )
        k2 = dt * f( r + 0.5_qp * k1, t + 0.5_qp * dt )
        k3 = dt * f( r + 0.5_qp * k2, t + 0.5_qp * dt )
        k4 = dt * f( r + k3         , t + dt          )

        rk4 = ( k1 + ( 2._qp * k2 ) + ( 2._qp * k3 ) + k4 ) / 6._qp

    END FUNCTION rk4
!**********************************************************************
END PROGRAM angiogenesis
