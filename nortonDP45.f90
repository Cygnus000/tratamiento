program norton
    use, intrinsic :: iso_fortran_env, only: qp=>real128
    implicit none
    !valores iniciales
    real(qp), parameter :: x10 = 30.0q0
    real(qp), parameter :: x20 = 0.0q0
    real(qp), parameter :: g = 0.1q0
    real(qp), parameter :: d = 0.04q0
    real(qp), parameter :: a = 0.1q0
    real(qp), parameter :: x_max = 40.0q0
    real(qp), parameter :: t0 = 0.0q0
    real(qp), parameter :: t_max = 100.0q0
    !parametros iniciales
    real(qp), parameter :: tinny = 0.001
    real(qp), parameter :: max_steps = 300
    integer , parameter :: N_equ = 2    ! Numero de ecuaciones
    real(qp) :: dt = 0.5q0
    real(qp) :: yscal = 0.0q0
    real(qp) :: step = 0.0q0
    real(qp) :: dt_next = 0.0q0
    real(qp) :: t, farmaco, x1, x2
    real(qp) :: r(N_equ), temp(N_equ)
!**********************************************************************
    t = t0                                          ! valores iniciales
    r = [ x10, x20 ]                                
!**********************************************************************
    do while(t < t_max .and. step<max_steps)              ! resolviendo
        step=step+1
        
        if ((t + dt) > t_max) then
               dt = t_max-t
        endif
        
        call adapt(r,t,farmaco,dt,dt_next)
        
        farmaco = 30.0q0*(sqrt(t)/(sqrt(10.0q0)+sqrt(t)))
        x1 = r(1)
        x2 = r(2)
        
        open(1,file='norton.dat')             ! llenando archivo
        write(1,*) t, x1+x2, x1, x2, farmaco
        print*,    t, x1+x2, x1, x2, farmaco
        
    end do
    close(1) 
    call system('gnuplot -c norton.p')
!**********************************************************************
contains
!**********************************************************************
    pure function f(r, t, far) ! Aqui se colocan las ecuaciones a resol
        real(qp), intent(in) :: r(N_equ) ! Valores
        real(qp), intent(in) :: t    ! Paso
        real(qp), intent(in) :: far
        real(qp)             :: f(N_equ)
        real(qp)             :: u,v

        u = r(1)
        v = r(2)

        f(1) = g*u*log(x_max/u)-a*far*u
        f(2) = a*far*u-d*v
        
    end function f
!**********************************************************************
    function dp45(r, t, far, dt, error_global)
    real(qp), intent(in)    :: r(N_equ) ! Valores
    real(qp), intent(in)    :: t    ! Paso
    real(qp), intent(in)    :: far  !funcion cantidad farmaco
    real(qp), intent(in)    :: dt   ! Tamano de paso
    real(qp), intent(inout) :: error_global
    real(qp) :: xerr(N_equ), dp45(N_equ)
    real(qp) :: k1(N_equ), k2(N_equ), k3(N_equ), k4(N_equ)
    real(qp) :: k5(N_equ), k6(N_equ), k7(N_equ)
    ! parametros de tiempo
    real(qp),parameter :: a2 =  1.0_qp / 5.0_qp
    real(qp),parameter :: a3 =  3.0_qp / 10.0_qp
    real(qp),parameter :: a4 =  4.0_qp / 5.0_qp
    real(qp),parameter :: a5 =  8.0_qp / 9.0_qp
    ! parametros de paso intermedio
    real(qp),parameter :: b21 = 1.0_qp / 5.0_qp
    real(qp),parameter :: b31 = 3.0_qp / 40.0_qp
    real(qp),parameter :: b32 = 9.0_qp / 40.0_qp
    real(qp),parameter :: b41 = 44.0_qp / 45.0_qp
    real(qp),parameter :: b42 = -56.0_qp / 15.0_qp
    real(qp),parameter :: b43 = 32.0_qp / 9.0_qp
    real(qp),parameter :: b51 = 19372.0_qp / 6561.0_qp
    real(qp),parameter :: b52 = -25360.0_qp / 2187.0_qp
    real(qp),parameter :: b53 = 64448.0_qp / 6561.0_qp
    real(qp),parameter :: b54 = -212.0_qp / 729.0_qp
    real(qp),parameter :: b61 = 9017.0_qp / 3168.0_qp
    real(qp),parameter :: b62 = -355.0_qp / 33.0_qp
    real(qp),parameter :: b63 = 46732.0_qp / 5247.0_qp
    real(qp),parameter :: b64 = 49.0_qp / 176.0_qp
    real(qp),parameter :: b65 = -5103.0_qp / 18656.0_qp
    real(qp),parameter :: b71 = 35.0_qp / 384.0_qp
    real(qp),parameter :: b73 = 500.0_qp / 1113.0_qp
    real(qp),parameter :: b74 = 125.0_qp / 192.0_qp
    real(qp),parameter :: b75 = -2187.0_qp / 6784.0_qp
    real(qp),parameter :: b76 = 11.0_qp / 84.0_qp
    ! parametros rk4
    real(qp),parameter :: c1  = b71
    real(qp),parameter :: c3  = b73
    real(qp),parameter :: c4  = b74
    real(qp),parameter :: c5  = b75
    real(qp),parameter :: c6  = b76
    ! parametros rk5
    real(qp),parameter :: d1  = 5179.0_qp / 57600.0_qp
    real(qp),parameter :: d3  = 7571.0_qp / 16695.0_qp
    real(qp),parameter :: d4  = 393.0_qp / 640.0_qp
    real(qp),parameter :: d5  = -92097.0_qp / 339200.0_qp
    real(qp),parameter :: d6  = 187.0_qp / 2100.0_qp
    real(qp),parameter :: d7  = 1.0_qp / 40.0_qp
    ! parametros error
    real(qp),parameter :: e1  = 71.0_qp/57600.0_qp     !c1 - d1
    real(qp),parameter :: e3  =-71.0_qp/16695.0_qp     !c3 - d3
    real(qp),parameter :: e4  = 71.0_qp/1920.0_qp      !c4 - d4
    real(qp),parameter :: e5  =-17253.0_qp/339200.0_qp !c5 - d5
    real(qp),parameter :: e6  = 22.0_qp/525.0_qp       !c6 - d6
    real(qp),parameter :: e7  =-1.0_qp/40.0_qp         !   - d7

        k1 = dt*f(r                                                        ,t        ,far)
        k2 = dt*f(r + (b21*k1                                             ),t + a2*dt,far)
        k3 = dt*f(r + (b31*k1 + b32*k2                                    ),t + a3*dt,far)
        k4 = dt*f(r + (b41*k1 + b42*k2 + b43*k3                           ),t + a4*dt,far)
        k5 = dt*f(r + (b51*k1 + b52*k2 + b53*k3 + b54*k4                  ),t + a5*dt,far)
        k6 = dt*f(r + (b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5         ),t + dt   ,far) ! paso calculado para rk4
        k7 = dt*f(r + (b71*k1 +          b73*k3 + b74*k4 + b75*k5 + b76*k6),t + dt   ,far) ! paso calculado para rk5

        dp45 = ( c1*k1 + c3*k3 + c4*k4 + c5*k5 + c6*k6 ) !rk4
        xerr = abs( e1*k1 + e3*k3 + e4*k4 + e5*k5 + e6*k6 + e7*k7 ) ! rk5-rk4
        error_global = maxval(xerr)

    end function dp45
!**********************************************************************
    subroutine adapt(r,t,far,dt,dt_next)
    real(qp), intent(inout) :: r(N_equ) ! Valores
    real(qp), intent(inout) :: t    ! Paso
    real(qp), intent(in)    :: far  !funcion cantidad farmaco
    real(qp), intent(inout) :: dt   ! Tamano de paso
    real(qp), intent(inout) :: dt_next
    real(qp) :: yscal = 0.0q0
    real(qp), parameter :: safety = 0.9q0
    real(qp), parameter :: econ = 0.000189q0
    real(qp), parameter :: eps = 0.00005q0
    real(qp) :: error_global = 0.0q0
    real(qp) :: ytemp(N_equ)
    real(qp) :: emax
    real(qp) :: dttemp
        
        do
            ytemp = dp45(r,t,far,dt,error_global)
            yscal = norm2(r)+norm2(dt*f(r, t, farmaco))+tinny
            emax = abs(error_global/yscal/eps)
            if (emax <= 1.0q0) exit
            dttemp=safety*dt*emax**(-0.25q0)
            dt=max(abs(dttemp),0.25q0*abs(dt))
            t=t+dt
        enddo
        if (emax > econ) then
            dt_next=safety*dt*(emax**(-0.2q0))
        else
            dt_next=5.0q0*dt
        endif
        t=t+dt
        r=r+ytemp
        
    end subroutine adapt
!**********************************************************************
end program norton
