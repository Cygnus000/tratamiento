program logkillH
!Resolucion de la ecuacion diferencial X'(t)=kg*X-kd*X con runge kutta de orden 2
implicit none
integer, parameter :: dp = 8
real(dp) t0,tmax,dt,x0,k_g,k_d,u0
real(dp), allocatable, dimension (:) :: t,x,z,u
integer i,j,N

!**********************************************************************
t0   = 0.0_dp
tmax = 100.0_dp
N    = 10000
x0   = 30.0_dp
k_g = 0.1_dp
k_d = 0.04_dp
u0  = 0.0_dp

allocate(t(0:N),x(0:N),z(0:N),u(0:N))
!**********************************************************************
dt = (tmax - t0) / dble(N) !llenando vector temporal
do i=0,N
  t(i) = t0 + dt * dble(i)
end do
!**********************************************************************
 x(0) = x0 !valores iniciales
 u(0) = u0
!**********************************************************************
do i=1,N !runge kutta
  do j=1,2
    if (j.eq.1) then
      z(i) = z(i-1) + 1.0_dp * dt
      u(i) = 30.0_dp*(sqrt(z(i-1))/(sqrt(10.0_dp)+sqrt(z(i-1))))
      x(i) = x(i-1) + ( k_g * x(i-1) - k_d * x(i-1) * u(i) ) * dt
    else
      z(i) = 0.5_dp * ( z(i-1) + (1.0_dp ) * dt + z(i) )
      u(i) = 30.0_dp*(sqrt(z(i-1))/(sqrt(10.0_dp)+sqrt(z(i-1))))
      x(i) = 0.5_dp * ( x(i-1) + ( k_g * x(i-1) - k_d * x(i-1) * u(i) ) * dt + x(i) )
    end if
  end do
end do
!**********************************************************************
open(1,file='logkill-H.dat') !llenando archivo
do i=0,N,1
  write(1,*) t(i),x(i),u(i)
end do
close(1)
call system('gnuplot -c logkill-H.gplot')
!**********************************************************************
end program logkillH
