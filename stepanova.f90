program stepanova
!Resolucion del modelo de stepanova
implicit none
integer, parameter :: dp = 8
real(dp) t0,tmax,dt,x0,z0,gama,kI,beta,delta,mu,kg,xm
real(dp), allocatable, dimension (:) :: t,x,z
integer i,j,N

!**********************************************************************
t0   = 0.0_dp
tmax = 100.0_dp
N    = 10000
x0   = 1.0_dp
z0 = 0.0_dp
kg = 0.6_dp
xm = 7.5_dp
gama = 1.0_dp
kI = 0.5_dp
beta = 0.3_dp
delta = 0.4_dp
mu = 0.1_dp

allocate(t(0:N),x(0:N),z(0:N))
!**********************************************************************
dt = (tmax - t0) / dble(N) !llenando vector temporal
do i=0,N
  t(i) = t0 + dt * dble(i)
end do
!**********************************************************************
 x(0) = x0 !valores iniciales
 z(0) = z0
!**********************************************************************
do i=1,N !runge kutta
  do j=1,2
    if (j.eq.1) then
      x(i) = x(i-1) + ( kg * x(i-1) * (1.0 - x(i-1)/xm) - gama * x(i-1) * z(i-1) ) * dt
      z(i) = z(i-1) + ( kI * ( x(i-1) - beta * x(i-1)**2 ) * z(i-1) - delta * z(i-1) + mu ) * dt
    else
      x(i) = 0.5_dp * ( x(i-1) + ( kg * x(i-1) * (1.0 - x(i-1)/xm) - gama * x(i-1) * z(i-1) ) * dt + x(i) )
      z(i) = 0.5_dp * ( z(i-1) + ( kI * ( x(i-1) - beta * x(i-1)**2 ) * z(i-1) - delta * z(i-1) + mu ) * dt + z(i) )
    end if
  end do
end do
!**********************************************************************
open(1,file='stepanova.dat') !llenando archivo
do i=0,N,1
  write(1,*) t(i),x(i),z(i)
end do
close(1)
call system('gnuplot -c stepanova.gplot')
!**********************************************************************
end program stepanova
