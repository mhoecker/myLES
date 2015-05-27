program test1dpoisson
 use onedpoissonsolver
 implicit none
 integer, parameter :: Nx = 1024
 real :: f(Nx)
 real :: y(Nx)
 real :: x(Nx)
 real, parameter :: dx=.01
 real, parameter :: dxsq=dx*dx
 real,parameter :: tol=1e-25
 real :: err
 integer, parameter :: Nmax=Nx*Nx
 integer :: i
 do i=1,Nx
  x(i) = dx*i-Nx*dx/2
  y(i) = 0
  f(i) = x(i)-Nx*dx/2
 end do
  y(1) = x(1)
  y(Nx) = x(Nx)
 call onedpoisson(dxsq,f,y,err,tol,Nmax)
 do i=1,Nx
  if((i.lt.2).or.(i.gt.Nx-1)) then
  print *, x(i),y(i),f(i),f(i)
  elseif((i.eq.2).or.(i.eq.Nx-1)) then
   print *,x(i),y(i),f(i),threepointder2(y,dxsq,i)
  else
!   print *,x(i),y(i),f(i),threepointder2(y,dxsq,i)
   print *,x(i),y(i),f(i),fivepointder2(y,dxsq,i)
  end if
 end do
end program
