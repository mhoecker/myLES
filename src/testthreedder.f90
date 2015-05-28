program testthreedder
 use threedderivatives
 implicit none
 integer, parameter :: N(3)=(/12,12,12/)
 real, parameter :: pi=4.0*atan(1.0)
 real, parameter :: dx(3)=(/2*pi/N(1),2*pi/N(2),2*pi/N(3)/)
 real, parameter :: dxsq(3)=(/dx(1)**2,dx(2)**2,dx(3)**2/)
 real :: w(N(1),N(2),N(3))
 real :: x(N(1))
 real :: y(N(2))
 real :: z(N(3))
 integer :: i,j,k
 do i=1,N(1)
  x(i) = (i-.5)*2*pi/N(1)
  do j=1,N(2)
   y(j) = (j-.5)*2*pi/N(1)
   do k=1,N(3)
    z(k) = (k-.5)*2*pi/N(1)
    w(i,j,k) = x(i)+ 2*y(j) +3*z(k)
   end do
  end do
 end do
 do i=1,N(1)
  do j=1,N(2)
   do k=1,N(3)
    print *, x(i),y(j),z(k),w(i,j,k), &
           ddx(w,dx,i,j,k),ddy(w,dx,i,j,k),ddz(w,dx,i,j,k), &
           ddxsq(w,dxsq,i,j,k),ddysq(w,dxsq,i,j,k),ddzsq(w,dxsq,i,j,k), &
           laplace(w,dxsq,i,j,k)
   end do
   print *
  end do
  print *
 end do
end program
