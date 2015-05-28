program idealfluid
 use threedderivatives
 use threedpoissonsolver
 integer :: N(3)
 real :: L(3)
 real :: dx(3)
 real, parameter :: pi=4.0*atan(1.0)
 real, allocatable :: x(:),y(:),z(:)
 real, allocatable :: u(:,:,:), du(:,:,:), duold(:,:,:)
 real, allocatable :: v(:,:,:), dv(:,:,:), dvold(:,:,:)
 real, allocatable :: w(:,:,:), dw(:,:,:), dwold(:,:,:)
 real, allocatable :: p(:,:,:), lapP(:,:,:)
 integer :: i,j,k,m,Nt
 real :: dt
! 
 N = (/16,16,16/)
 L = (/1,1,1/)
 Nt = 512
 dt = 0.001
 do i=1,3
  dx(i) = L(i)/N(i)
 end do
 call makecoordinates(N,L)
 call makeflowfield(N)

 call fillflow
 do m=1,Nt
  call stepfields(dt)
  call restepfields(dt)
 end do
 call dumpflowfield
  
 call unmakecoordinates
 call unmakeflowfields

contains
 
 subroutine stepfields(dt)
  real, intent(in) :: dt
  integer :: i,j,k
  do i=2,N(1)-1
   do j=2,N(2)-1
    do k=2,N(3)-1
     du(i,j,k) = -u(i,j,k)*ddx(u,dx,i,j,k)-v(i,j,k)*ddy(u,dx,i,j,k)-w(i,j,k)*ddz(u,dx,i,j,k)
     dv(i,j,k) = -u(i,j,k)*ddx(v,dx,i,j,k)-v(i,j,k)*ddy(v,dx,i,j,k)-w(i,j,k)*ddz(v,dx,i,j,k)
     dw(i,j,k) = -u(i,j,k)*ddx(w,dx,i,j,k)-v(i,j,k)*ddy(w,dx,i,j,k)-w(i,j,k)*ddz(w,dx,i,j,k)
    end do
   end do
  end do
!  P = 
  u = u+du*dt
  v = v+dv*dt
  w = w+dw*dt
 end subroutine

 subroutine restepfields(dt)
  real :: dt
  integer :: i,j,k
  duold = du
  dvold = dv
  dwold = dw
  do i=2,N(1)-1
   do j=2,N(2)-1
    do k=2,N(3)-1
     du(i,j,k) = -u(i,j,k)*ddx(u,dx,i,j,k)-v(i,j,k)*ddy(u,dx,i,j,k)-w(i,j,k)*ddz(u,dx,i,j,k)
     dv(i,j,k) = -u(i,j,k)*ddx(v,dx,i,j,k)-v(i,j,k)*ddy(v,dx,i,j,k)-w(i,j,k)*ddz(v,dx,i,j,k)
     dw(i,j,k) = -u(i,j,k)*ddx(w,dx,i,j,k)-v(i,j,k)*ddy(w,dx,i,j,k)-w(i,j,k)*ddz(w,dx,i,j,k)
    end do
   end do
  end do
!  P =
  u = u+(du-duold)*dt/2
  v = v+(dv-dvold)*dt/2
  w = w+(dw-dwold)*dt/2
 end subroutine

 subroutine fillflow
  integer :: i,j,k
  do i=1,N(1)
   do j=1,N(2)
    do k=1,N(3)
     u(i,j,k) = ((L(1)/(2*pi)))*sin(2*pi*x(i)/L(1))*cos(2*pi*y(j)/L(2))
     v(i,j,k) = -((L(2)/(2*pi)))*cos(2*pi*x(i)/L(1))*sin(2*pi*y(j)/L(2))
     
     v(i,j,k) = v(i,j,k)+((L(2)/(2*pi)))*cos(2*pi*z(k)/L(3))*sin(2*pi*y(j)/L(2))
     w(i,j,k) = -((L(3)/(2*pi)))*cos(2*pi*y(j)/L(2))*sin(2*pi*z(k)/L(3))

     u(i,j,k) = u(i,j,k)-((L(1)/(2*pi)))*sin(2*pi*x(i)/L(1))*cos(2*pi*z(k)/L(3))
     w(i,j,k) = w(i,j,k)+((L(3)/(2*pi)))*cos(2*pi*x(i)/L(1))*sin(2*pi*z(k)/L(3))
    end do
   end do
  end do
 end subroutine

 subroutine dumpflowfield
  integer :: i,j,k
  do i=1,N(1)
   do j=1,N(2)
    do k=1,N(3)
     print *,x(i),y(j),z(k),u(i,j,k),v(i,j,k),w(i,j,k),div(u,v,w,dx,i,j,k)
    end do
   end do
  end do
 end subroutine

 subroutine makeflowfield(Nd)
  integer, intent(in) :: Nd(3)
  allocate(u(Nd(1),Nd(2),Nd(3)))
  u(:,:,:) = 0
  allocate(du(Nd(1),Nd(2),Nd(3)))
  du(:,:,:) = 0
  allocate(duold(Nd(1),Nd(2),Nd(3)))
  duold(:,:,:) = 0
  allocate(v(Nd(1),Nd(2),Nd(3)))
  v(:,:,:) = 0
  allocate(dv(Nd(1),Nd(2),Nd(3)))
  dv(:,:,:) = 0
  allocate(dvold(Nd(1),Nd(2),Nd(3)))
  dvold(:,:,:) = 0
  allocate(w(Nd(1),Nd(2),Nd(3)))
  w(:,:,:) = 0
  allocate(dw(Nd(1),Nd(2),Nd(3)))
  dw(:,:,:) = 0
  allocate(dwold(Nd(1),Nd(2),Nd(3)))
  dwold(:,:,:) = 0
  allocate(p(Nd(1),Nd(2),Nd(3)))
  p(:,:,:) = 0 
  allocate(lapp(Nd(1),Nd(2),Nd(3)))
  lapp(:,:,:) = 0
end subroutine

 subroutine makecoordinates(Nd,Ld)
  integer, intent(in) :: Nd(3)
  real, intent(in) :: Ld(3)
  integer :: i
  allocate(x(Nd(1)))
  do i=1,N(1)
   x(i) = (i-.5*(Nd(1)+1))*Ld(1)/N(1)
  end do
  allocate(y(Nd(2)))
  do i=1,N(2)
   y(i)=(i-.5*(Nd(2)+1))*Ld(2)/Nd(2)
  end do
  allocate(z(Nd(3)))
  do i=1,Nd(3)
   z(i) = (i-.5*(Nd(3)+1))*Ld(3)/N(3)
  end do
 end subroutine

 subroutine unmakecoordinates
  deallocate(x)
  deallocate(y)
  deallocate(z)
 end subroutine

 subroutine unmakeflowfields()
  deallocate(u)
  deallocate(du)
  deallocate(duold)
  deallocate(v)
  deallocate(dv)
  deallocate(dvold)
  deallocate(w)
  deallocate(dw)
  deallocate(dwold)
  deallocate(P)
  deallocate(lapP)
 end subroutine

 subroutine setBCbox(u,v,w)
  real, intent(in) :: u(:,:,:),v(:,:,:),w(:,:,:)
  integer :: N(3),i
  do i=1,3
   N(i) = size(u,i)
  end do
 end subroutine
end program
