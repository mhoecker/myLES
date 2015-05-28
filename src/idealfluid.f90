program idealfluid
 use threedderivatives
 use threedpoissonsolver
 integer :: N(3)
 integer :: i,j,k,m,Nt
 real    :: L(3)
 real    :: dx(3), dxsq(3)
 real    :: dt
 real    :: Perr
 real    :: Ptol
 real, parameter   :: pi = 4.0*atan(1.0)
 real, allocatable :: x(:),y(:),z(:)
 real, allocatable :: u(:,:,:), du(:,:,:), duold(:,:,:)
 real, allocatable :: v(:,:,:), dv(:,:,:), dvold(:,:,:)
 real, allocatable :: w(:,:,:), dw(:,:,:), dwold(:,:,:)
 real, allocatable :: p(:,:,:), lapP(:,:,:)
 N = (/64,64,64/)
 L = (/128,128,128/)
 Nt = 10
 dt = 0.1
 Ptol = N(1)*N(2)*N(3)*.01
 do i=1,3
  dx(i) = L(i)/N(i)
  dxsq(i) = dx(i)**2
 end do
 call makecoordinates
 call makeflowfield
 call fillflow

 do m=1,Nt
  call predictcorrect(dt,Ptol,(N(1)+N(2)+N(3))**2,Perr)
 end do

 call dumpflowfield
 call unmakecoordinates
 call unmakeflowfields

contains

 subroutine predictcorrect(dt,tol,Niter,errnorm)
  real, intent(in) :: dt
  real, intent(in) :: tol
  integer, intent(in) :: Niter
  real, intent(out) :: errnorm
  call stepfields(dt,tol,Niter,errnorm)
  call restepfields(dt,tol,Niter,errnorm)
 end subroutine

 subroutine stepfields(dt,tol,Niter,errnorm)
  real, intent(in) :: dt
  real, intent(in) :: tol
  integer, intent(in) :: Niter
  real, intent(out) :: errnorm
  call calculatediffields()
  call PoissonPsolver(errnorm,tol,Niter)
  u = u+du*dt
  v = v+dv*dt
  w = w+dw*dt
 end subroutine

 subroutine restepfields(dt,tol,Niter,errnorm)
  real, intent(in) :: dt
  real, intent(in) :: tol
  integer, intent(in) :: Niter
  real, intent(out) :: errnorm
  duold = du
  dvold = dv
  dwold = dw
  call calculatediffields()
  call PoissonPsolver(errnorm,tol,Niter)
  u = u+(du-duold)*dt/2
  v = v+(dv-dvold)*dt/2
  w = w+(dw-dwold)*dt/2
 end subroutine

 subroutine calculatediffields()
  integer :: i,j,k
  do k=2,N(3)-1
   do j=2,N(2)-1
    do i=2,N(1)-1
     du(i,j,k) = -u(i,j,k)*ddx(u,dx,i,j,k)-v(i,j,k)*ddy(u,dx,i,j,k)-w(i,j,k)*ddz(u,dx,i,j,k)
     dv(i,j,k) = -u(i,j,k)*ddx(v,dx,i,j,k)-v(i,j,k)*ddy(v,dx,i,j,k)-w(i,j,k)*ddz(v,dx,i,j,k)
     dw(i,j,k) = -u(i,j,k)*ddx(w,dx,i,j,k)-v(i,j,k)*ddy(w,dx,i,j,k)-w(i,j,k)*ddz(w,dx,i,j,k)
    end do
   end do
  end do
 end subroutine

 subroutine PoissonPsolver(errnorm,tol,Niter)
  real, intent(in) :: tol
  integer, intent(in) :: Niter
  real, intent(out) :: errnorm
  do k=1,N(3)
   do j=1,N(2)
    do i=1,N(1)
     LapP(i,j,k) = div(du,dv,dw,dx,i,j,k)
    end do
   end do
  end do
  call threedpoisson(dxsq,LapP,P,errnorm,tol,Niter)
  do k=2,N(3)-1
   do j=2,N(2)-1
    do i=2,N(1)-1
     du(i,j,k) = du(i,j,k)-ddx(P,dx,i,j,k)
     dv(i,j,k) = dv(i,j,k)-ddy(P,dx,i,j,k)
     dw(i,j,k) = dw(i,j,k)-ddz(P,dx,i,j,k)
    end do
   end do
  end do
 end subroutine

 subroutine fillflow
  integer :: i,j,k
  do k=1,N(3)
   do j=1,N(2)
    do i=1,N(1)
     u(i,j,k) = u(i,j,k)+((L(1)/(2*pi)))*sin(2*pi*x(i)/L(1))*cos(2*pi*y(j)/L(2))
     v(i,j,k) = v(i,j,k)-((L(2)/(2*pi)))*cos(2*pi*x(i)/L(1))*sin(2*pi*y(j)/L(2))

     v(i,j,k) = v(i,j,k)+((L(2)/(2*pi)))*cos(2*pi*z(k)/L(3))*sin(2*pi*y(j)/L(2))
     w(i,j,k) = w(i,j,k)-((L(3)/(2*pi)))*cos(2*pi*y(j)/L(2))*sin(2*pi*z(k)/L(3))

     u(i,j,k) = u(i,j,k)-((L(1)/(2*pi)))*sin(2*pi*x(i)/L(1))*cos(2*pi*z(k)/L(3))
     w(i,j,k) = w(i,j,k)+((L(3)/(2*pi)))*cos(2*pi*x(i)/L(1))*sin(2*pi*z(k)/L(3))
    end do
   end do
  end do
 end subroutine

 subroutine dumpflowfield
  integer :: i,j,k
  do k=1,N(3)
   do j=1,N(2)
    do i=1,N(1)
     print *,x(i),y(j),z(k),u(i,j,k),v(i,j,k),w(i,j,k),div(u,v,w,dx,i,j,k)
    end do
   end do
  end do
 end subroutine

 subroutine makeflowfield
  allocate(u(N(1),N(2),N(3)))
  u(:,:,:) = 0
  allocate(du(N(1),N(2),N(3)))
  du(:,:,:) = 0
  allocate(duold(N(1),N(2),N(3)))
  duold(:,:,:) = 0
  allocate(v(N(1),N(2),N(3)))
  v(:,:,:) = 0
  allocate(dv(N(1),N(2),N(3)))
  dv(:,:,:) = 0
  allocate(dvold(N(1),N(2),N(3)))
  dvold(:,:,:) = 0
  allocate(w(N(1),N(2),N(3)))
  w(:,:,:) = 0
  allocate(dw(N(1),N(2),N(3)))
  dw(:,:,:) = 0
  allocate(dwold(N(1),N(2),N(3)))
  dwold(:,:,:) = 0
  allocate(p(N(1),N(2),N(3)))
  p(:,:,:) = 0
  allocate(lapp(N(1),N(2),N(3)))
  lapp(:,:,:) = 0
end subroutine

 subroutine makecoordinates
  integer :: i
  allocate(x(N(1)))
  do i=1,N(1)
   x(i) = (i-.5*(N(1)+1))*L(1)/N(1)
  end do
  allocate(y(N(2)))
  do i=1,N(2)
   y(i)=(i-.5*(N(2)+1))*L(2)/N(2)
  end do
  allocate(z(N(3)))
  do i=1,N(3)
   z(i) = (i-.5*(N(3)+1))*L(3)/N(3)
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
