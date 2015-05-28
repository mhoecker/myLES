module threedpoissonsolver
use threedderivatives
contains
! Assumes end points are fixed values
! solves d^2w/dx^2+d^2w/dy^2+d^2w/dz^2 = f
! uses three and five point stencils to solve the DE
subroutine threedpoisson(dxsq,f,w,errnorm,tol,Nmax)
 implicit none
 integer, intent(in) :: Nmax
 real, intent(in) :: tol
 real, intent(out) :: errnorm
 real, intent(in) :: dxsq(3)
 real, intent(inout) :: w(:,:,:)
 real, intent(in) :: f(:,:,:)
 integer :: N
 errnorm = 2*tol
 N=0
 if(size(w).ne.size(f)) then
  print *, "jacobi: ERROR Forcing and value unequal size"
 end if
 do while((errnorm.gt.tol).and.(N.lt.Nmax))
  N = N+1
  call rescalejacobi(dxsq,f,w,errnorm)
 end do
end subroutine

recursive subroutine rescalejacobi(dxsq,f,w,errnorm)
 real, intent(out) :: errnorm
 real, intent(in) :: dxsq(3)
 real, intent(inout) :: w(:,:,:)
 real, intent(in) :: f(:,:,:)
 real, allocatable :: v(:,:,:), g(:,:,:)
 real :: dysq(3)
 integer :: n(3)
 integer :: m(3)
 integer :: i
! Find the array sizes
 do i=1,3
  n(i) = size(f,i)
  call roughhalf(n(i),dxsq(i),m(i),dysq(i))
 end do
 allocate(v(m(1),m(2),m(3)))
 allocate(g(m(1),m(2),m(3)))
 call interp3(f,g)
 call interp3(w,v)
 if((n(1).ne.m(1)).or.(n(2).ne.m(2)).or.(n(3).ne.m(3))) then
  call rescalejacobi(dysq,g,v,errnorm)
  call interp3(v,w)
  call threedjacobi(dxsq,f,w,errnorm)
  call interp3(w,v)
  call rescalejacobi(dysq,g,v,errnorm)
 else
  call threedjacobi(dysq,g,v,errnorm)
 end if
 call interp3(v,w)
 call threedjacobi(dxsq,f,w,errnorm)
 deallocate(v)
 deallocate(g)
end subroutine

subroutine interp3(a,b)
 real, intent(in) :: a(:,:,:)
 real, intent(out) :: b(:,:,:)
 integer :: n(3)
 integer :: m(3)
 integer :: i,j,k,l, t(3,2)
 real :: w(3,2)
 do i=1,3
  n(i) = size(a,i)
  m(i) = size(b,i)
 end do
 do k=1,m(3)
  call interpparams(k,m(3),n(3),t(3,:),w(3,:))
  do j=1,m(2)
   call interpparams(j,m(2),n(2),t(2,:),w(2,:))
   do i=1,m(1)
    call interpparams(i,m(1),n(1),t(1,:),w(1,:))
    b(i,j,k) = 0
    b(i,j,k) = b(i,j,k)+w(1,1)*w(2,1)*w(3,1)*a(t(1,1),t(2,1),t(3,1))
    b(i,j,k) = b(i,j,k)+w(1,1)*w(2,1)*w(3,2)*a(t(1,1),t(2,1),t(3,2))
    b(i,j,k) = b(i,j,k)+w(1,1)*w(2,2)*w(3,1)*a(t(1,1),t(2,2),t(3,1))
    b(i,j,k) = b(i,j,k)+w(1,1)*w(2,2)*w(3,2)*a(t(1,1),t(2,2),t(3,2))
    b(i,j,k) = b(i,j,k)+w(1,2)*w(2,1)*w(3,1)*a(t(1,2),t(2,1),t(3,1))
    b(i,j,k) = b(i,j,k)+w(1,2)*w(2,1)*w(3,2)*a(t(1,2),t(2,1),t(3,2))
    b(i,j,k) = b(i,j,k)+w(1,2)*w(2,2)*w(3,1)*a(t(1,2),t(2,2),t(3,1))
    b(i,j,k) = b(i,j,k)+w(1,2)*w(2,2)*w(3,2)*a(t(1,2),t(2,2),t(3,2))
   end do
  end do
 end do
end subroutine

subroutine interpparams(i,m,n,t,w)
 integer, intent(in) :: i,m,n
 integer, intent(out) :: t(2)
 real, intent(out) :: w(2)
  if(m.ne.n) then
   w(1) = (i-1.0)*(n-1.0)/(m-1.0)
   t(1) = 1+floor(w(1))
   t(2) = 1+ceiling(w(1))
   w(2) = w(1)-floor(w(1))
   w(1) = 1-w(2)
  else
   w(1) = 1
   t(1) = i
   w(2) = 0
   t(2) = i
  end if
  if((t(1).lt.1).or.(t(1).gt.n)) then
   print *, i,m,n,t(1),t(2),w(1),w(2)
  elseif((t(2).lt.1).or.(t(2).gt.n)) then
   print *, i,m,n,t(1),t(2),w(1),w(2)
  end if
end subroutine

subroutine roughhalf(n,dx2,m,dy2)
 integer, intent(in) :: n
 real, intent(in) :: dx2
 integer, intent(out) :: m
 real, intent(out) :: dy2
 if((modulo(n,2).eq.0).and.(n.gt.4)) then
  m = 1+n/2
 elseif((modulo(n,2).eq.1).and.(n.gt.3)) then
  m = (n+1)/2
 else
  m = n
 end if
 dy2 = dx2*((n-1.0)/(m-1.0))**2
end subroutine

subroutine threedjacobi(dxsq,f,w,errnorm)
 real, intent(out) :: errnorm
 real, intent(in) :: dxsq(3)
 real, intent(inout) :: w(:,:,:)
 real, intent(in) :: f(:,:,:)
 integer :: nx,ny,nz
 integer :: i,j,k

 nx = size(f,1)
 ny = size(f,2)
 nz = size(f,3)
 errnorm = 0
 do k=2,nz-1
  do j=2,ny-1
   do i=2,nx-1
    w(i,j,k) = jacobi(f,w,dxsq,i,j,k)
    errnorm = errnorm+abs(f(i,j,k)-laplace(w,dxsq,i,j,k))
   end do
  end do
 end do
end subroutine

real function jacobi(f,w,dxsq,i,j,k)
 real, intent(in) :: f(:,:,:)
 real, intent(in) :: w(:,:,:)
 real, intent(in) :: dxsq(3)
 integer, intent(in) :: i,j,k
 integer :: nx,ny,nz
 logical :: onedge,nearedge
 nx = size(f,1)
 ny = size(f,2)
 nz = size(f,3)
! Test if on an edge
 onedge = (i.eq.1).or.(i.eq.nx)
 onedge = onedge.or.(j.eq.1).or.(j.eq.ny)
 onedge = onedge.or.(k.eq.1).or.(k.eq.nz)
! test if near an edge
 nearedge = (i.eq.2).or.(i.eq.nx-1)
 nearedge = nearedge.or.(j.eq.2).or.(j.eq.ny-1)
 nearedge = nearedge.or.(k.eq.2).or.(k.eq.nz-1)
! Don’t change edge points
 if(onedge) then
  jacobi = w(i,j,k)
! Near Edge
 else if(nearedge) then
  jacobi = -f(i,j,k)*dxsq(1)*dxsq(2)*dxsq(3)
  jacobi = jacobi+(w(i-1,j  ,k  )+w(i+1,j  ,k  ))*dxsq(2)*dxsq(3)
  jacobi = jacobi+(w(i  ,j-1,k  )+w(i  ,j+1,k  ))*dxsq(1)*dxsq(3)
  jacobi = jacobi+(w(i  ,j  ,k-1)+w(i  ,j  ,k+1))*dxsq(1)*dxsq(2)
  jacobi = jacobi/(2*(dxsq(1)*dxsq(2)+dxsq(2)*dxsq(3)+dxsq(3)*dxsq(1)))
! Interior
 else
  jacobi = -6*f(i,j,k)*dxsq(1)*dxsq(2)*dxsq(3)
  jacobi = jacobi+(16*(w(i-1,j,k)+w(i+1,j,k))-w(i-2,j,k)-w(i+2,j,k))*dxsq(2)*dxsq(3)
  jacobi = jacobi+(16*(w(i,j-1,k)+w(i,j+1,k))-w(i,j-2,k)-w(i,j+2,k))*dxsq(1)*dxsq(2)
  jacobi = jacobi+(16*(w(i,j,k-1)+w(i,j,k+1))-w(i,j,k-2)-w(i,j,k+2))*dxsq(1)*dxsq(2)
  jacobi = jacobi/(30*(dxsq(1)*dxsq(2)+dxsq(2)*dxsq(3)+dxsq(3)*dxsq(1)))
 end if
end function

end module
