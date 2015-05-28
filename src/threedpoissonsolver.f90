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
 do while((errnorm.gt.tol).and.(N.lt.Nmax))
  N = N+1
  call threedjacobi(dxsq,f,w,errnorm)
 end do
end subroutine

subroutine threedjacobi(dxsq,f,w,errnorm)
 real, intent(out) :: errnorm
 real, intent(in) :: dxsq(3)
 real, intent(inout) :: w(:,:,:)
 real, intent(in) :: f(:,:,:)
 integer :: nx,ny,nz
 integer :: i,j,k

 if(size(w).ne.size(f)) then
  print *, "jacobi: ERROR Forcing and value unequal size"
 end if
 nx = size(f,1)
 ny = size(f,2)
 nz = size(f,3)
 errnorm = 0
 do i=2,nx-1
  do j=2,ny-1
   do k=2,nz-1
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
