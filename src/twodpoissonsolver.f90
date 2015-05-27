module twodpoissonsolver
contains
! Assumes end points are fixed values
! solves d^2z/dx^2+d^2z/dy^2 = f
! uses three and five  point stencils to solve the DE
subroutine twodpoisson(dxsq,f,z,errnorm,tol,Nmax)
 implicit none
 integer, intent(in) :: Nmax
 real, intent(in) :: tol
 real, intent(out) :: errnorm
 real, intent(in) :: dxsq(2)
 real, intent(inout) :: z(:,:)
 real, intent(in) :: f(:,:)
 integer :: N
 errnorm = 2*tol
 N=0
 do while((errnorm.gt.tol).and.(N.lt.Nmax))
  N = N+1
  call twodjacobi(dxsq,f,z,errnorm)
 end do
end subroutine

subroutine twodjacobi(dxsq,f,z,errnorm)
 real, intent(out) :: errnorm
 real, intent(in) :: dxsq(2)
 real, intent(inout) :: z(:,:)
 real, intent(in) :: f(:,:)
 integer :: nx,ny
 integer :: i,j
 
 if(size(z).ne.size(f)) then
  print *, "onedjacobi: ERROR Forcing and value unequal size"
 end if
 nx = size(z,1)
 ny = size(z,2)
 errnorm = 0
 do i=2,nx-1
  do j=2,ny-1
   z(i,j) = jacobi(f,z,dxsq,i,j)
   errnorm = errnorm+abs(f(i,j)-laplace(z,dxsq,i,j))
  end do
 end do
end subroutine

real function laplace(z,dxsq,i,j)
 real, intent(in) :: z(:,:)
 real, intent(in) :: dxsq(2)
 integer, intent(in) :: i,j
 integer :: nx,ny
 nx = size(z,1)
 ny = size(z,2)
! x derrivative
! use 1 sided 3 point stencil for border points
 if(i.eq.1) then
  laplace = (z(i+2,j)+z(i,j)-2*z(i+1,j))/dxsq(1)
 else if(i.eq.nx) then
  laplace = (z(i,j)+z(i-2,j)-2*z(i-1,j))/dxsq(1)
! use 3 point centered stencil near edge
 else if((i.eq.2).or.(i.eq.nx-1)) then
  laplace = (z(i+1,j)+z(i-1,j)-2*z(i,j))/dxsq(1)
! use 5 point centered stencil for interior points
 else
  laplace = (-30*z(i,j)+16*(z(i+1,j)+z(i-1,j))-z(i+2,j)-z(i-2,j))/(6*dxsq(1))
 end if
! y derrivative
! use 1 sided 3 point stencil for border points
 if(j.eq.1) then
  laplace = laplace+(z(i,j+2)+z(i,j)-2*z(i,j+1))/dxsq(2)
 else if(j.eq.ny) then
  laplace = laplace+(z(i,j)+z(i,j-2)-2*z(i,j-1))/dxsq(2)
! use 3 point centered stencil near edge
 else if((j.eq.2).or.(j.eq.ny-1)) then
  laplace = laplace+(z(i,j+1)+z(i,j-1)-2*z(i,j))/dxsq(2)
! use 5 point centered stencil for interior points
 else
  laplace = laplace+(-30*z(i,j)+16*(z(i,j+1)+z(i,j-1))-z(i,j+2)-z(i,j-2))/(6*dxsq(2))
 end if
end function

real function jacobi(f,z,dxsq,i,j)
 real, intent(in) :: f(:,:)
 real, intent(in) :: z(:,:)
 real, intent(in) :: dxsq(2)
 integer, intent(in) :: i,j
 integer :: nx,ny
 nx = size(z,1)
 ny = size(z,2)
! Don’t change edge points
 if((i.eq.1).or.(i.eq.nx).or.(j.eq.1).or.(j.eq.ny)) then
  jacobi = z(i,j)
! Near Edge
 else if(((i.eq.2).or.(i.eq.nx-1)).or.((j.eq.2).or.(j.eq.ny-1))) then
!  print *, "Corner",i,j
  jacobi = -f(i,j)*dxsq(1)*dxsq(2)
  jacobi = jacobi+(z(i-1,j)+z(i+1,j))*dxsq(2)
  jacobi = jacobi+(z(i,j-1)+z(i,j+1))*dxsq(1)
  jacobi = jacobi/(2*(dxsq(1)+dxsq(2)))
! Interior
 else
  jacobi = -6*f(i,j)*dxsq(1)*dxsq(2)
  jacobi = jacobi+(16*(z(i-1,j)+z(i+1,j))-z(i-2,j)-z(i+2,j))*dxsq(2)
  jacobi = jacobi+(16*(z(i,j-1)+z(i,j+1))-z(i,j-2)-z(i,j+2))*dxsq(1)
  jacobi = jacobi/(30*(dxsq(1)+dxsq(2)))
 end if
end function

end module
