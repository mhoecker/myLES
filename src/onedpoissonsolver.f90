module onedpoissonsolver
contains
! Assumes end points are fixed values
! solves d^2y/dx^2 = f
! uses three point stencil to solve the DE
subroutine onedpoisson(dxsq,f,y,errnorm,tol,Nmax)
 implicit none
 integer, intent(in) :: Nmax
 real, intent(in) :: tol
 real, intent(out) :: errnorm
 real, intent(in) :: dxsq
 real, intent(inout) :: y(:)
 real, intent(in) :: f(:)
 integer :: N
 
 errnorm = 2*tol
 N=0
 do while((errnorm.gt.tol).and.(N.lt.Nmax))
  N = N+1
  call onedjacobi(dxsq,f,y,errnorm)
 end do
end subroutine

subroutine onedjacobi(dxsq,f,y,errnorm)
 real, intent(out) :: errnorm
 real, intent(in) :: dxsq
 real, intent(inout) :: y(:)
 real, intent(in) :: f(:)
 integer :: nx
 integer :: i,j
 
 if(size(y).ne.size(f)) then
  print *, "onedjacobi: ERROR Forcing and value unequal size"
 end if
 nx = size(y)
 errnorm = 0.0
! Edge cases try 1
 i=2
 j=nx-1
 y(i) = threejacobi(f,y,dxsq,i)
 y(j) = threejacobi(f,y,dxsq,j)
 errnorm = errnorm+abs(f(i)-threepointder2(y,dxsq,i))
 errnorm = errnorm+abs(f(j)-threepointder2(y,dxsq,j))
! use five point stencil for interior
 do i=3,floor(nx/2.0)
   j = nx+1-i
   y(i) = fivejacobi(f,y,dxsq,i)
   y(j) = fivejacobi(f,y,dxsq,j)
   errnorm = errnorm+abs(f(i)-fivepointder2(y,dxsq,i))
   errnorm = errnorm+abs(f(j)-fivepointder2(y,dxsq,j))
 end do
 if(modulo(nx,2).eq.1) then
  i=ceiling(nx/2.0)
  y(i) = fivejacobi(f,y,dxsq,i)
  errnorm = errnorm+abs(f(i)-fivepointder2(y,dxsq,i))
 end if

end subroutine

real function threepointder2(y,dxsq,i)
 real, intent(in) :: y(:)
 real, intent(in) :: dxsq
 integer, intent(in) :: i
 threepointder2 = (y(i+1)+y(i-1)-2*y(i))/(dxsq)
end function

real function fivepointder2(y,dxsq,i)
 real, intent(in) :: y(:)
 real, intent(in) :: dxsq
 integer, intent(in) :: i
 fivepointder2 = -30*y(i)+16*(y(i+1)+y(i-1))-y(i+2)-y(i-2)
 fivepointder2 = fivepointder2/(6*dxsq)
end function

real function fivejacobi(f,y,dxsq,i)
 real, intent(in) :: f(:)
 real, intent(in) :: y(:)
 real, intent(in) :: dxsq
 integer, intent(in) :: i
 fivejacobi = (-y(i-2)-y(i+2)+16*(y(i-1)+y(i+1))-6*f(i)*dxsq)/30
end function

real function threejacobi(f,y,dxsq,i)
 real, intent(in) :: f(:)
 real, intent(in) :: y(:)
 real, intent(in) :: dxsq
 integer, intent(in) :: i
 threejacobi = (y(i-1)+y(i+1)-f(i)*dxsq)/2
end function


end module
