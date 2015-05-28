module threedderivatives

contains
! d /dx derrivative
real function ddx(w,dx,i,j,k)
 real, intent(in) :: w(:,:,:)
 real, intent(in) :: dx(3)
 integer, intent(in) :: i,j,k
 integer :: nx
 nx = size(w,1)
! use 1 sided 3 point stencil for border points
 if(i.eq.1) then
  ddx = (4*w(i+1,j,k)-w(i+2,j,k)-3*w(i,j,k))/(2*dx(1))
 else if(i.eq.nx) then
  ddx = (4*w(i-1,j,k)-w(i-2,j,k)-3*w(i,j,k))/(-2*dx(1))
! use 3 point centered stencil near edge
 else if((i.eq.2).or.(i.eq.nx-1)) then
  ddx = (w(i+1,j,k)-w(i-1,j,k))/(2*dx(1))
! use 5 point centered stencil for interior points
 else
  ddx = (8*(w(i+1,j,k)-w(i-1,j,k))-w(i+2,j,k)+w(i-2,j,k))/(12*dx(1))
 end if
end function

! d /dy derrivative
real function ddy(w,dx,i,j,k)
 real, intent(in) :: w(:,:,:)
 real, intent(in) :: dx(3)
 integer, intent(in) :: i,j,k
 integer :: ny
 ny = size(w,2)
! use 1 sided 3 point stencil for border points
 if(j.eq.1) then
  ddy = (4*w(i,j+1,k)-w(i,j+2,k)-3*w(i,j,k))/(2*dx(2))
 else if(j.eq.ny) then
  ddy = (4*w(i,j-1,k)-w(i,j-2,k)-3*w(i,j,k))/(-2*dx(2))
! use 3 point centered stencil near edge
 else if((j.eq.2).or.(j.eq.ny-1)) then
  ddy = (w(i,j+1,k)-w(i,j-1,k))/(2*dx(2))
! use 5 point centered stencil for interior points
 else
  ddy = (8*(w(i,j+1,k)-w(i,j-1,k))-w(i,j+2,k)+w(i,j-2,k))/(12*dx(2))
 end if
end function

! d /dz derrivative
real function ddz(w,dx,i,j,k)
 real, intent(in) :: w(:,:,:)
 real, intent(in) :: dx(3)
 integer, intent(in) :: i,j,k
 integer :: nz
 nz = size(w,3)
! use 1 sided 3 point stencil for border points
 if(k.eq.1) then
  ddz = (4*w(i,j,k+1)-w(i,j,k+2)-3*w(i,j,k))/(2*dx(3))
 else if(k.eq.nz) then
  ddz = (4*w(i,j,k-1)-w(i,j,k-2)-3*w(i,j,k))/(-2*dx(3))
! use 3 point centered stencil near edge
 else if((k.eq.2).or.(k.eq.nz-1)) then
  ddz = (w(i,j,k+1)-w(i,j,k-1))/(2*dx(3))
! use 5 point centered stencil for interior points
 else
  ddz = (8*(w(i,j,k+1)-w(i,j,k-1))-w(i,j,k+2)+w(i,j,k-2))/(12*dx(3))
 end if
end function

! d^2 /dx^2 derrivative
real function ddxsq(w,dxsq,i,j,k)
 real, intent(in) :: w(:,:,:)
 real, intent(in) :: dxsq(3)
 integer, intent(in) :: i,j,k
 integer :: nx
 nx = size(w,1)
! use 1 sided 3 point stencil for border points
 if(i.eq.1) then
  ddxsq = (w(i+2,j,k)+w(i,j,k)-2*w(i+1,j,k))/dxsq(1)
 else if(i.eq.nx) then
  ddxsq = (w(i,j,k)+w(i-2,j,k)-2*w(i-1,j,k))/dxsq(1)
! use 3 point centered stencil near edge
 else if((i.eq.2).or.(i.eq.nx-1)) then
  ddxsq = (w(i+1,j,k)+w(i-1,j,k)-2*w(i,j,k))/dxsq(1)
! use 5 point centered stencil for interior points
 else
  ddxsq = (-30*w(i,j,k)+16*(w(i+1,j,k)+w(i-1,j,k))-w(i+2,j,k)-w(i-2,j,k))/(6*dxsq(1))
 end if
end function

! d^2 /dy^2 derrivative
real function ddysq(w,dxsq,i,j,k)
 real, intent(in) :: w(:,:,:)
 real, intent(in) :: dxsq(3)
 integer, intent(in) :: i,j,k
 integer :: ny
 ny = size(w,2)
! y derrivative
! use 1 sided 3 point stencil for border points
 if(j.eq.1) then
  ddysq = (w(i,j+2,k)+w(i,j,k)-2*w(i,j+1,k))/dxsq(2)
 else if(j.eq.ny) then
  ddysq = (w(i,j,k)+w(i,j-2,k)-2*w(i,j-1,k))/dxsq(2)
! use 3 point centered stencil near edge
 else if((j.eq.2).or.(j.eq.ny-1)) then
  ddysq = (w(i,j+1,k)+w(i,j-1,k)-2*w(i,j,k))/dxsq(2)
! use 5 point centered stencil for interior points
 else
  ddysq = (-30*w(i,j,k)+16*(w(i,j+1,k)+w(i,j-1,k))-w(i,j+2,k)-w(i,j-2,k))/(6*dxsq(2))
 end if
end function

! d^2 /dz^2 derrivative
real function ddzsq(w,dxsq,i,j,k)
 real, intent(in) :: w(:,:,:)
 real, intent(in) :: dxsq(3)
 integer, intent(in) :: i,j,k
 integer :: nz
 nz = size(w,3)
! z derrivative
! use 1 sided 3 point stencil for border points
 if(k.eq.1) then
  ddzsq = (w(i,j,k+2)+w(i,j,k)-2*w(i,j,k+1))/dxsq(3)
 else if(k.eq.nz) then
  ddzsq = (w(i,j,k)+w(i,j,k-2)-2*w(i,j,k-1))/dxsq(3)
! use 3 point centered stencil near edge
 else if((k.eq.2).or.(k.eq.nz-1)) then
  ddzsq = (w(i,j,k+1)+w(i,j,k-1)-2*w(i,j,k))/dxsq(3)
! use 5 point centered stencil for interior points
 else
  ddzsq = (-30*w(i,j,k)+16*(w(i,j,k+1)+w(i,j,k-1))-w(i,j,k+2)-w(i,j,k-2))/(6*dxsq(3))
 end if
end function

! Summ d^2w/dx^2, d^2w/dy^2, and d^2w/dz^2 derrivatives
real function laplace(w,dxsq,i,j,k)
 real, intent(in) :: w(:,:,:)
 real, intent(in) :: dxsq(3)
 integer, intent(in) :: i,j,k
 integer :: nx,ny,nz
 nz = size(w,3)
 laplace = ddxsq(w,dxsq,i,j,k)
 laplace = laplace+ ddysq(w,dxsq,i,j,k)
 laplace = laplace+ ddzsq(w,dxsq,i,j,k)
end function


end module
