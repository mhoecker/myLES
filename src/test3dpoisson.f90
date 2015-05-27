program test2dpoisson
 use threedpoissonsolver
 implicit none
 integer, parameter :: Nx = 10
 integer, parameter :: Ny = 10
 integer, parameter :: Nz = 10
 real :: f(Nx,Ny,Nz)
 real :: x(Nx)
 real :: y(Ny)
 real :: z(Nz)
 real :: w(Nx,Ny,Nz)
 real, parameter :: dx = 1.0
 real, parameter :: dy = 1.0
 real, parameter :: dz = 1.0
 real, parameter :: pi = 4*atan(1.0)
 real,parameter :: tol=1e-5
 integer, parameter :: Nmax=(Nx+Ny+Nz)**3
 real :: dxsq(3)
 real :: err
 integer :: i,j,k
 dxsq(1) = dx*dx
 dxsq(2) = dy*dy
 dxsq(3) = dz*dz
 do i=1,Nx
  x(i) = (i-.5*(Nx+1))*dx
  do j=1,Ny
   y(j) = (j-.5*(Ny+1))*dy
   do k=1,Nz
    z(k) = (k-.5*(Nz+1))*dz
    
    f(i,j,k) = (2*y(j)/dy-Ny)/Ny
    f(i,j,k) = f(i,j,k)*(2*y(j)/dy+Ny)/Ny
    f(i,j,k) = f(i,j,k)*(2*x(i)/dx-Nx)/Nx
    f(i,j,k) = f(i,j,k)*(2*x(i)/dx+Nx)/Nx
    f(i,j,k) = f(i,j,k)*(2*z(k)/dz+Nz)/Nz
    f(i,j,k) = f(i,j,k)*(2*z(k)/dz-Nz)/Nz
    f(i,j,k) = f(i,j,k)
    if((i.eq.1).or.(i.eq.Nx).or.(j.eq.1).or.(j.eq.Ny).or.(k.eq.1).or.(k.eq.Nz)) then
     w(i,j,k) = z(k)/(Nz*dz)+cos(2*pi*x(i)/(Nx*dx))+cos(4*pi*y(j)/(Ny*dy))
    else
     w(i,j,k) = 1
    end if
   end do
  end do
 end do
 call threedpoisson(dxsq,f,w,err,tol,Nmax)
 do i=1,Nx
  do j=1,Ny
   do k=1,Nz
    print *,x(i),y(j),z(k),w(i,j,k),f(i,j,k),laplace(w,dxsq,i,j,k),jacobi(f,w,dxsq,i,j,k)
   end do
   print *
  end do
  print *
 end do
end program
