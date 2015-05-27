program test2dpoisson
 use twodpoissonsolver
 implicit none
 integer, parameter :: Nx = 32
 integer, parameter :: Ny = 32
 real :: f(Nx,Ny)
 real :: x(Nx)
 real :: y(Ny)
 real :: z(Nx,Ny)
 real, parameter :: dx = 1
 real, parameter :: dy = 2
 real,parameter :: tol=1e-5
 integer, parameter :: Nmax=Nx*Ny*floor(sqrt(1.0*Nx*Ny))
 real :: dxsq(2)
 real :: err
 integer :: i,j
 dxsq(1) = dx*dx
 dxsq(2) = dy*dy
 do i=1,Nx
  x(i) = (i-.5*(Nx+1))*dx
  do j=1,Ny
   y(j) = (j-.5*(Ny+1))*dy
   f(i,j) = (2*y(j)/dy-Ny)/Ny
   f(i,j) = f(i,j)*(2*y(j)/dy+Ny)/Ny
   f(i,j) = f(i,j)*(2*x(i)/dx-Nx)/Nx
   f(i,j) = f(i,j)*(2*x(i)/dx+Nx)/Nx
   if((i.eq.1).or.(i.eq.Nx).or.(j.eq.1).or.(j.eq.Ny)) then
    z(i,j) = 0
   else
    z(i,j) = 1
   end if
  end do
 end do
 call twodpoisson(dxsq,f,z,err,tol,Nmax)
 do i=1,Nx
  do j=1,Ny
   print *,x(i),y(j),z(i,j),f(i,j),laplace(z,dxsq,i,j),jacobi(f,z,dxsq,i,j)
  end do
  print *
 end do
end program
