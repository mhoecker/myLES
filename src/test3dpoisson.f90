program test3dpoisson
 use threedpoissonsolver
 implicit none
 integer, parameter :: Nx = 32
 integer, parameter :: Ny = 32
 integer, parameter :: Nz = 32
 real :: f(Nx,Ny,Nz)
 real :: x(Nx)
 real :: y(Ny)
 real :: z(Nz)
 real :: w(Nx,Ny,Nz)
 real, parameter :: dx(3) = (/1,2,3/)
 real, parameter :: pi = 4*atan(1.0)
 real,parameter :: tol=1e-6
 integer, parameter :: Nmax=1
 real :: dxsq(3)
 real :: errnorm
 integer :: i,j,k
 dxsq(:) = dx(:)**2
 do i=1,Nx
  x(i) = (i-.5*(Nx+1))*dx(1)
  do j=1,Ny
   y(j) = (j-.5*(Ny+1))*dx(2)
   do k=1,Nz
    z(k) = (k-.5*(Nz+1))*dx(3)

    f(i,j,k) = (2*y(j)/dx(2)-Ny)/Ny
    f(i,j,k) = f(i,j,k)*(2*y(j)/dx(2)+Ny)/Ny
    f(i,j,k) = f(i,j,k)*(2*x(i)/dx(1)-Nx)/Nx
    f(i,j,k) = f(i,j,k)*(2*x(i)/dx(1)+Nx)/Nx
    f(i,j,k) = f(i,j,k)*(2*z(k)/dx(3)+Nz)/Nz
    f(i,j,k) = f(i,j,k)*(2*z(k)/dx(3)-Nz)/Nz
    f(i,j,k) = 0
    if((i.eq.1).or.(i.eq.Nx).or.(j.eq.1).or.(j.eq.Ny).or.(k.eq.1).or.(k.eq.Nz)) then
     w(i,j,k) = z(k)/(Nz*dx(3))
    else
     w(i,j,k) = 1
    end if
   end do
  end do
 end do
 call threedpoisson(dxsq,f,w,errnorm,tol,Nmax)
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
