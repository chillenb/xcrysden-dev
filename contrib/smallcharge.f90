!------------------------------------------------------------------------------
! Small Readcharge program
! Conversion from .xsf file (XCrySDen) to x,y,z,rho(x,y,z) format file.
! Arnaud Ponchel- Laboratoire d'Etudes des Microstructures - ONERA - 01/2000 -
!------------------------------------------------------------------------------

       Program readcharge

! nx,ny,nz: number of charge densities points in the 3 directions
! a,b,c: lattice parameters
! x,y,z: cartesian coordinates of a point in the 3D grid
! ix,iy,iz: incrementations of x,y,z

! Must entry the dimension of the array

       integer :: l,n,u,nx,ny,nz
       real,dimension(524800) :: t
       real :: a,b,c,x,y,z,ix,iy,iz

       open(96,file='RHO.inp',status='unknown')
       open(97,file='RHO.out',status='new')

! Definition of parameters

       u=0
       n=1
       a=7.996/0.529
       b=7.996/0.529
       c=8.158/0.529
       nx=80
       ny=80
       nz=82
       l=nx*ny*nz-5
       ix=a/nx
       iy=b/ny
       iz=c/nz

! Reading of the rho values in the array t(i)
! From the case.xsf file (only 3D-grid values)
! (6 columns file)

       do while (n<l)
       read (96,*) t(n),t(n+1),t(n+2),t(n+3),t(n+4),t(n+5)
       n=n+6
       end do

! Writing of the rho values in the new format x,y,z,rho(x,y,z)

       do x=0,a,ix
       do y=0,b,iy
       do z=0,c,iz

       u=u+1
       write (97,*) x,y,z,t(u)

       end do
       end do
       end do

       close(96)
       close(97)
       end program readcharge
