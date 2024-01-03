program density
implicit none

integer :: i,j,k,l,c1,ig,ig1,ig2,ig3,ok,dt,m,iz,t,n,num1,b,nbinx,nbiny,nbinz,natom
integer, parameter :: nstep10=2000!number of frames(200ns traj)
integer, parameter :: np=12744, nbin=25,nmyr=49,nhbr=281
integer, dimension (1000000) :: sl_no,id,sl_no1,id1
real, parameter :: dgrid=0.05, rcut=0.6,ts=100,xcenter=7.0,ycenter=7.0,zcenter=12.0
real:: lower_cut, upper_cut
character(len=*), parameter :: gro='(i5,2a5,i5,3f8.3)',last='(3f10.5)'
character (len=5), dimension (1000000) :: resi,atom,res1,atom1
real, dimension (1000000) :: x,y,z
real*8 :: dx, dy, dz, dr, dist, x_ref, y_ref, z_ref, x_com, y_com, z_com, r, cnt, junk,boxx,boxy,boxz
real*8::zmem,zmem_center, SAPI1z,SAPI2z,SAPI3z,SAPI4z,SAPI5z,SAPI6z,SAPI7z,SAPI8z,SAPI9z,SAPI10z
real*8::dist1,dist2,dist3,dist4,dist5,dist6,dist7,dist8,dist9,dist10,prob
real*8::dimer_com_x,dimer_com_y,dimer_com_z
integer,dimension(:,:),allocatable::den
CHARACTER(len=32) :: INPUT, OUTPUT
CHARACTER(len=5) :: lipid

CALL get_command_argument(1, INPUT)
CALL get_command_argument(2, OUTPUT)

write(*,*)'Enter a value of lower and upper z-coordinate cutoff, total number of atoms in traj file, lipid type(a5)'!cutoff to select lipid headgroups
read(*,*) lower_cut, upper_cut, natom, lipid

open(10, file=INPUT)
open(100, file=OUTPUT)

allocate(den(1000,1000),stat=ok)

do i=-100,500
do j=-100,500

den(i,j)=0

end do
end do

do k=1,nstep10
if(mod(k,50)==0)print*,"k",k
num1=0
zmem=0.0
zmem_center=0.0

dimer_com_x=0.0
dimer_com_y=0.0
dimer_com_z=0.0

 read(10,*)
 read(10,*)


 do i=1,natom

  read(10,gro)sl_no(i),resi(i),atom(i),id(i),x(i),y(i),z(i)

 end do

  read(10,last)boxx,boxy,boxz

  
 nbinx=boxx/dgrid
 nbiny=boxy/dgrid
 nbinz=boxz/dgrid

!!centering MA6 

 do i=8768,9048 !!HBR domain at MA3-MA3 interface

  dimer_com_x=dimer_com_x+x(i)
  dimer_com_y=dimer_com_y+y(i)
  dimer_com_z=dimer_com_z+z(i)
  
end do

  dimer_com_x=dimer_com_x/nhbr
  dimer_com_y=dimer_com_y/nhbr
  dimer_com_z=dimer_com_z/nhbr

 dx=xcenter-dimer_com_x
 dy=ycenter-dimer_com_y
 dz=zcenter-dimer_com_z

  write(200,*)"Title in water t=   ",k*ts," step= ",k
  write(200,*)natom

 do i=1,natom


x(i)=x(i)+dx
y(i)=y(i)+dy

if(x(i)<0.0)x(i)=x(i)+boxx
if(y(i)<0.0)y(i)=y(i)+boxy

if(x(i)>boxx)x(i)=x(i)-boxx
if(y(i)>boxy)y(i)=y(i)-boxy

if (resi(i)==lipid.and.z(i)>lower_cut.and.z(i)<upper_cut)then
nbinx=nint(x(i)/dgrid)+1 
nbiny=nint(y(i)/dgrid)+1         

den(nbinx,nbiny)=den(nbinx,nbiny)+1
       
end if

 end do

  write(200,last)boxx,boxy,boxz

end do

do i=-100,500
do j=-100,500

prob=real(den(i,j))/(nstep10*(dgrid**2)*(upper_cut-lower_cut))
write(100,*)i*dgrid,j*dgrid,prob

end do
end do


end program density

