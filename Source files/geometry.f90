module dam

contains

subroutine create_domain(p,np,length,x0,length_f,height,y0,height_f,dp,dp1)
use particles
!use TRANS__genmod
implicit none
type(particle),intent(out)::p(:)
integer,intent(out)::np
double precision,intent(in)::length,height,dp,length_f,height_f,dp1,x0,y0
integer::i,j,nn
integer::nx,nz,nx_f,nz_f
double precision::trans


nx=int(length/dp1)
nz=int(height/dp1)
nx_f=int(length_f/dp)
nz_f=int(height_f/dp)


p(:)%coord%x(1)=0.0!initialize all attributes
p(:)%coord%x(2)=0.0
nn=0
!4 layers of dummy fluid particles

!left boundary particles
 !creating 1st layer
do i=1,nz
  nn=nn+1
   p(nn)%coord%x(1)=0.0
   p(nn)%coord%x(2)=(i-1)*dp1
   p(nn)%id=2
 end do

 !creating 2nd layer

 do i=1,nz-1
    nn=nn+1
     p(nn)%coord%x(1)=trans(p(nn-nz+1)%coord%x(1),dp1)
     p(nn)%coord%x(2)=(i)*(dp1)
     p(nn)%id=2
 end do
 
   ! creating 3rd layer
    do i=1,nz-2
    nn=nn+1
     p(nn)%coord%x(1)=trans(p(nn-nz+2)%coord%x(1),dp1)
     p(nn)%coord%x(2)=trans(p(nn-nz+2)%coord%x(2),0.0d0)
     p(nn)%id=2
 end do

 !creating fourth layer
   do i=1,nz-3
    nn=nn+1
     p(nn)%coord%x(1)=trans(p(nn-nz+3)%coord%x(1),dp1)
     p(nn)%coord%x(2)=trans(p(nn-nz+3)%coord%x(2),0.0d0)
     p(nn)%id=2
 end do

 !creation of right boundary
 
 !creation of 1st layers
   do i=1,nz
    nn=nn+1
   p(nn)%coord%x(1)=length
   p(nn)%coord%x(2)=(i-1)*dp1
   p(nn)%id=2
 end do
  
  

    
  !creation of second layer
    do i=1,nz-1
     nn=nn+1
   p(nn)%coord%x(1)=trans(p(nn-nz+1)%coord%x(1),-dp1)
   p(nn)%coord%x(2)=p(nn-nz+1)%coord%x(2)
   p(nn)%id=2
   end do
 
    !creation of 3rd layer
    do i=1,nz-2
     nn=nn+1
   p(nn)%coord%x(1)=trans(p(nn-nz+2)%coord%x(1),-dp1)
   p(nn)%coord%x(2)=p(nn-nz+2)%coord%x(2)
   p(nn)%id=2
   end do

   !creation of fourth layer
    do i=1,nz-3
     nn=nn+1
   p(nn)%coord%x(1)=trans(p(nn-nz+3)%coord%x(1),-dp1)
   p(nn)%coord%x(2)=p(nn-nz+3)%coord%x(2)
   p(nn)%id=2
   end do

  !create bottom boundary
  
  !creating 1st layer
         
      do j=1,nx-1
       nn=nn+1
      p(nn)%coord%x(1)=j*dp1
      p(nn)%coord%x(2)=0.0d0
      p(nn)%id=2
    end do
     
!creating 2nd layer
   
    do j=1,nx-2
       nn=nn+1
    p(nn)%coord%x(1)=trans(p(nn-nx+1)%coord%x(1),0.0d0)
    p(nn)%coord%x(2)=trans(p(nn-nx+1)%coord%x(2),dp1)
    p(nn)%id=2
      end do
     


 !creation of 3rd layer
      

    do j=1,nx-3
       nn=nn+1
    p(nn)%coord%x(1)=trans(p(nn-nx+2)%coord%x(1),0.0d0)
    p(nn)%coord%x(2)=trans(p(nn-nx+2)%coord%x(2),dp1)
    p(nn)%id=2
      end do

  !creation of fourth layer
     do j=1,nx-4
       nn=nn+1
    p(nn)%coord%x(1)=trans(p(nn-nx+3)%coord%x(1),0.0d0)
    p(nn)%coord%x(2)=trans(p(nn-nx+3)%coord%x(2),dp1)
    p(nn)%id=2
      end do
     
!creation of fluid particles
 do i=1,nz_f
  do j=1,nx_f
  nn=nn+1
  p(nn)%coord%x(1)=x0+((j-1)*dp)+(4.0*dp1)
  p(nn)%coord%x(2)=y0+((i-1)*dp)+(4.0*dp1)
  p(nn)%id=1
  end do
    end do
   
   np=nn

 end subroutine
 

  end module  
  

 function trans(x,dx)
 double precision,intent(in)::x
 double precision,intent(in)::dx
 double precision::trans
 trans=x+dx
 end function

     


