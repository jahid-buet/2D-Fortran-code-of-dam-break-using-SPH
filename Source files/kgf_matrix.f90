module ptcl_matrix
contains
subroutine matrix(visco_type,kernal_type,length,height,p,ntotal,h) !compute inverse matrix necessary for kgf-sph method
use particles
use kernal
implicit none
double precision,intent(in)::h
integer,intent(in)::ntotal
type(particle),intent(inout)::p(:) !particle type
integer,intent(in)::kernal_type,visco_type
double precision,intent(in)::length,height
double precision::wij,dwdx,dwdy
double precision,allocatable::a11(:),a12(:),a13(:),a21(:),a22(:),a23(:),a31(:),a32(:),a33(:) !variables to store matix component for i particle
double precision::rij,rij_x,rij_z,det
integer::i,j,k,ii,lx,lz,ii_neigh
integer::ncx,ncz,nct
double precision::xmax,zmax!max domain size
integer,allocatable::ibox(:),ll(:),np(:)
integer::icell,kcell
integer::ndx(4),ndz(4) !no.of neighbour cells
double precision::kf
double precision::L(3,3)!kgf_matrix






if(kernal_type.eq.1)then
  kf=3.0d0
else
 kf=2.0d0
end if


xmax=length+kf*h !assuming square cell (rectangular cell is also possible,in this case(xmax.ne.zmax))
zmax=height+kf*h

ncx=int(xmax/(kf*h))
ncz=int(zmax/(kf*h))
nct=ncz*ncx






allocate(ibox(nct),ll(ntotal),np(nct))
allocate(a11(ntotal),a12(ntotal),a13(ntotal),a21(ntotal),a22(ntotal),a23(ntotal),a31(ntotal),a32(ntotal),a33(ntotal))

!initialize matrix component
  do i=1,ntotal
  a11(i)=0.0
  a12(i)=0.0
  a13(i)=0.0
  a21(i)=0.0
  a22(i)=0.0
  a23(i)=0.0
  a31(i)=0.0
  a32(i)=0.0
  a33(i)=0.0
  end do

  !initialize matrix

         do i=1,ntotal    
             p(i)%mat(1,1)=0.0d0
             p(i)%mat(1,2)=0.0d0
             p(i)%mat(1,3)=0.0d0
             p(i)%mat(2,1)=0.0d0
             p(i)%mat(2,2)=0.0d0
             p(i)%mat(2,3)=0.0d0
             p(i)%mat(3,1)=0.0d0
             p(i)%mat(3,2)=0.0d0
             p(i)%mat(3,3)=0.0d0
      
          end do

    !initialize cells
       
      do i=1,nct
        ibox(i)=0
        np(i)=0
     end do

  do i=1,ntotal
    ll(i)=0
   end do

     !sorting particles into cells

      do k=1,ntotal
      icell=int(p(k)%coord%x(1)/(kf*h))+1
      kcell=int(p(k)%coord%x(2)/(kf*h))+1
      ii=icell+(kcell-1)*ncx
      ll(k)=ibox(ii)
      ibox(ii)=k
      np(ii)=np(ii)+1 !no.of particles in each cell
      end do

       
       do lz=1,ncz
          do lx=1,ncx
           ii=lx+(lz-1)*ncx
       if(np(ii).gt.0)then!there are particles in the cell
         i=ibox(ii)

        do while(i.gt.0)
         j=ll(i)
           do while(j.gt.0) !loop in same cell
              !if(p(j)%id.eq.1)then !fluid particles only
                if(p(j)%id*p(i)%id.lt.4)then!fluid+boundary particles 
             rij_x=p(j)%coord%x(1)-p(i)%coord%x(1)
             rij_z=p(j)%coord%x(2)-p(i)%coord%x(2)
             rij=rij_x**2+rij_z**2
          
          if(rij.le.kf*kf*h*h)then
           call ComKer(kernal_type,rij,rij_x,rij_z,h,wij,dwdx,dwdy)
              
              a11(i)=a11(i)+(p(j)%mass/p(j)%dens)*wij
              a11(j)=a11(j)+(p(i)%mass/p(i)%dens)*wij             
              a12(i)=a12(i)+((p(j)%mass/p(j)%dens)*rij_x*wij)
              a12(j)=a12(j)-((p(i)%mass/p(i)%dens)*rij_x*wij)

              a13(i)=a13(i)+((p(j)%mass/p(j)%dens)*rij_z*wij)  
              a13(j)=a13(j)-((p(i)%mass/p(i)%dens)*rij_z*wij)
              
              a21(i)=a21(i)+((p(j)%mass/p(j)%dens)*rij_x*wij)
              a21(j)=a21(j)-((p(i)%mass/p(i)%dens)*rij_x*wij)  
              
              a22(i)=a22(i)+((p(j)%mass/p(j)%dens)*rij_x*rij_x*wij)
              a22(j)=a22(j)+((p(i)%mass/p(i)%dens)*rij_x*rij_x*wij)
              
              a23(i)=a23(i)+ ((p(j)%mass/p(j)%dens)*rij_x*rij_z*wij)
              a23(j)=a23(j)+ ((p(i)%mass/p(i)%dens)*rij_x*rij_z*wij)

              a31(i)=a31(i)+((p(j)%mass/p(j)%dens)*rij_z*wij)
              a31(j)=a31(j)-((p(i)%mass/p(i)%dens)*rij_z*wij)
              
              a32(i)=a32(i)+ ((p(j)%mass/p(j)%dens)*rij_x*rij_z*wij)               
              a32(j)=a32(j)+ ((p(i)%mass/p(i)%dens)*rij_x*rij_z*wij)
               
              a33(i)=a33(i)+ ((p(j)%mass/p(j)%dens)*rij_z*rij_z*wij)
              a33(j)=a33(j)+ ((p(i)%mass/p(i)%dens)*rij_z*rij_z*wij)
              
             end if
              end if
         
           j=ll(j)
    
            end do
           
       !now for neighbour cells

          ndx=(/1,0,1,-1/)
          ndz=(/0,1,1, 1/)
              
       do k=1,4

        icell=lx+ndx(k)
        kcell=lz+ndz(k)
        !PERIOD bc

       if(icell.lt.1) cycle

       if(icell.gt.ncx) cycle
        
       if(kcell.lt.1) cycle
        
       if(kcell.gt.ncz) cycle
        
       ii_neigh=icell+(kcell-1)*ncx

           if(np(ii_neigh).gt.0)then
           j=ibox(ii_neigh)
         do while(j.gt.0)
              !if(p(j)%id.eq.1)then
            if(p(j)%id*p(i)%id.lt.4)then
             rij_x=p(j)%coord%x(1)-p(i)%coord%x(1)
             rij_z=p(j)%coord%x(2)-p(i)%coord%x(2)
             rij=rij_x**2+rij_z**2

             if(rij.le.kf*kf*h*h)then
            call ComKer(kernal_type,rij,rij_x,rij_z,h,wij,dwdx,dwdy)
                  
              a11(i)=a11(i)+(p(j)%mass/p(j)%dens)*wij
              a11(j)=a11(j)+(p(i)%mass/p(i)%dens)*wij             
              a12(i)=a12(i)+((p(j)%mass/p(j)%dens)*rij_x*wij)
              a12(j)=a12(j)-((p(i)%mass/p(i)%dens)*rij_x*wij)

              a13(i)=a13(i)+((p(j)%mass/p(j)%dens)*rij_z*wij)  
              a13(j)=a13(j)-((p(i)%mass/p(i)%dens)*rij_z*wij)
              
              a21(i)=a21(i)+((p(j)%mass/p(j)%dens)*rij_x*wij)
              a21(j)=a21(j)-((p(i)%mass/p(i)%dens)*rij_x*wij)  
              
              a22(i)=a22(i)+((p(j)%mass/p(j)%dens)*rij_x*rij_x*wij)
              a22(j)=a22(j)+((p(i)%mass/p(i)%dens)*rij_x*rij_x*wij)
              
              a23(i)=a23(i)+ ((p(j)%mass/p(j)%dens)*rij_x*rij_z*wij)
              a23(j)=a23(j)+ ((p(i)%mass/p(i)%dens)*rij_x*rij_z*wij)

              a31(i)=a31(i)+((p(j)%mass/p(j)%dens)*rij_z*wij)
              a31(j)=a31(j)-((p(i)%mass/p(i)%dens)*rij_z*wij)
              
              a32(i)=a32(i)+ ((p(j)%mass/p(j)%dens)*rij_x*rij_z*wij)               
              a32(j)=a32(j)+ ((p(i)%mass/p(i)%dens)*rij_x*rij_z*wij)
               
              a33(i)=a33(i)+ ((p(j)%mass/p(j)%dens)*rij_z*rij_z*wij)
              a33(j)=a33(j)+ ((p(i)%mass/p(i)%dens)*rij_z*rij_z*wij)
              
           
                end if
                  end if
           j=ll(j)
            end do
                end if
              end do
                  
          i=ll(i)

          end do 

              end if

               end do
            
                      end do


 !computing inverse of matrix for fluid and boundary

      do i=1,ntotal
        !if(p(i)%id.eq.1)then !matrix for fluid particles only
        L(1,1)=a11(i)
        L(1,2)=a12(i)
        L(1,3)=a13(i)
        L(2,1)=a21(i)
        L(2,2)=a22(i)
        L(2,3)=a23(i)
        L(3,1)=a31(i)
        L(3,2)=a32(i)
        L(3,3)=a33(i)

      
        !matrix inversion
  
       det=(L(1,1)*L(2,2)*L(3,3)- L(1,1)*L(2,3)*L(3,2)&
              - L(1,2)*L(2,1)*L(3,3)+ L(1,2)*L(2,3)*L(3,1)&
              + L(1,3)*L(2,1)*L(3,2)- L(1,3)*L(2,2)*L(3,1))
              
       det=1.0/det       
           if(abs(det).gt.0.0d0)then !COMPUTING 3 BY 3 INVERSE MATRIX FOR FLUID PARTICLES
            
  
           p(i)%mat(1,1) = +det * (L(2,2)*L(3,3) - L(2,3)*L(3,2))
           p(i)%mat(2,1) = -det * (L(2,1)*L(3,2) - L(2,3)*L(3,3))
           p(i)%mat(3,1) = +det * (L(2,1)*L(3,2) - L(2,2)*L(3,1))
           p(i)%mat(1,2) = -det * (L(1,2)*L(3,3) - L(1,3)*L(3,2))
           p(i)%mat(2,2) = +det * (L(1,1)*L(3,3) - L(1,3)*L(3,1))
           p(i)%mat(3,2) = -det * (L(1,1)*L(3,2) - L(1,2)*L(3,1))
           p(i)%mat(1,3) = +det * (L(1,2)*L(2,3) - L(1,3)*L(2,2))
           p(i)%mat(2,3) = -det * (L(1,1)*L(2,3) - L(1,3)*L(2,1))
           p(i)%mat(3,3) = +det * (L(1,1)*L(2,2) - L(1,2)*L(2,1))

         else
             p(i)%mat(1,1)=1.0d0
             p(i)%mat(1,2)=0.0d0
             p(i)%mat(1,3)=0.0d0
             p(i)%mat(2,1)=0.0d0
             p(i)%mat(2,2)=1.0d0
             p(i)%mat(2,3)=0.0d0
             p(i)%mat(3,1)=0.0d0
             p(i)%mat(3,2)=0.0d0
             p(i)%mat(3,3)=1.0d0
          end if
            !end if
              end do

!$$$$$$   do i=1,ntotal
!$$$$$$      if(p(i)%id.gt.1)then!boundary particles
!$$$$$$              p(i)%mat(1,1)=1.0d0
!$$$$$$              p(i)%mat(1,2)=0.0d0
!$$$$$$              p(i)%mat(1,3)=0.0d0
!$$$$$$              p(i)%mat(2,1)=0.0d0
!$$$$$$              p(i)%mat(2,2)=1.0d0
!$$$$$$              p(i)%mat(2,3)=0.0d0
!$$$$$$              p(i)%mat(3,1)=0.0d0
!$$$$$$              p(i)%mat(3,2)=0.0d0
!$$$$$$              p(i)%mat(3,3)=1.0d0
!$$$$$$        end if
!$$$$$$    end do







                      

         end subroutine
     end module



























!initialize matrix component
  !ai11=0.0
  !ai12=0.0
  !ai13=0.0
  !ai21=0.0
  !ai22=0.0
  !ai23=0.0
  !ai31=0.0
  !ai32=0.0
  !ai33=0.0


    !do i=1,ntotal

      !do j=1,ntotal
          !if(i.ne.j)then

     !rij_x=p(j)%coord%x(1)-p(i)%coord%x(1)
     !rij_z=p(j)%coord%x(3)-p(i)%coord%x(3)
     !rij=sqrt(rij_x**2+rij_z**2)

        !if(rij.le.2.0*h)then
           !rij_x=p(j)%coord%x(1)-p(i)%coord%x(1)
           !rij_z=p(j)%coord%x(3)-p(i)%coord%x(3)

       !matrix formulation for i particle

      !ai11=ai11+(p(j)%mass/p(j)%dens)*wij
      !ai12=ai12+((p(j)%mass/p(j)%dens)*rij_x*wij)
      !ai13=ai13+((p(j)%mass/p(j)%dens)*rij_z*wij)  
      !ai21=ai21+((p(j)%mass/p(j)%dens)*rij_x*wij)
      !ai22=ai22+((p(j)%mass/p(j)%dens)*rij_x*rij_x*wij)
     ! ai23=ai23+ ((p(j)%mass/p(j)%dens)*rij_x*rij_z*wij)
      !ai31=ai31+((p(j)%mass/p(j)%dens)*rij_z*wij)
      !ai32=ai32+ ((p(j)%mass/p(j)%dens)*rij_x*rij_z*wij)
      !ai33=ai33+ ((p(j)%mass/p(j)%dens)*rij_z*rij_z*wij)
   !end if
     !end if
        !end do


  !p(i)%mat(1,1)=ai11
  !p(i)%mat(1,2)=ai12
  !p(i)%mat(1,3)=ai13
  !p(i)%mat(2,1)=ai21
  !p(i)%mat(2,2)=ai22
  !p(i)%mat(2,3)=ai23
  !p(i)%mat(3,1)=ai31
  !p(i)%mat(3,2)=ai32
  !p(i)%mat(3,3)=ai33
   
  !matrix inversion
  
!det= 1.0/(p(i)%mat(1,1)*p(i)%mat(2,2)*p(i)%mat(3,3) - p(i)%mat(1,1)*p(i)%mat(2,3)*p(i)%mat(3,2)&
             ! - p(i)%mat(1,2)*p(i)%mat(2,1)*p(i)%mat(3,3) + p(i)%mat(1,2)*p(i)%mat(2,3)*p(i)%mat(3,1)&
             ! + p(i)%mat(1,3)*p(i)%mat(2,1)*p(i)%mat(3,2) - p(i)%mat(1,3)*p(i)%mat(2,2)*p(i)%mat(3,1))


 
  !calculate of inverse

    ! p(i)%mat(1,1) = +det * (p(i)%mat(2,2)*p(i)%mat(3,3) - p(i)%mat(2,3)*p(i)%mat(3,2))
    ! p(i)%mat(2,1) = -det * (p(i)%mat(2,1)*p(i)%mat(3,3) - p(i)%mat(2,3)*p(i)%mat(3,1))
     !p(i)%mat(3,1) = +det * (p(i)%mat(2,1)*p(i)%mat(3,2) - p(i)%mat(2,2)*p(i)%mat(3,1))
     !p(i)%mat(1,2) = -det * (p(i)%mat(1,2)*p(i)%mat(3,3) - p(i)%mat(1,3)*p(i)%mat(3,2))
     !p(i)%mat(2,2) = +det * (p(i)%mat(1,1)*p(i)%mat(3,3) - p(i)%mat(1,3)*p(i)%mat(3,1))
    ! p(i)%mat(3,2) = -det * (p(i)%mat(1,1)*p(i)%mat(3,2) - p(i)%mat(1,2)*p(i)%mat(3,1))
     !p(i)%mat(1,3) = +det * (p(i)%mat(1,2)*p(i)%mat(2,3) - p(i)%mat(1,3)*p(i)%mat(2,2))
     !p(i)%mat(2,3) = -det * (p(i)%mat(1,1)*p(i)%mat(2,3) - p(i)%mat(1,3)*p(i)%mat(2,1))
     !p(i)%mat(3,3) = +det * (p(i)%mat(1,1)*p(i)%mat(2,2) - p(i)%mat(1,2)*p(i)%mat(2,1))
 
   !end do
    !end subroutine
     !end module

   