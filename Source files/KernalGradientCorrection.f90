module kernal_gradient_correction
contains
subroutine KGC(p,ntotal,h,dp,length,height,kernal_type)
use particles
use kernal
implicit none
double precision,intent(in)::h
integer,intent(in)::ntotal
type(particle),intent(inout)::p(:) !particle type
integer,intent(in)::kernal_type
double precision,intent(in)::length,height
double precision,intent(in)::dp
double precision::wij,dwdx,dwdy
double precision,allocatable::a11(:),a12(:),a21(:),a22(:)!local variables to store matix component for i particle
double precision::rij,rij_x,rij_y,det
double precision::L(2,2) !MATRIX COMPONENT
integer,save::i,j,k,ii,lx,ly,ii_neigh
integer::ncx,ncy,nct
double precision::xmax,ymax!max domain size
integer,allocatable::ibox(:),ll(:),np(:)
integer::icell,kcell,jj,j1,j2
integer::ndx(4),ndz(4) !no.of neighbour cells
double precision::kf


if(kernal_type.eq.1)then
  kf=3.0d0
else
 kf=2.0d0
end if

   xmax=length+kf*h !maximum domain extension in x and y direction for cell-division
   ymax=height+kf*h



ncx=nint(xmax/(kf*h))
ncy=nint(ymax/(kf*h))
nct=ncy*ncx

 allocate(ibox(nct),np(nct),ll(ntotal))
 allocate(a11(ntotal),a12(ntotal),a21(ntotal),a22(ntotal))

  ibox=0
  do i=1,nct
   np(i)=0
  end do
  
  do i=1,ntotal
    ll(i)=0
   end do 

!initialize matrix component
 do i=1,ntotal
  a11(i)=0.0
  a12(i)=0.0
  a21(i)=0.0
  a22(i)=0.0
  end do

 !initialize matrix as identity matrix
do i=1,ntotal
 p(i)%mat(1,1)=1.0d0
 p(i)%mat(1,2)=0.0d0
 p(i)%mat(2,1)=0.0d0
 p(i)%mat(2,2)=1.0d0
end do

!$$$$$$ !sorting particles into each cell
!$$$$$$      do k=1,ntotal
!$$$$$$       icell=int(p(k)%coord%x(1)/(3.0*h))+1
!$$$$$$       kcell=int(p(k)%coord%x(2)/(3.0*h))+1
!$$$$$$       ii=icell+(kcell-1)*ncx
!$$$$$$       np(ii)=np(ii)+1 
!$$$$$$       ibox(ii,np(ii))=k      
!$$$$$$       end do
 





  !sorting particles into cells
      do k=1,ntotal      
      icell=int(p(k)%coord%x(1)/(kf*h))+1
      kcell=int(p(k)%coord%x(2)/(kf*h))+1              
      ii=icell+(kcell-1)*ncx
      ll(k)=ibox(ii)
      ibox(ii)=k
      np(ii)=np(ii)+1 !no.of particles in each cell      
      end do




       do ly=1,ncy
          do lx=1,ncx
           ii=lx+(ly-1)*ncx
       if(np(ii).gt.0)then!there are particles in the cell
         i=ibox(ii)
                         
         do while(i.gt.0)
        !if(p(i)%id.eq.1)then!fluid particles only
           j=ll(i)
           do while(j.gt.0) !loop in same cell
             if(p(j)%id.eq.1)then
            
       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x**2+rij_y**2

        if(rij.le.kf*kf*h*h)then
         

            call ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)
            !dwdx=wij*(-rij_x)
            !dwdy=wij*(-rij_y)

            a11(i)=a11(i)+(p(j)%mass/p(j)%dens)*(-rij_x)*dwdx
            a11(j)=a11(j)+(p(i)%mass/p(i)%dens)*(-rij_x)*dwdx
            a12(i)=a12(i)+(p(j)%mass/p(j)%dens)*(-rij_y)*dwdx
            a12(j)=a12(j)+(p(i)%mass/p(i)%dens)*(-rij_y)*dwdx
            a21(i)=a21(i)+(p(j)%mass/p(j)%dens)*(-rij_x)*dwdy
            a21(j)=a21(j)+(p(i)%mass/p(i)%dens)*(-rij_x)*dwdy
            a22(i)=a22(i)+(p(j)%mass/p(j)%dens)*(-rij_y)*dwdy
            a22(j)=a22(j)+(p(i)%mass/p(i)%dens)*(-rij_y)*dwdy
              end if

            end if
           j=ll(j)
         end do

       !now for neighbour cells

          ndx=(/1,0,1,-1/)
          ndz=(/0,1,1, 1/)

        do k=1,4
        icell=lx+ndx(k)
        kcell=ly+ndz(k)

        !PERIOD bc

       

        if(icell.lt.1) cycle
        
        if(icell.gt.ncx) cycle

        if(kcell.lt.1) cycle
      
        if(kcell.gt.ncy) cycle

       ii_neigh=icell+(kcell-1)*ncx
 
        if(np(ii_neigh).gt.0)then
           j=ibox(ii_neigh)
         do while(j.gt.0)
             if(p(j)%id.eq.1)then

       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x**2+rij_y**2

        if(rij.le.kf*kf*h*h)then
           
          call ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)
     
            a11(i)=a11(i)+(p(j)%mass/p(j)%dens)*(-rij_x)*dwdx
            a11(j)=a11(j)+(p(i)%mass/p(i)%dens)*(-rij_x)*dwdx
            a12(i)=a12(i)+(p(j)%mass/p(j)%dens)*(-rij_y)*dwdx
            a12(j)=a12(j)+(p(i)%mass/p(i)%dens)*(-rij_y)*dwdx
            a21(i)=a21(i)+(p(j)%mass/p(j)%dens)*(-rij_x)*dwdy
            a21(j)=a21(j)+(p(i)%mass/p(i)%dens)*(-rij_x)*dwdy
            a22(i)=a22(i)+(p(j)%mass/p(j)%dens)*(-rij_y)*dwdy
            a22(j)=a22(j)+(p(i)%mass/p(i)%dens)*(-rij_y)*dwdy
     
            end if
     
          end if

       j=ll(j)
   
     end do

       end if

         end do
     
           !end if

          i=ll(i)
          end do 
                      
              end if
               end do
                     end do

    !computing inverse of matrix for fluid and boundary

      do i=1,ntotal
        if(p(i)%id.eq.1)then !matrix for fluid paarticles only
         L(1,1)=a11(i)
         L(1,2)=a12(i)
         L(2,1)=a21(i)
         L(2,2)=a22(i)

      
        det=L(1,1)*L(2,2)-L(1,2)*L(2,1)
           if(abs(det).gt.0.0d0)then
  
          p(i)%mat(1,1)=(1.0/det)*L(2,2)
          p(i)%mat(1,2)=(1.0/det)*(-L(1,2))
          p(i)%mat(2,1)=(1.0/det)*(-L(2,1))
          p(i)%mat(2,2)=(1.0/det)*L(1,1)

         else
             p(i)%mat(1,1)=1.0
             p(i)%mat(1,2)=0.0d0
             p(i)%mat(2,1)=0.0d0
             p(i)%mat(2,2)=1.0d0
          end if
            end if
              end do

  do i=1,ntotal
     if(p(i)%id.gt.1)then!boundary particles
             p(i)%mat(1,1)=1.0
             p(i)%mat(1,2)=0.0d0
             p(i)%mat(2,1)=0.0d0
             p(i)%mat(2,2)=1.0d0
       end if
   end do
       








!$$$$$$         do ly=1,ncy
!$$$$$$           do lx=1,ncx
!$$$$$$            ii=lx+(ly-1)*ncx
!$$$$$$        if(np(ii).gt.0)then!there are particles in the cell
!$$$$$$           do j1=1,np(ii)
!$$$$$$         i=ibox(ii,j1)       
!$$$$$$       if(p(i)%id.eq.1)then !fluid particles only
!$$$$$$                                                               !id=1 is fluid partilces
!$$$$$$                                                                                                                            !id=2 boundary particles
!$$$$$$          do j2=1,np(ii)                      
!$$$$$$          j=ibox(ii,j2)
!$$$$$$         if(i.ne.j.and.p(j)%id.eq.1)then !fluid particles only                             
!$$$$$$        rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
!$$$$$$        rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
!$$$$$$        rij=rij_x*rij_x +rij_y*rij_y
!$$$$$$         if(sqrt(rij).le.3.0*h)then
!$$$$$$                
!$$$$$$        call ComKer(rij,rij_x,rij_y,h,wij,dwdx,dwdy) 
!$$$$$$               
!$$$$$$         a11(i)=a11(i)+(p(j)%mass/p(j)%dens)*dwdx*(-rij_x)
!$$$$$$         a12(i)=a12(i)+(p(j)%mass/p(j)%dens)*dwdx*(-rij_y)
!$$$$$$         a21(i)=a21(i)+(p(j)%mass/p(j)%dens)*dwdy*(-rij_x)
!$$$$$$         a22(i)=a22(i)+(p(j)%mass/p(j)%dens)*dwdy*(-rij_y)
!$$$$$$      
!$$$$$$             end if
!$$$$$$               end if
!$$$$$$              
!$$$$$$               end do
!$$$$$$ 
!$$$$$$   !searcing over neighbour cells
!$$$$$$         
!$$$$$$           ndx=(/1,0,0,1,-1,-1,1,-1/)
!$$$$$$           ndz=(/0,1,-1,1,1,0,-1,-1/)
!$$$$$$               
!$$$$$$        do k=1,8
!$$$$$$ 
!$$$$$$         icell=lx+ndx(k)
!$$$$$$         kcell=ly+ndz(k)
!$$$$$$         !PERIOD bc
!$$$$$$ 
!$$$$$$        if(icell.lt.1) cycle
!$$$$$$ 
!$$$$$$        if(icell.gt.ncx) cycle
!$$$$$$         
!$$$$$$        if(kcell.lt.1) cycle
!$$$$$$         
!$$$$$$        if(kcell.gt.ncy) cycle
!$$$$$$         
!$$$$$$        ii_neigh=icell+(kcell-1)*ncx
!$$$$$$             
!$$$$$$         if(np(ii_neigh).gt.0)then
!$$$$$$           do jj=1,np(ii_neigh)
!$$$$$$              j=ibox(ii_neigh,jj)
!$$$$$$              
!$$$$$$            if(p(j)%id.eq.1)then !fluid particles only
!$$$$$$            rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
!$$$$$$            rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
!$$$$$$            rij=rij_x**2+rij_y**2
!$$$$$$       
!$$$$$$        if(sqrt(rij).le.3.0*h)then          
!$$$$$$       call ComKer(rij,rij_x,rij_y,h,wij,dwdx,dwdy)
!$$$$$$ 
!$$$$$$         a11(i)=a11(i)+(p(j)%mass/p(j)%dens)*dwdx*(-rij_x)
!$$$$$$         a12(i)=a12(i)+(p(j)%mass/p(j)%dens)*dwdx*(-rij_y)
!$$$$$$         a21(i)=a21(i)+(p(j)%mass/p(j)%dens)*dwdy*(-rij_x)
!$$$$$$         a22(i)=a22(i)+(p(j)%mass/p(j)%dens)*dwdy*(-rij_y)
!$$$$$$ 
!$$$$$$              
!$$$$$$           end if 
!$$$$$$              end if           
!$$$$$$                 end do
!$$$$$$                    end if
!$$$$$$                       end do
!$$$$$$ 
!$$$$$$ 
!$$$$$$                  end if
!$$$$$$                    
!$$$$$$          end do       
!$$$$$$               end if         
!$$$$$$                  end do
!$$$$$$                        end do
!$$$$$$ 
!$$$$$$     !computing inverse of matrix
!$$$$$$  
!$$$$$$       do i=1,ntotal
!$$$$$$          if(p(i)%id.eq.1)then !matrix for fluid paarticles only
!$$$$$$           L(1,1)=a11(i)
!$$$$$$           L(1,2)=a12(i)
!$$$$$$           L(2,1)=a21(i)
!$$$$$$           L(2,2)=a22(i)
!$$$$$$  
!$$$$$$        
!$$$$$$          det=L(1,1)*L(2,2)-L(1,2)*L(2,1)
!$$$$$$             if(det.ne.0.0d0)then
!$$$$$$    
!$$$$$$            p(i)%mat(1,1)=(1.0/det)*L(2,2)
!$$$$$$            p(i)%mat(1,2)=(1.0/det)*(-L(1,2))
!$$$$$$            p(i)%mat(2,1)=(1.0/det)*(-L(2,1))
!$$$$$$            p(i)%mat(2,2)=(1.0/det)*L(1,1)
!$$$$$$  
!$$$$$$           else
!$$$$$$               p(i)%mat(1,1)=1.0
!$$$$$$               p(i)%mat(1,2)=0.0d0
!$$$$$$               p(i)%mat(2,1)=0.0d0
!$$$$$$               p(i)%mat(2,2)=1.0d0
!$$$$$$            end if
!$$$$$$              end if
!$$$$$$                end do










           end subroutine
                end module