module Correction
contains
subroutine xsph_corr(p,ntotal,h,corrvx,corrvy)
use particles
use kernal
implicit none
type(particle),intent(inout)::p(:) !particle type
integer,intent(in)::ntotal
double precision,intent(in)::h
double precision,intent(out)::corrvx(:),corrvy(:) !correct velocity
double precision::wij,dwdx,dwdy !kernal function
integer,save::i,j,k,ii,lx,ly,ii_neigh
double precision::rij,rij_x,rij_y
integer,save::ncx,ncy,nct
double precision::xmax,ymax!max domain size
integer,allocatable::ibox(:,:),ll(:),np(:)
integer,save::icell,kcell
integer,save::ndx(8),ndz(8) !no.of neighbour cells
double precision::s1,s2,v1,v2,v3
integer,parameter::max_part=200
integer::j1,j2,jj


xmax=10.0 !assuming square cell-rectangular cell is also possible
ymax=10.0

ncx=nint(xmax/(3.0*h))
ncy=nint(ymax/(3.0*h))
nct=ncy*ncx


allocate(ibox(nct,max_part),ll(ntotal),np(nct))



    !initialize cells
      
  ibox=0
  
   do i=1,nct
   np(i)=0
    end do 

  do i=1,ntotal
    corrvx(i)=0.0d0
    corrvy(i)=0.0d0
  end do
    
     !sorting particles into cells
     do k=1,ntotal
       icell=int(p(k)%coord%x(1)/(3.0*h))+1
       kcell=int(p(k)%coord%x(2)/(3.0*h))+1
       ii=icell+(kcell-1)*ncx
       np(ii)=np(ii)+1 
       ibox(ii,np(ii))=k      
      end do


    
    !sorting particles into cells
!$$$$$$       do k=1,ntotal      
!$$$$$$       icell=int(p(k)%coord%x(1)/(3.0*h))+1
!$$$$$$       kcell=int(p(k)%coord%x(2)/(3.0*h))+1              
!$$$$$$       ii=icell+(kcell-1)*ncx
!$$$$$$       ll(k)=ibox(ii)
!$$$$$$       ibox(ii)=k
!$$$$$$       np(ii)=np(ii)+1 !no.of particles in each cell      
!$$$$$$       end do


 !loop over all cells

       do ly=1,ncy
          do lx=1,ncx
           ii=lx+(ly-1)*ncx
       if(np(ii).gt.0)then!there are particles in the cell
          do j1=1,np(ii)
        i=ibox(ii,j1)       
    if(p(i)%id.eq.1)then!fluid particles only
      s1=0.0
      s2=0.0

     do j2=1,np(ii)           
         j=ibox(ii,j2)
         !if(i.ne.j)then

        if(i.ne.j.and.p(j)%id.eq.1)then !fluid particles only                              
       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x**2+rij_y**2
       if(rij.le.9.0*h*h)then
               
       call ComKer(rij,rij_x,rij_y,h,wij,dwdx,dwdy)
             
         v1=2.0*p(j)%mass/(p(i)%dens+p(j)%dens)
         v2=(p(j)%vel%x(1)-p(i)%vel%x(1))*wij
         v3=(p(j)%vel%x(2)-p(i)%vel%x(2))*wij
         s1=s1+v1*v2
         s2=s2+v1*v3

            end if
              end if
             
              end do

  !searcing over neighbour cells
        
          ndx=(/1,0,0,1,-1,-1,1,-1/)
          ndz=(/0,1,-1,1,1,0,-1,-1/)
              
       do k=1,8

        icell=lx+ndx(k)
        kcell=ly+ndz(k)
        !PERIOD bc

       if(icell.lt.1) cycle

       if(icell.gt.ncx) cycle
        
       if(kcell.lt.1) cycle
        
       if(kcell.gt.ncy) cycle
        
       ii_neigh=icell+(kcell-1)*ncx
            
        if(np(ii_neigh).gt.0)then
          do jj=1,np(ii_neigh)
             j=ibox(ii_neigh,jj)
            if(p(j)%id.eq.1)then !fluid particles only
             
       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x**2+rij_y**2
       if(rij.le.9.0*h*h)then
          
      call ComKer(rij,rij_x,rij_y,h,wij,dwdx,dwdy)

         v1=2.0*p(j)%mass/(p(i)%dens+p(j)%dens)
         v2=(p(j)%vel%x(1)-p(i)%vel%x(1))*wij
         v3=(p(j)%vel%x(2)-p(i)%vel%x(2))*wij
         s1=s1+v1*v2
         s2=s2+v1*v3
    
           
          end if
          
              end if
            
            end do
                end if
                     end do

                          
   corrvx(i)=p(i)%vel%x(1)+0.5*s1
   corrvy(i)=p(i)%vel%x(2)+0.5*s2
   
  end if
                   
    end do       
      end if         
         end do
             end do


















!$$$$$$        do ly=1,ncy
!$$$$$$           do lx=1,ncx
!$$$$$$            ii=lx+(ly-1)*ncx
!$$$$$$        if(ibox(ii).gt.0)then!there are particles in the cell
!$$$$$$          i=ibox(ii)
!$$$$$$                          
!$$$$$$          do while(i.gt.0)
!$$$$$$         
!$$$$$$            j=ll(i)
!$$$$$$            do while(j.gt.0) !loop in same cell
!$$$$$$            
!$$$$$$        rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
!$$$$$$        rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
!$$$$$$        rij=rij_x**2+rij_y**2
!$$$$$$ 
!$$$$$$         if(rij.le.9.0*h*h)then
!$$$$$$           if(p(j)%id.eq.1)then !fluid particles only
!$$$$$$             
!$$$$$$             call ComKer(rij,rij_x,rij_y,h,wij,dwdx,dwdy)         
!$$$$$$             v1=2.0*p(j)%mass/(p(i)%dens+p(j)%dens)
!$$$$$$             v2=(p(j)%vel%x(1)-p(i)%vel%x(1))*wij 
!$$$$$$             v3=(p(j)%vel%x(2)-p(i)%vel%x(2))*wij
!$$$$$$             corrvx(i)=corrvx(i)+0.5*v1*v2
!$$$$$$             corrvy(i)=corrvy(i)+0.5*v1*v3
!$$$$$$             corrvx(j)=corrvx(j)-0.5*v1*v2
!$$$$$$             corrvy(j)=corrvy(j)-0.5*v1*v3
!$$$$$$      
!$$$$$$              end if
!$$$$$$             end if
!$$$$$$            j=ll(j)
!$$$$$$          end do
!$$$$$$ 
!$$$$$$        !now for neighbour cells
!$$$$$$ 
!$$$$$$           ndx=(/1,0,1,-1/)
!$$$$$$           ndz=(/0,1,1, 1/)
!$$$$$$ 
!$$$$$$         do k=1,4
!$$$$$$         icell=lx+ndx(k)
!$$$$$$         kcell=ly+ndz(k)
!$$$$$$ 
!$$$$$$         !PERIOD bc      
!$$$$$$ 
!$$$$$$         if(icell.lt.1) cycle
!$$$$$$         
!$$$$$$         if(icell.gt.ncx) cycle
!$$$$$$ 
!$$$$$$         if(kcell.lt.1) cycle
!$$$$$$       
!$$$$$$         if(kcell.gt.ncy) cycle
!$$$$$$ 
!$$$$$$        ii_neigh=icell+(kcell-1)*ncx
!$$$$$$  
!$$$$$$         if(np(ii_neigh).gt.0)then
!$$$$$$            j=ibox(ii_neigh)
!$$$$$$          do while(j.gt.0)
!$$$$$$ 
!$$$$$$        rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
!$$$$$$        rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
!$$$$$$        rij=rij_x**2+rij_y**2
!$$$$$$ 
!$$$$$$         if(rij.le.9.0*h*h)then
!$$$$$$            if(p(j)%id.eq.1)then !fluid particles only
!$$$$$$           call ComKer(rij,rij_x,rij_y,h,wij,dwdx,dwdy)
!$$$$$$           
!$$$$$$             v1=2.0*p(j)%mass/(p(i)%dens+p(j)%dens)
!$$$$$$             v2=(p(j)%vel%x(1)-p(i)%vel%x(1))*wij 
!$$$$$$             v3=(p(j)%vel%x(2)-p(i)%vel%x(2))*wij
!$$$$$$             corrvx(i)=corrvx(i)+0.5*v1*v2
!$$$$$$             corrvy(i)=corrvy(i)+0.5*v1*v3
!$$$$$$             corrvx(j)=corrvx(j)-0.5*v1*v2
!$$$$$$             corrvy(j)=corrvy(j)-0.5*v1*v3
!$$$$$$              end if
!$$$$$$                 end if
!$$$$$$ 
!$$$$$$        j=ll(j)
!$$$$$$    
!$$$$$$      end do
!$$$$$$ 
!$$$$$$        end if
!$$$$$$ 
!$$$$$$          end do
!$$$$$$      
!$$$$$$   
!$$$$$$         i=ll(i)
!$$$$$$           end do 
!$$$$$$          
!$$$$$$               end if
!$$$$$$                end do
!$$$$$$                      end do
!$$$$$$ 
!$$$$$$ 
!$$$$$$   do i=1,ntotal
!$$$$$$    if(p(i)%id.eq.1)then
!$$$$$$      corrvx(i)=corrvx(i)+p(i)%vel%x(1)
!$$$$$$      corrvy(i)=corrvy(i)+p(i)%vel%x(2)
!$$$$$$      end if
!$$$$$$     end do





     end subroutine
   end module
