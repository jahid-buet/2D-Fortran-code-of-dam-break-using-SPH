module pressure

contains
subroutine Compute_press_fluid(p,ntotal,rho,coefsound,hswl)
use particles  
implicit none
integer,intent(in)::ntotal
type(particle),intent(inout)::p(:)
double precision,intent(in)::rho !reference density of water
double precision::cs0 !artificial soundspped
double precision,intent(in)::hswl !still water level
integer,parameter::ceff=7
double precision,parameter::g=9.81
double precision,intent(in)::coefsound
integer::i
double precision::B


 !initialize pressures between particles(for all particles)
 do i=1,ntotal
  p(i)%press=0.0
 end do

!compuing pressure betwn both fluid particles particles
   do i=1,ntotal
   if(p(i)%id.eq.1)then !fluid particles only  
     !compute soundspeed
     cs0=coefsound*sqrt(g*hswl)
     B=(cs0*cs0*rho/7.0)
    p(i)%press=B*((p(i)%dens/rho)**ceff-1.0d0) !dynamic pressure between fluid particles
   end if
     end do
 end subroutine
     
 

  
              
  !compute pressure and velocity of boundary particles by extrapolation
  subroutine Compute_press_and_vel_bound(p,ntotal,h,dp,pr_vel,ep_vel,length,height,kernal_type)
  use particles
  use kernal 
  implicit none
  integer,intent(in)::ntotal
  type(particle),intent(inout)::p(:)
  double precision,intent(in)::dp !particles spacing
  type(vector_2d),intent(out)::ep_vel(:)!extrapolated wall velocity
  type(vector_2d),intent(in)::pr_vel(:) !prescribed wall velocity
  integer::i,j
  integer,intent(in)::kernal_type
  double precision,intent(in)::h
  double precision,intent(in)::length,height
  double precision::gra(2)
  integer,allocatable::box(:,:),np(:)
  integer,save::ncellx(8),ncelly(8) !no.of neighbour cells in x and y direction
  double precision::xmax,ymax!max domain size
  integer::ix,iy,neigh_cell,i1,i2
  integer::cellx,celly,cell_total
  integer::j1,j2,jj,sc_index,k
  double precision::wij,dwdx,dwdy  
  double precision::sum1,sum2,sum3,wij_sum
  double precision::ex_var1,ex_var2
  double precision,allocatable::pa1(:),pa2(:),pa3(:),wij_ex(:)
  double precision::rij_x,rij_y,rij
  double precision::kf 
  integer,parameter::max_part=400


   

  


  if(kernal_type.eq.1)then
   kf=3.0d0
   else
   kf=2.0d0
  end if

   xmax=length+kf*h !maximum domain extension in x and y direction for cell-division
   ymax=height+kf*h 

  cellx=nint(xmax/(kf*h)) !no. of cell in x direction
  celly=nint(ymax/(kf*h)) !no. of cell in y direction
  cell_total=cellx*celly  !total no. of cell

  allocate(box(cell_total,max_part),np(cell_total))
  allocate(pa1(ntotal),pa2(ntotal),pa3(ntotal),wij_ex(ntotal))
 !initialize gravity vector
  gra(1)=0.0
  gra(2)=-9.81
 


!initialize kernal value
wij=0.0
dwdx=0.0
dwdy=0.0



!initialize cell parameter 
  box=0

  do i=1,cell_total
   np(i)=0
  end do
  
 
 !initialize pressure of boundary
   do i=1,ntotal
   if(p(i)%id.gt.1)then
      p(i)%press=0.0d0
      end if
    end do



 
!sorting particles into each cell
     do k=1,ntotal
      ix=int(p(k)%coord%x(1)/(kf*h))+1
      iy=int(p(k)%coord%x(2)/(kf*h))+1
      sc_index=ix+(iy-1)*cellx
      np(sc_index)=np(sc_index)+1 
      box(sc_index,np(sc_index))=k      
      end do



 !linked-cell list
       do i2=1,celly
          do i1=1,cellx
           sc_index=i1+(i2-1)*cellx
       if(np(sc_index).gt.0)then!there are particles in the cell
          do j1=1,np(sc_index)
        i=box(sc_index,j1)       
      if(p(i)%id.gt.1)then !boundary particles only
        sum1=0.0
        sum2=0.0
        sum3=0.0
        ex_var1=0.0
        ex_var2=0.0
        wij_sum=0.0                                                   !id=1 is fluid partilces
                                                                      !id >1 is boundary particles
         do j2=1,np(sc_index)                      
         j=box(sc_index,j2)
        if(p(j)%id.eq.1)then !fluid particles only                             
       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x*rij_x +rij_y*rij_y
        if(sqrt(rij).le.kf*h)then
               
       call ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)        
  

    wij_sum=wij_sum+wij 
    sum1=sum1+(p(j)%press*wij)
    sum2=sum2+(p(j)%dens*rij_x*wij)
    sum3=sum3+(p(j)%dens*rij_y*wij)
    ex_var1=ex_var1+(p(j)%vel%x(1)*wij)
    ex_var2=ex_var2+(p(j)%vel%x(2)*wij)
     
            end if
              end if
             
              end do

   !searcing over neighbour cells(all 8 cells must be searched for extrapolation of pressure values)
        
          ncellx=(/1,0,0,1,-1,-1,1,-1/)
          ncelly=(/0,1,-1,1,1,0,-1,-1/)
              
       do k=1,8

        ix=i1+ncellx(k)
        iy=i2+ncelly(k)
        !PERIOD bc

       if(ix.lt.1) cycle

       if(ix.gt.cellx) cycle
        
       if(iy.lt.1) cycle
        
       if(iy.gt.celly) cycle
        
       neigh_cell=ix+(iy-1)*cellx
            
        if(np(neigh_cell).gt.0)then
          do jj=1,np(neigh_cell)
             j=box(neigh_cell,jj)
             
           if(p(j)%id.eq.1)then !fluid particles only
           rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
           rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
           rij=rij_x**2+rij_y**2
      
       if(sqrt(rij).le.kf*h)then          
      call ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)


    wij_sum=wij_sum+wij 

    sum1=sum1+(p(j)%press*wij)
    sum2=sum2+(p(j)%dens*rij_x*wij)
    sum3=sum3+(p(j)%dens*rij_y*wij)
    ex_var1=ex_var1+(p(j)%vel%x(1)*wij)
    ex_var2=ex_var2+(p(j)%vel%x(2)*wij)
    
           
          end if 
             end if           
                end do
                   end if
                      end do


              

       
       
          if(wij_sum.gt.0.0)then
         !extrapolate boundary pressure
        p(i)%press=(sum1+(gra(1)*sum2)+(gra(2)*sum3))/wij_sum
         !extrapolate boundary velocity
         ep_vel(i)%x(1)=2.0*pr_vel(i)%x(1)-(ex_var1/wij_sum)
         ep_vel(i)%x(2)=2.0*pr_vel(i)%x(2)-(ex_var2/wij_sum)
          end if
   

  end if
                   
    end do       
      end if         
         end do
             end do




       end subroutine
  end module        