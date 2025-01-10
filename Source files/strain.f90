module CompStrainRate
use particles
use kernal
contains
subroutine Strain_Rate_Slip(kernal_type,visco_type,p,ntotal,h,dp,length,height)
implicit none
type(particle),intent(inout)::p(:)
double precision,intent(in)::h,dp
integer,intent(in)::kernal_type,visco_type
double precision::wij,dwdx,dwdy
integer,intent(in)::ntotal
double precision,intent(in)::length,height
double precision::rij_x,rij_y,rij
integer,save::ncellx(4),ncelly(4) !no.of neighbour cells in x and y direction
double precision::xmax,ymax!max domain size
double precision::kf
integer,allocatable::box(:),ll(:),np(:)
integer::ix,iy
integer::i,j,k,sc_index,i1,i2,neigh_cell
integer::cellx,celly,cell_total!no.of neighbour cells
double precision,allocatable::sxx(:),syy(:),sxy(:)
double precision,parameter::c_s=0.12
double precision,parameter::c_t=0.080d0
double precision::k_sps,v_t    

    
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

allocate(sxx(ntotal),syy(ntotal),sxy(ntotal))

allocate(box(cell_total),ll(ntotal),np(cell_total))


    !initialize cells
       
      do i=1,cell_total
        box(i)=0
        np(i)=0
       end do

  do i=1,ntotal
    ll(i)=0
   end do





     !sorting particles into cells
      do k=1,ntotal      
      ix=int(p(k)%coord%x(1)/(kf*h))+1
      iy=int(p(k)%coord%x(2)/(kf*h))+1              
      sc_index=ix+(iy-1)*cellx
      ll(k)=box(sc_index)
      box(sc_index)=k
      np(sc_index)=np(sc_index)+1 !no.of particles in each cell      
      end do

      !computation of strain-rate using cell-linked list
    !initialize  stress tensor
     do i=1,ntotal
      p(i)%txx=0.0d0
      p(i)%txy=0.0d0
      p(i)%tyy=0.0d0
     end do
   !initialize strain rate
   do i=1,ntotal    
   sxx(i)=0.0d0
   syy(i)=0.0d0
   sxy(i)=0.0d0 
   end do


        !loop over all cells

        do i2=1,celly
          do i1=1,cellx
           sc_index=ix+(iy-1)*cellx
       if(np(sc_index).gt.0)then!there are particles in the cell
         i=box(sc_index)
                         
         do while(i.gt.0)
        if(p(i)%id.eq.1)then !fluid particle only

           j=ll(i)
           do while(j.gt.0) !loop in same cell
          if(p(j)%id.eq.1)then !interaction with fluid only
       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x**2+rij_y**2

        if(rij.le.kf*kf*h*h)then
          call ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)
          
        
     if(visco_type.eq.1)then!laminar+sps 
       !compute sps-strain rate tensor
      sxx(i)=sxx(i)+(p(j)%mass/p(j)%dens)*(p(j)%vel%x(1)-p(i)%vel%x(1))*dwdx
      
      sxy(i)=sxy(i)+(0.5d0*(p(j)%mass/p(j)%dens)*(p(j)%vel%x(1)-p(i)%vel%x(1))*dwdy&
                   +0.5d0*(p(j)%mass/p(j)%dens)*(p(j)%vel%x(2)-p(i)%vel%x(2))*dwdx)
                   
      syy(i)=syy(i)+(p(j)%mass/p(j)%dens)*(p(j)%vel%x(2)-p(i)%vel%x(2))*dwdy
      

      sxx(j)=sxx(j)+(p(i)%mass/p(i)%dens)*(p(j)%vel%x(1)-p(i)%vel%x(1))*dwdx
      
      sxy(j)=sxy(j)+(0.5d0*(p(i)%mass/p(i)%dens)*(p(j)%vel%x(1)-p(i)%vel%x(1))*dwdy&
                   +0.50d0*(p(i)%mass/p(i)%dens)*(p(j)%vel%x(2)-p(i)%vel%x(2))*dwdx)
                   
      syy(j)=syy(j)+(p(i)%mass/p(i)%dens)*(p(j)%vel%x(2)-p(i)%vel%x(2))*dwdy
      
         end if
         

     if(visco_type.eq.3)then !orginal viscosity
        
    !compute strain xx for each particles
    sxx(i)=sxx(i)+2.0*((p(i)%mat(2,1)*(p(j)%vel%x(1)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(2,2)*(p(j)%vel%x(1)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(2,3)*(p(j)%vel%x(1)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))&

                       -(2./3.)*((p(i)%mat(2,1)*(p(j)%vel%x(1)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(2,2)*(p(j)%vel%x(1)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(2,3)*(p(j)%vel%x(1)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))&

                       -(2./3.)*((p(i)%mat(3,1)*(p(j)%vel%x(2)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(3,2)*(p(j)%vel%x(2)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(2)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))


   
    sxx(j)=sxx(j)+2.0*((p(j)%mat(2,1)*(p(i)%vel%x(1)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(2,2)*(p(i)%vel%x(1)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(2,3)*(p(i)%vel%x(1)*wij*rij_y*p(i)%mass/p(i)%dens)))&

                       -(2./3.)*((p(j)%mat(2,1)*(p(i)%vel%x(1)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(2,2)*(p(i)%vel%x(1)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(2,3)*(p(i)%vel%x(1)*wij*rij_y*p(i)%mass/p(i)%dens)))&

                       -(2./3.)*((p(j)%mat(3,1)*(p(i)%vel%x(2)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(3,2)*(p(i)%vel%x(2)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(3,3)*(p(i)%vel%x(2)*wij*rij_y*p(i)%mass/p(i)%dens)))


  



 !compute strain zz component

 syy(i)=syy(i)+2.0*((p(i)%mat(3,1)*(p(j)%vel%x(2)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(3,2)*(p(j)%vel%x(2)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(2)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))&

                       -(2./3.)*((p(i)%mat(2,1)*(p(j)%vel%x(1)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(2,2)*(p(j)%vel%x(1)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(2,3)*(p(j)%vel%x(1)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))&

                       -(2./3.)*((p(i)%mat(3,1)*(p(j)%vel%x(2)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(3,2)*(p(j)%vel%x(2)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(2)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))





syy(j)=syy(j)+2.0*((p(j)%mat(3,1)*(p(i)%vel%x(2)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(3,2)*(p(i)%vel%x(2)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(3,3)*(p(i)%vel%x(2)*wij*rij_y*p(i)%mass/p(i)%dens)))&

                       -(2./3.)*((p(j)%mat(2,1)*(p(i)%vel%x(1)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(2,2)*(p(i)%vel%x(1)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(2,3)*(p(i)%vel%x(1)*wij*rij_y*p(i)%mass/p(i)%dens)))&

                       -(2./3.)*((p(j)%mat(3,1)*(p(i)%vel%x(2)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(3,2)*(p(i)%vel%x(2)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(3,3)*(p(i)%vel%x(2)*wij*rij_y*p(i)%mass/p(i)%dens)))

                       


!compute strain xz component

sxy(i)=sxy(i)+(p(i)%mat(3,1)*(p(j)%vel%x(1)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(3,2)*(p(j)%vel%x(1)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(1)*wij*(-rij_y)*p(j)%mass/p(j)%dens))&

                       +(p(i)%mat(2,1)*(p(j)%vel%x(2)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(2,2)*(p(j)%vel%x(2)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(2,3)*(p(j)%vel%x(2)*wij*(-rij_y)*p(j)%mass/p(j)%dens))



sxy(j)=sxy(j)+(p(j)%mat(3,1)*(p(i)%vel%x(1)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(3,2)*(p(i)%vel%x(1)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(1)*wij*rij_y*p(j)%mass/p(j)%dens))&

                       +(p(j)%mat(2,1)*(p(i)%vel%x(2)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(2,2)*(p(i)%vel%x(2)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(2,3)*(p(i)%vel%x(2)*wij*rij_y*p(i)%mass/p(i)%dens))


           end if
               


            end if
             end if
           j=ll(j)
         end do

       !now for neighbour cells

          ncellx=(/1,0,1,-1/)
          ncelly=(/0,1,1, 1/)

        do k=1,4
        ix=i1+ncellx(k)
        iy=i2+ncelly(k)

        !PERIOD bc

       

        if(ix.lt.1) cycle
        
        if(ix.gt.cellx) cycle

        if(iy.lt.1) cycle
      
        if(iy.gt.celly) cycle

       neigh_cell=ix+(iy-1)*cellx
 
        if(np(neigh_cell).gt.0)then
           j=box(neigh_cell)
         do while(j.gt.0)
          if(p(j)%id.eq.1)then !intraction with fluid only
       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x**2+rij_y**2

        if(rij.le.kf*kf*h*h)then
          call ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)
         
           
       if(visco_type.eq.1)then!laminar+sps 
      !compute sps-strain rate tensor
      sxx(i)=sxx(i)+(p(j)%mass/p(j)%dens)*(p(j)%vel%x(1)-p(i)%vel%x(1))*dwdx
      
      sxy(i)=sxy(i)+(0.5d0*(p(j)%mass/p(j)%dens)*(p(j)%vel%x(1)-p(i)%vel%x(1))*dwdy&
                   +0.5d0*(p(j)%mass/p(j)%dens)*(p(j)%vel%x(2)-p(i)%vel%x(2))*dwdx)
                   
      syy(i)=syy(i)+(p(j)%mass/p(j)%dens)*(p(j)%vel%x(2)-p(i)%vel%x(2))*dwdy
      

      sxx(j)=sxx(j)+(p(i)%mass/p(i)%dens)*(p(j)%vel%x(1)-p(i)%vel%x(1))*dwdx
      
      sxy(j)=sxy(j)+(0.5d0*(p(i)%mass/p(i)%dens)*(p(j)%vel%x(1)-p(i)%vel%x(1))*dwdy&
                   +0.50d0*(p(i)%mass/p(i)%dens)*(p(j)%vel%x(2)-p(i)%vel%x(2))*dwdx)
                   
      syy(j)=syy(j)+(p(i)%mass/p(i)%dens)*(p(j)%vel%x(2)-p(i)%vel%x(2))*dwdy
       end if



           if(visco_type.eq.3)then !orginal viscosity
           
     !compute strain xx for each particles
    sxx(i)=sxx(i)+2.0*((p(i)%mat(2,1)*(p(j)%vel%x(1)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(2,2)*(p(j)%vel%x(1)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(2,3)*(p(j)%vel%x(1)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))&

                       -(2./3.)*((p(i)%mat(2,1)*(p(j)%vel%x(1)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(2,2)*(p(j)%vel%x(1)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(2,3)*(p(j)%vel%x(1)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))&

                       -(2./3.)*((p(i)%mat(3,1)*(p(j)%vel%x(2)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(3,2)*(p(j)%vel%x(2)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(2)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))


   
    sxx(j)=sxx(j)+2.0*((p(j)%mat(2,1)*(p(i)%vel%x(1)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(2,2)*(p(i)%vel%x(1)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(2,3)*(p(i)%vel%x(1)*wij*rij_y*p(i)%mass/p(i)%dens)))&

                       -(2./3.)*((p(j)%mat(2,1)*(p(i)%vel%x(1)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(2,2)*(p(i)%vel%x(1)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(2,3)*(p(i)%vel%x(1)*wij*rij_y*p(i)%mass/p(i)%dens)))&

                       -(2./3.)*((p(j)%mat(3,1)*(p(i)%vel%x(2)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(3,2)*(p(i)%vel%x(2)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(3,3)*(p(i)%vel%x(2)*wij*rij_y*p(i)%mass/p(i)%dens)))


  



 !compute strain zz component

 syy(i)=syy(i)+2.0*((p(i)%mat(3,1)*(p(j)%vel%x(2)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(3,2)*(p(j)%vel%x(2)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(2)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))&

                       -(2./3.)*((p(i)%mat(2,1)*(p(j)%vel%x(1)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(2,2)*(p(j)%vel%x(1)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(2,3)*(p(j)%vel%x(1)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))&

                       -(2./3.)*((p(i)%mat(3,1)*(p(j)%vel%x(2)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(3,2)*(p(j)%vel%x(2)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(2)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))





syy(j)=syy(j)+2.0*((p(j)%mat(3,1)*(p(i)%vel%x(2)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(3,2)*(p(i)%vel%x(2)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(3,3)*(p(i)%vel%x(2)*wij*rij_y*p(i)%mass/p(i)%dens)))&

                       -(2./3.)*((p(j)%mat(2,1)*(p(i)%vel%x(1)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(2,2)*(p(i)%vel%x(1)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(2,3)*(p(i)%vel%x(1)*wij*rij_y*p(i)%mass/p(i)%dens)))&

                       -(2./3.)*((p(j)%mat(3,1)*(p(i)%vel%x(2)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(3,2)*(p(i)%vel%x(2)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(3,3)*(p(i)%vel%x(2)*wij*rij_y*p(i)%mass/p(i)%dens)))

                       


!compute strain xz component

sxy(i)=sxy(i)+(p(i)%mat(3,1)*(p(j)%vel%x(1)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(3,2)*(p(j)%vel%x(1)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(1)*wij*(-rij_y)*p(j)%mass/p(j)%dens))&

                       +(p(i)%mat(2,1)*(p(j)%vel%x(2)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(2,2)*(p(j)%vel%x(2)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(2,3)*(p(j)%vel%x(2)*wij*(-rij_y)*p(j)%mass/p(j)%dens))



sxy(j)=sxy(j)+(p(j)%mat(3,1)*(p(i)%vel%x(1)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(3,2)*(p(i)%vel%x(1)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(1)*wij*rij_y*p(j)%mass/p(j)%dens))&

                       +(p(j)%mat(2,1)*(p(i)%vel%x(2)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(2,2)*(p(i)%vel%x(2)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(2,3)*(p(i)%vel%x(2)*wij*rij_y*p(i)%mass/p(i)%dens))



            end if


     
          end if
            end if

       j=ll(j)
   
     end do

       end if

         end do
     
           end if

         i=ll(i)
          end do 
         
              end if
               end do
                     end do

do i=1,ntotal
  if(p(i)%id.eq.1.and.visco_type.eq.1)then !fluid particles only and laminar+sps viscosity
    v_t=(c_s*dp)**2*sqrt(2.0*sxx(i)*sxx(i)+2.0*syy(i)*syy(i)+4.0*sxy(i)*sxy(i))
    k_sps=(v_t/(c_t*dp))**2 
 !compute stress tensor xx for each particles
p(i)%txx=2.0*v_t*sxx(i)*p(i)%dens-(2.0/3.0)*k_sps*p(i)%dens

!compute stress tensor xy component
p(i)%txy=2.0*v_t*sxy(i)*p(i)%dens

!compute stress tensor yy component
p(i)%tyy=2.0*v_t*syy(i)*p(i)%dens-(2.0/3.0)*k_sps*p(i)%dens

 end if
 
if(p(i)%id.eq.1.and.visco_type.eq.3)then
  p(i)%hxx=sxx(i)
  p(i)%hxy=sxy(i)
  p(i)%hyy=syy(i)
  end if
   
        
      
 end do

    end subroutine



subroutine Strain_Rate_NoSlip(kernal_type,visco_type,p,ntotal,h,dp,ep_vel,length,height)
implicit none
type(particle),intent(inout)::p(:)
double precision,intent(in)::h,dp
integer,intent(in)::kernal_type,visco_type
double precision::wij,dwdx,dwdy
integer,intent(in)::ntotal
double precision,intent(in)::length,height
type(vector_2d),intent(in)::ep_vel(:)
double precision::rij_x,rij_y,rij
integer,save::ncellx(4),ncelly(4) !no.of neighbour cells in x and y direction
double precision::xmax,ymax!max domain size
double precision::kf
integer,allocatable::box(:),ll(:),np(:)
integer::ix,iy
integer::i,j,k,sc_index,i1,i2,neigh_cell
integer::cellx,celly,cell_total!no.of neighbour cells
double precision,allocatable::sxx(:),syy(:),sxy(:)
double precision,parameter::c_s=0.12
double precision,parameter::c_t=0.080d0
double precision::k_sps,v_t    





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

allocate(sxx(ntotal),syy(ntotal),sxy(ntotal))

allocate(box(cell_total),ll(ntotal),np(cell_total))


    !initialize cells
       
      do i=1,cell_total
        box(i)=0
        np(i)=0
       end do

  do i=1,ntotal
    ll(i)=0
   end do


     !sorting particles into cells
      do k=1,ntotal      
      ix=int(p(k)%coord%x(1)/(kf*h))+1
      iy=int(p(k)%coord%x(2)/(kf*h))+1              
      sc_index=ix+(iy-1)*cellx
      ll(k)=box(sc_index)
      box(sc_index)=k
      np(sc_index)=np(sc_index)+1 !no.of particles in each cell      
      end do

      !computation of strain-rate using cell-linked list

  !initialize stress-tensor
  do i=1,ntotal
    p(i)%txx=0.0d0
    p(i)%txy=0.0d0
    p(i)%tyy=0.0d0
  end do
 !initialize strain rate tensor
   do i=1,ntotal    
   sxx(i)=0.0d0
   syy(i)=0.0d0
   sxy(i)=0.0d0 
   end do

 !initialize wall velocity to extrapolated velocity
  do i=1,ntotal
    if(p(i)%id.gt.1)then
      p(i)%vel%x(1)=ep_vel(i)%x(1)
      p(i)%vel%x(2)=ep_vel(i)%x(2)
      end if
      end do
          
        !loop over all cells

         do i2=1,celly
          do i1=1,cellx
           sc_index=i1+(i2-1)*cellx
       if(box(sc_index).gt.0)then!there are particles in the cell
         i=box(sc_index)
                         
         do while(i.gt.0)
        if(p(i)%id.eq.1)then !fluid particle only

           j=ll(i)
           do while(j.gt.0) !loop in same cell
          
       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x**2+rij_y**2

        if(rij.le.kf*kf*h*h)then
          call ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)
          
        
        if(visco_type.eq.1)then
       !compute sps-strain rate tensor
      sxx(i)=sxx(i)+(p(j)%mass/p(j)%dens)*(p(j)%vel%x(1)-p(i)%vel%x(1))*dwdx
      
      sxy(i)=sxy(i)+(0.5d0*(p(j)%mass/p(j)%dens)*(p(j)%vel%x(1)-p(i)%vel%x(1))*dwdy&
                   +0.5d0*(p(j)%mass/p(j)%dens)*(p(j)%vel%x(2)-p(i)%vel%x(2))*dwdx)
                   
      syy(i)=syy(i)+(p(j)%mass/p(j)%dens)*(p(j)%vel%x(2)-p(i)%vel%x(2))*dwdy
      

      sxx(j)=sxx(j)+(p(i)%mass/p(i)%dens)*(p(j)%vel%x(1)-p(i)%vel%x(1))*dwdx
      
      sxy(j)=sxy(j)+(0.5d0*(p(i)%mass/p(i)%dens)*(p(j)%vel%x(1)-p(i)%vel%x(1))*dwdy&
                   +0.50d0*(p(i)%mass/p(i)%dens)*(p(j)%vel%x(2)-p(i)%vel%x(2))*dwdx)
                   
      syy(j)=syy(j)+(p(i)%mass/p(i)%dens)*(p(j)%vel%x(2)-p(i)%vel%x(2))*dwdy
       end if
       

        if(visco_type.eq.3)then !orginal viscosity
            !compute strain xx for each particles
    sxx(i)=sxx(i)+2.0*((p(i)%mat(2,1)*(p(j)%vel%x(1)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(2,2)*(p(j)%vel%x(1)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(2,3)*(p(j)%vel%x(1)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))&

                       -(2./3.)*((p(i)%mat(2,1)*(p(j)%vel%x(1)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(2,2)*(p(j)%vel%x(1)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(2,3)*(p(j)%vel%x(1)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))&

                       -(2./3.)*((p(i)%mat(3,1)*(p(j)%vel%x(2)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(3,2)*(p(j)%vel%x(2)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(2)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))


   
    sxx(j)=sxx(j)+2.0*((p(j)%mat(2,1)*(p(i)%vel%x(1)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(2,2)*(p(i)%vel%x(1)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(2,3)*(p(i)%vel%x(1)*wij*rij_y*p(i)%mass/p(i)%dens)))&

                       -(2./3.)*((p(j)%mat(2,1)*(p(i)%vel%x(1)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(2,2)*(p(i)%vel%x(1)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(2,3)*(p(i)%vel%x(1)*wij*rij_y*p(i)%mass/p(i)%dens)))&

                       -(2./3.)*((p(j)%mat(3,1)*(p(i)%vel%x(2)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(3,2)*(p(i)%vel%x(2)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(3,3)*(p(i)%vel%x(2)*wij*rij_y*p(i)%mass/p(i)%dens)))


  



 !compute strain zz component

 syy(i)=syy(i)+2.0*((p(i)%mat(3,1)*(p(j)%vel%x(2)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(3,2)*(p(j)%vel%x(2)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(2)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))&

                       -(2./3.)*((p(i)%mat(2,1)*(p(j)%vel%x(1)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(2,2)*(p(j)%vel%x(1)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(2,3)*(p(j)%vel%x(1)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))&

                       -(2./3.)*((p(i)%mat(3,1)*(p(j)%vel%x(2)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(3,2)*(p(j)%vel%x(2)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(2)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))





syy(j)=syy(j)+2.0*((p(j)%mat(3,1)*(p(i)%vel%x(2)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(3,2)*(p(i)%vel%x(2)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(3,3)*(p(i)%vel%x(2)*wij*rij_y*p(i)%mass/p(i)%dens)))&

                       -(2./3.)*((p(j)%mat(2,1)*(p(i)%vel%x(1)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(2,2)*(p(i)%vel%x(1)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(2,3)*(p(i)%vel%x(1)*wij*rij_y*p(i)%mass/p(i)%dens)))&

                       -(2./3.)*((p(j)%mat(3,1)*(p(i)%vel%x(2)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(3,2)*(p(i)%vel%x(2)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(3,3)*(p(i)%vel%x(2)*wij*rij_y*p(i)%mass/p(i)%dens)))

                       


!compute strain xz component

sxy(i)=sxy(i)+(p(i)%mat(3,1)*(p(j)%vel%x(1)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(3,2)*(p(j)%vel%x(1)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(1)*wij*(-rij_y)*p(j)%mass/p(j)%dens))&

                       +(p(i)%mat(2,1)*(p(j)%vel%x(2)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(2,2)*(p(j)%vel%x(2)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(2,3)*(p(j)%vel%x(2)*wij*(-rij_y)*p(j)%mass/p(j)%dens))



sxy(j)=sxy(j)+(p(j)%mat(3,1)*(p(i)%vel%x(1)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(3,2)*(p(i)%vel%x(1)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(1)*wij*rij_y*p(j)%mass/p(j)%dens))&

                       +(p(j)%mat(2,1)*(p(i)%vel%x(2)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(2,2)*(p(i)%vel%x(2)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(2,3)*(p(i)%vel%x(2)*wij*rij_y*p(i)%mass/p(i)%dens))
           



          end if




 


            end if
             
           j=ll(j)
         end do

       !now for neighbour cells

         ncellx=(/1,0,1,-1/)
         ncelly=(/0,1,1, 1/)

        do k=1,4
        ix=i1+ncellx(k)
        iy=i2+ncelly(k)
        !PERIOD bc

       

        if(ix.lt.1) cycle
        
        if(ix.gt.cellx) cycle

        if(iy.lt.1) cycle
      
        if(iy.gt.celly) cycle

       neigh_cell=ix+(iy-1)*cellx
 
        if(np( neigh_cell).gt.0)then
           j=box( neigh_cell)
         do while(j.gt.0)
          
       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x**2+rij_y**2

        if(rij.le.kf*kf*h*h)then
          call ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)
         
           
       if(visco_type.eq.1)then!laminar+sps
      !compute sps-strain rate tensor
      sxx(i)=sxx(i)+(p(j)%mass/p(j)%dens)*(p(j)%vel%x(1)-p(i)%vel%x(1))*dwdx
      
      sxy(i)=sxy(i)+(0.5d0*(p(j)%mass/p(j)%dens)*(p(j)%vel%x(1)-p(i)%vel%x(1))*dwdy&
                   +0.5d0*(p(j)%mass/p(j)%dens)*(p(j)%vel%x(2)-p(i)%vel%x(2))*dwdx)
                   
      syy(i)=syy(i)+(p(j)%mass/p(j)%dens)*(p(j)%vel%x(2)-p(i)%vel%x(2))*dwdy
      

      sxx(j)=sxx(j)+(p(i)%mass/p(i)%dens)*(p(j)%vel%x(1)-p(i)%vel%x(1))*dwdx
      
      sxy(j)=sxy(j)+(0.5d0*(p(i)%mass/p(i)%dens)*(p(j)%vel%x(1)-p(i)%vel%x(1))*dwdy&
                   +0.50d0*(p(i)%mass/p(i)%dens)*(p(j)%vel%x(2)-p(i)%vel%x(2))*dwdx)
                   
      syy(j)=syy(j)+(p(i)%mass/p(i)%dens)*(p(j)%vel%x(2)-p(i)%vel%x(2))*dwdy
           end if

     if(visco_type.eq.3)then!orginal viscosity
       !compute strain xx for each particles
    sxx(i)=sxx(i)+2.0*((p(i)%mat(2,1)*(p(j)%vel%x(1)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(2,2)*(p(j)%vel%x(1)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(2,3)*(p(j)%vel%x(1)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))&

                       -(2./3.)*((p(i)%mat(2,1)*(p(j)%vel%x(1)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(2,2)*(p(j)%vel%x(1)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(2,3)*(p(j)%vel%x(1)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))&

                       -(2./3.)*((p(i)%mat(3,1)*(p(j)%vel%x(2)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(3,2)*(p(j)%vel%x(2)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(2)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))


   
    sxx(j)=sxx(j)+2.0*((p(j)%mat(2,1)*(p(i)%vel%x(1)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(2,2)*(p(i)%vel%x(1)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(2,3)*(p(i)%vel%x(1)*wij*rij_y*p(i)%mass/p(i)%dens)))&

                       -(2./3.)*((p(j)%mat(2,1)*(p(i)%vel%x(1)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(2,2)*(p(i)%vel%x(1)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(2,3)*(p(i)%vel%x(1)*wij*rij_y*p(i)%mass/p(i)%dens)))&

                       -(2./3.)*((p(j)%mat(3,1)*(p(i)%vel%x(2)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(3,2)*(p(i)%vel%x(2)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(3,3)*(p(i)%vel%x(2)*wij*rij_y*p(i)%mass/p(i)%dens)))


  



 !compute strain zz component

 syy(i)=syy(i)+2.0*((p(i)%mat(3,1)*(p(j)%vel%x(2)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(3,2)*(p(j)%vel%x(2)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(2)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))&

                       -(2./3.)*((p(i)%mat(2,1)*(p(j)%vel%x(1)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(2,2)*(p(j)%vel%x(1)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(2,3)*(p(j)%vel%x(1)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))&

                       -(2./3.)*((p(i)%mat(3,1)*(p(j)%vel%x(2)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(3,2)*(p(j)%vel%x(2)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(2)*wij*(-rij_y)*p(j)%mass/p(j)%dens)))





syy(j)=syy(j)+2.0*((p(j)%mat(3,1)*(p(i)%vel%x(2)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(3,2)*(p(i)%vel%x(2)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(3,3)*(p(i)%vel%x(2)*wij*rij_y*p(i)%mass/p(i)%dens)))&

                       -(2./3.)*((p(j)%mat(2,1)*(p(i)%vel%x(1)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(2,2)*(p(i)%vel%x(1)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(2,3)*(p(i)%vel%x(1)*wij*rij_y*p(i)%mass/p(i)%dens)))&

                       -(2./3.)*((p(j)%mat(3,1)*(p(i)%vel%x(2)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(3,2)*(p(i)%vel%x(2)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(3,3)*(p(i)%vel%x(2)*wij*rij_y*p(i)%mass/p(i)%dens)))

                       


!compute strain xz component

sxy(i)=sxy(i)+(p(i)%mat(3,1)*(p(j)%vel%x(1)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(3,2)*(p(j)%vel%x(1)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(1)*wij*(-rij_y)*p(j)%mass/p(j)%dens))&

                       +(p(i)%mat(2,1)*(p(j)%vel%x(2)*wij*p(j)%mass/p(j)%dens))+(p(i)%mat(2,2)*(p(j)%vel%x(2)*wij*(-rij_x)*p(j)%mass/p(j)%dens))&
                       + (p(i)%mat(2,3)*(p(j)%vel%x(2)*wij*(-rij_y)*p(j)%mass/p(j)%dens))



sxy(j)=sxy(j)+(p(j)%mat(3,1)*(p(i)%vel%x(1)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(3,2)*(p(i)%vel%x(1)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(i)%mat(3,3)*(p(j)%vel%x(1)*wij*rij_y*p(j)%mass/p(j)%dens))&

                       +(p(j)%mat(2,1)*(p(i)%vel%x(2)*wij*p(i)%mass/p(i)%dens))+(p(j)%mat(2,2)*(p(i)%vel%x(2)*wij*rij_x*p(i)%mass/p(i)%dens))&
                       + (p(j)%mat(2,3)*(p(i)%vel%x(2)*wij*rij_y*p(i)%mass/p(i)%dens))
       

           end if




     
          end if
            

       j=ll(j)
   
     end do

       end if

         end do
     
           end if

         i=ll(i)
          end do 
         
              end if
               end do
                     end do

do i=1,ntotal
  if(p(i)%id.eq.1.and.visco_type.eq.1)then !fluid particles only
    v_t=(c_s*dp)**2*sqrt(2.0*sxx(i)*sxx(i)+2.0*syy(i)*syy(i)+4.0*sxy(i)*sxy(i))
    k_sps=(v_t/(c_t*dp))**2 
 !compute stress tensor xx for each particles
p(i)%txx=2.0*v_t*sxx(i)*p(i)%dens-(2.0/3.0)*k_sps*p(i)%dens

!compute stress tensor xy component
p(i)%txy=2.0*v_t*sxy(i)*p(i)%dens

!compute stress tensor yy component
p(i)%tyy=2.0*v_t*syy(i)*p(i)%dens-(2.0/3.0)*k_sps*p(i)%dens

 end if


  if(p(i)%id.eq.1.and.visco_type.eq.3)then
  p(i)%hxx=sxx(i)
  p(i)%hxy=sxy(i)
  p(i)%hyy=syy(i)
  end if

    
 end do

 end subroutine


   
      end module
 
                       
                       

                     
   
    
