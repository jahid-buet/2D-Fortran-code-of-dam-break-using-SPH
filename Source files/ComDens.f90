module density
contains
subroutine Compute_density_fluid(p,ntotal,h,dp,pr_vel,length,height,kernal_type,coeff,hswl,drodt,var_max) !COMPUTING DENSITY BY SPH METHOD
use particles
use kernal
implicit none
type(particle),intent(inout)::p(:) !particle type
integer,intent(in)::ntotal
double precision,intent(in)::h
double precision,intent(in)::length,height
double precision,intent(in)::dp
double precision,intent(out)::drodt(:)
type(vector_2d),intent(in)::pr_vel(:)!prescribed velocity of wall/boundary particles
double precision,intent(out)::var_max(:)!variables for computing timestep
integer,intent(in)::kernal_type

!local variables
double precision::wij,dwdx,dwdy !kernal function
integer::i,j,k,sc_index,i1,i2,neigh_cell
double precision::rij,rij_x,rij_y
integer::cellx,celly,cell_total
double precision::xmax,ymax!max domain size
integer,allocatable::box(:),ll(:),np(:)
integer::ix,iy !cell index in x and y(or z) direction
integer,save::ncellx(4),ncelly(4) !no.of neighbour cells in x and y direction
double precision::var1,var2,var3,var4,var5,var6
double precision,allocatable::sum1(:),sum2(:),sum3(:),sum4(:),sum5(:)
integer::j1,j2,jj
double precision,parameter::delta=0.2
double precision,intent(in)::coeff,hswl
double precision::cs
double precision,parameter::g=9.81
double precision::kf
double precision,allocatable::t1(:)!variables for computing timestep
  


if(kernal_type.eq.1)then
  kf=3.0d0
else
 kf=2.0d0
end if

xmax=length+kf*h !assuming square cell-rectangular cell is also possible
ymax=height+kf*h





cellx=nint(xmax/(kf*h)) !no. of cell in x direction
celly=nint(ymax/(kf*h)) !no. of cell in y direction
cell_total=cellx*celly  !total no. of cell




cs=coeff*sqrt(g*hswl)

allocate(box(cell_total),ll(ntotal),np(cell_total))
allocate(sum1(ntotal),sum2(ntotal),sum3(ntotal),sum4(ntotal),sum5(ntotal))
allocate(t1(ntotal))
!initialize kernal value
wij=0.0
dwdx=0.0
dwdy=0.0

    
   box=0
   
    do i=1,cell_total
     np(i)=0
    end do

   do i=1,ntotal
    ll(i)=0
   end do  
    
!initialize rate of density change
   do i=1,ntotal 
   drodt(i)=0.0d0 !for sph sum
   sum1(i)=0.0
   sum2(i)=0.0
   sum3(i)=0.0d0
   sum4(i)=0.0d0
   sum5(i)=0.0d0
   var_max(i)=0.0d0
   t1(i)=0.0d0
   end do
!initialize wall/boundary particles velocity to prescribed wall velocity
do i=1,ntotal
if(p(i)%id.gt.1)then
 p(i)%vel%x(1)=pr_vel(i)%x(1)
 p(i)%vel%x(2)=pr_vel(i)%x(2)
end if
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
  





 !computation of density using cell-linked list

  !call matrix(p,ntotal,h) 

       do i2=1,celly
          do i1=1,cellx
           sc_index=i1+(i2-1)*cellx
       if(np(sc_index).gt.0)then!there are particles in the cell
         i=box(sc_index)
                         
         do while(i.gt.0)
           
           j=ll(i)
           do while(j.gt.0) !loop in same cell
            !if(p(j)%id*p(i)%id.lt.4)then
       
       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x**2+rij_y**2

         if(sqrt(rij).lt.kf*h)then
         call ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)
   
      
     var1=(p(i)%vel%x(1)-p(j)%vel%x(1))*dwdx+(p(i)%vel%x(2)-p(j)%vel%x(2))*dwdy
     var2=rij_x*dwdx+rij_y*dwdy
     var3=p(i)%dens-p(j)%dens
     var4=rij+(0.01*h)**2
     var5=(p(i)%vel%x(1)-p(j)%vel%x(1))*rij_x+(p(i)%vel%x(2)-p(j)%vel%x(2))*rij_y
     sum1(i)=sum1(i)+((p(j)%mass/p(j)%dens)*var1)
     sum1(j)=sum1(j)+((p(i)%mass/p(i)%dens)*var1)
     sum2(i)=sum2(i)+(((var2*var3)/var4)*p(j)%mass/p(j)%dens)
     sum2(j)=sum2(j)-(((var2*var3)/var4)*p(i)%mass/p(i)%dens)
     
     !t1(j)=abs((h*var5)/(rij+0.01*h*h))

    !cspm correction to density
     !sum1(i)=sum1(i)+((p(j)%mass/p(j)%dens)*dwdx*(p(j)%vel%x(1)-p(i)%vel%x(1)))
     !sum1(j)=sum1(j)+((p(i)%mass/p(i)%dens)*dwdx*(p(j)%vel%x(1)-p(i)%vel%x(1)))
     !sum3(i)=sum3(i)+((p(j)%mass/p(j)%dens)*dwdx*(-rij_x))
     !sum3(j)=sum3(j)+((p(i)%mass/p(i)%dens)*dwdx*(-rij_x))

     !sum4(i)=sum4(i)+((p(j)%mass/p(j)%dens)*dwdy*(p(j)%vel%x(2)-p(i)%vel%x(2)))
     !sum4(j)=sum4(j)+((p(i)%mass/p(i)%dens)*dwdy*(p(j)%vel%x(2)-p(i)%vel%x(2)))
     !sum5(i)=sum5(i)+((p(j)%mass/p(j)%dens)*dwdy*(-rij_y))
     !sum5(j)=sum5(j)+((p(i)%mass/p(i)%dens)*dwdy*(-rij_y))



            end if

              !end if
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
      !if(p(j)%id*p(i)%id.lt.4)then
       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x**2+rij_y**2

         if(sqrt(rij).lt.kf*h)then
          call ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)
       
      
       
     var1=(p(i)%vel%x(1)-p(j)%vel%x(1))*dwdx+(p(i)%vel%x(2)-p(j)%vel%x(2))*dwdy
     var2=rij_x*dwdx+rij_y*dwdy
     var3=p(i)%dens-p(j)%dens
     var4=rij+(0.01*h)**2
     var5=(p(i)%vel%x(1)-p(j)%vel%x(1))*rij_x+(p(i)%vel%x(2)-p(j)%vel%x(2))*rij_y
     sum1(i)=sum1(i)+((p(j)%mass/p(j)%dens)*var1)
     sum1(j)=sum1(j)+((p(i)%mass/p(i)%dens)*var1)
     sum2(i)=sum2(i)+(((var2*var3)/var4)*p(j)%mass/p(j)%dens)
     sum2(j)=sum2(j)-(((var2*var3)/var4)*p(i)%mass/p(i)%dens)

      !t1(j)=abs((h*var5)/(rij+0.01*h*h)) 


      !cspm correction to density(1st order accuracy)
     !sum1(i)=sum1(i)+((p(j)%mass/p(j)%dens)*dwdx*(p(j)%vel%x(1)-p(i)%vel%x(1)))
     !sum1(j)=sum1(j)+((p(i)%mass/p(i)%dens)*dwdx*(p(j)%vel%x(1)-p(i)%vel%x(1)))
     !sum3(i)=sum3(i)+((p(j)%mass/p(j)%dens)*dwdx*(-rij_x))
     !sum3(j)=sum3(j)+((p(i)%mass/p(i)%dens)*dwdx*(-rij_x))

     !sum4(i)=sum4(i)+((p(j)%mass/p(j)%dens)*dwdy*(p(j)%vel%x(2)-p(i)%vel%x(2)))
     !sum4(j)=sum4(j)+((p(i)%mass/p(i)%dens)*dwdy*(p(j)%vel%x(2)-p(i)%vel%x(2)))
     !sum5(i)=sum5(i)+((p(j)%mass/p(j)%dens)*dwdy*(-rij_y))
     !sum5(j)=sum5(j)+((p(i)%mass/p(i)%dens)*dwdy*(-rij_y))
     
    
          end if

            !end if

       j=ll(j)
   
     end do

       end if

         end do
        !var_max(i)=maxval(t1)  !very time consuming          
        i=ll(i)
          end do 
         
              end if
               end do
                     end do


 do i=1,ntotal
if(p(i)%id.gt.1)then
  drodt(i)=0.0
elseif(p(i)%id.eq.1)then!fluid particles only
 drodt(i)=p(i)%dens*sum1(i)  +delta*h*cs*sum2(i) ! density diffusion term(molteni& Colagrossi 2009)
 end if
end do


end subroutine




   !compute density of boundary particles
    subroutine Compute_density_bound(p,ntotal,rho0,coefsound,hswl)
    use particles
    implicit none
    type(particle),intent(inout)::p(:)
    double precision,intent(in)::rho0!density of reference fluid particles
    double precision::cs
    double precision,parameter::g=9.81
    double precision,intent(in)::hswl
    double precision,intent(in)::coefsound
    double precision::ref_press
    integer::i
    double precision::var1,pow
    integer,intent(in)::ntotal
   
 !initialize density of boundary partucles
 do i=1,ntotal
 if(p(i)%id.gt.1)then
  p(i)%dens=rho0
  end if
  end do

 

!compute density of boundary particles
          cs=coefsound*sqrt(hswl*g) 
          ref_press=(rho0*cs*cs)/7.0
          pow=1.0/7.0

        do i=1,ntotal
         if(p(i)%id.gt.1)then !boundary particles only
                 
        var1=(p(i)%press+ref_press)
                 
        p(i)%dens=rho0*((var1/ref_press)**pow)
        if(p(i)%dens.lt.rho0)then
          p(i)%dens=rho0!this prevents sticking of fluid particles to boundary wall
          end if
                  
         end if
  
         end do

  end subroutine


    end module

 

      

     













