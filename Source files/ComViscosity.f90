module CompStress
use particles
use kernal
contains
subroutine viscosity_Slip(visco_type,kernal_type,p,ntotal,nu,alpha,coeff,h,dp,rho,length,height,hswl)
implicit none
type(particle),intent(inout)::p(:)
double precision,intent(in)::h
integer,intent(in)::ntotal
double precision,intent(in)::length,height
double precision,intent(in)::dp
!type(vector_2d),intent(in)::ep_vel(:) !extrapolated wall velocity
double precision,intent(in)::nu,hswl !kinematic viscosity
double precision,intent(in)::coeff
integer,intent(in)::visco_type !viscosity type laminar , artificial or orginal
integer,intent(in)::kernal_type

!local variables
integer::i,j,k,sc_index,i1,i2,neigh_cell
double precision::rij_x,rij_y,rij,wij,dwdx,dwdy
double precision,allocatable::dwdx_corr(:),dwdy_corr(:)
integer::cellx,celly,cell_total
double precision::xmax,ymax!max domain size
double precision::mu=1d-3 !dynamic viscosity of fluid
double precision,intent(in)::rho !reference fluid density
integer,allocatable::box(:),ll(:),np(:)
integer::ix,iy
integer,save::ncellx(4),ncelly(4) !no.of neighbour cells in x and y direction
double precision::var1,var2,var3,var4,var5,var6
double precision::tij
double precision,allocatable::sum1(:),sum2(:)
integer::j1,j2,jj
double precision,allocatable::f_viscx(:),f_viscy(:)
double precision,parameter::g=9.81
double precision,intent(in)::alpha
double precision,allocatable::c(:)
double precision::phi
double precision::kf
double precision,allocatable::t_forcex(:),t_forcey(:)
double precision,allocatable::fxx(:),fyy(:),fxy(:)

  
   
 

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



allocate(sum1(ntotal),sum2(ntotal))
allocate(box(cell_total),ll(ntotal),np(cell_total),f_viscx(ntotal),f_viscy(ntotal))
allocate(c(ntotal))
allocate(dwdx_corr(ntotal),dwdy_corr(ntotal))
allocate(t_forcex(ntotal),t_forcey(ntotal))
allocate(fxx(ntotal),fyy(ntotal),fxy(ntotal))

do i=1,ntotal
if(p(i)%id.eq.1)then!only for fluid particles
 c(i)=coeff*sqrt(g*hswl) 
else
c(i)=0.0
end if
end do


    !initialize cells
       
   box=0  
        
   do i=1,cell_total
   np(i)=0
    end do
    
  do i=1,ntotal
    ll(i)=0
   end do


!initialize laminar-viscosity forces
do i=1,ntotal
 f_viscx(i)=0.0d0
 f_viscy(i)=0.0d0
sum1(i)=0.0
sum2(i)=0.0
dwdx_corr(i)=0.0d0
dwdy_corr(i)=0.0d0
end do


!initialize turbulence-viscosity forces
do i=1,ntotal
  t_forcex(i)=0.0d0
  t_forcey(i)=0.0d0
  end do
  
    DO i=1,ntotal
    fxx(i)=0.0d0
    fxy(i)=0.0d0
    fyy(i)=0.0d0
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
 
     !computation of acceleration due to viscous interaction using cell-linked list

  !call matrix(p,ntotal,h) 

       do i2=1,celly
          do i1=1,cellx
           sc_index=ix+(iy-1)*cellx
       if(np(sc_index).gt.0)then!there are particles in the cell
         i=box(sc_index)
                         
         do while(i.gt.0)

          if(p(i)%id.eq.1)then !fluid particles only
   
           j=ll(i)
           do while(j.gt.0) !loop in same cell
      if(p(j)%id.eq.1)then!fluid particles only for free-slip bc     
       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x**2+rij_y**2

        if(rij.le.kf*kf*h*h)then
       call ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)
       !kernal gradient correction
       dwdx_corr(i)=p(i)%mat(1,1)*dwdx+p(i)%mat(1,2)*dwdy
       dwdy_corr(i)=p(i)%mat(2,1)*dwdx+p(i)%mat(2,2)*dwdy
       dwdx_corr(j)=p(j)%mat(1,1)*dwdx+p(j)%mat(1,2)*dwdy
       dwdy_corr(j)=p(j)%mat(2,1)*dwdx+p(j)%mat(2,2)*dwdy

     if(visco_type.eq.1)then
              
     !laminar +sps turbulence viscosity
      
       var1=p(i)%vel%x(1)-p(j)%vel%x(1) 
       var2=p(i)%vel%x(2)-p(j)%vel%x(2)
       var3= rij_x*dwdx_corr(i)+ rij_y*dwdy_corr(i)
       var6= rij_x*dwdx_corr(j)+ rij_y*dwdy_corr(j)
       var4=(p(i)%dens+p(j)%dens)
       var5=rij+0.01*h**2
                       
    
     !laminar viscosity  
      f_viscx(i)=f_viscx(i)+((4.0*p(j)%mass*nu*rij_x*var3)/(var4*var5))*var1
      f_viscy(i)=f_viscy(i)+((4.0*p(j)%mass*nu*rij_y*var3)/(var4*var5))*var2
      f_viscx(j)=f_viscx(j)-((4.0*p(j)%mass*nu*rij_x*var6)/(var4*var5))*var1
      f_viscy(j)=f_viscy(j)-((4.0*p(j)%mass*nu*rij_y*var6)/(var4*var5))*var2

      !sps turbulence model
t_forcex(i)=t_forcex(i)+(p(j)%mass*(p(i)%txx/p(i)%dens**2+p(j)%txx/p(j)%dens**2)*dwdx_corr(i))&
                       +(p(j)%mass*(p(i)%txy/p(i)%dens**2+p(j)%txy/p(j)%dens**2)*dwdy_corr(i))

t_forcey(i)=t_forcey(i)+(p(j)%mass*(p(i)%txy/p(i)%dens**2+p(j)%txy/p(j)%dens**2)*dwdx_corr(i))&
                       +(p(j)%mass*(p(i)%tyy/p(i)%dens**2+p(j)%tyy/p(j)%dens**2)*dwdy_corr(i))

t_forcex(j)=t_forcex(j)-(p(i)%mass*(p(i)%txx/p(i)%dens**2+p(j)%txx/p(j)%dens**2)*dwdx_corr(j))&
                       -(p(i)%mass*(p(i)%txy/p(i)%dens**2+p(j)%txy/p(j)%dens**2)*dwdy_corr(j))

t_forcey(j)=t_forcey(j)-(p(i)%mass*(p(i)%txy/p(i)%dens**2+p(j)%txy/p(j)%dens**2)*dwdx_corr(j))&
                       -(p(i)%mass*(p(i)%tyy/p(i)%dens**2+p(j)%tyy/p(j)%dens**2)*dwdy_corr(j))
      

    end if

   if(visco_type.eq.2)then
    !artificial viscosity
      var1=(p(i)%vel%x(1)-p(j)%vel%x(1))*rij_x + (p(i)%vel%x(2)-p(j)%vel%x(2))*rij_y    
      var2=rij+(0.01*h)**2
      tij=var1/var2
         
        !artificial viscosity computation  
       sum1(i)=sum1(i)+((p(j)%mass/p(j)%dens)*tij*dwdx_corr(i))
       sum2(i)=sum2(i)+((p(j)%mass/p(j)%dens)*tij*dwdy_corr(i))
       sum1(j)=sum1(j)-((p(i)%mass/p(i)%dens)*tij*dwdx_corr(j))
       sum2(j)=sum2(j)-((p(i)%mass/p(i)%dens)*tij*dwdy_corr(j))
        end if


       if(visco_type.eq.3)then !orginal viscosity
                   
         !computing acc due to viscous stress
   
    fxx(i)=fxx(i)+p(i)%mat(2,1)*(p(j)%hxx*wij*(p(j)%mass/p(j)%dens))+p(i)%mat(2,2)*(p(j)%hxx*wij*(-rij_x)*(p(j)%mass/p(j)%dens))&
                 +p(i)%mat(2,3)*(p(j)%hxx*wij*(-rij_y)*(p(j)%mass/p(j)%dens))


    fxx(j)=fxx(j)+p(j)%mat(2,1)*(p(i)%hxx*wij*(p(i)%mass/p(i)%dens))+p(j)%mat(2,2)*(p(i)%hxx*wij*rij_x*(p(i)%mass/p(i)%dens))&
                 +p(j)%mat(2,3)*(p(i)%hxx*wij*rij_y*(p(i)%mass/p(i)%dens))

                                  

    fxy(i)=fxy(i)+p(i)%mat(3,1)*(p(j)%hxy*wij*(p(j)%mass/p(j)%dens))+p(i)%mat(3,2)*(p(j)%hxy*wij*(-rij_x)*(p(j)%mass/p(j)%dens))&
                 +p(i)%mat(3,3)*(p(j)%hxy*wij*(-rij_y)*(p(j)%mass/p(j)%dens))

    fxy(j)=fxy(j)+p(j)%mat(3,1)*(p(i)%hxy*wij*(p(j)%mass/p(j)%dens))+p(i)%mat(3,2)*(p(i)%hxy*wij*rij_x*(p(i)%mass/p(i)%dens))&
                 +p(j)%mat(3,3)*(p(i)%hxy*wij*rij_y*(p(i)%mass/p(i)%dens)) 

                 

    fyy(i)=fyy(i)+p(i)%mat(3,1)*(p(j)%hyy*wij*(p(j)%mass/p(j)%dens))+p(i)%mat(3,2)*(p(j)%hyy*wij*(-rij_x)*(p(j)%mass/p(j)%dens))&
                 +p(i)%mat(3,3)*(p(j)%hyy*wij*(-rij_y)*(p(j)%mass/p(j)%dens))

    fyy(j)=fyy(j)+p(j)%mat(3,1)*(p(i)%hyy*wij*(p(i)%mass/p(i)%dens))+p(j)%mat(3,2)*(p(i)%hyy*wij*rij_x*(p(i)%mass/p(i)%dens))&
                 +p(j)%mat(3,3)*(p(i)%hyy*wij*rij_y*(p(i)%mass/p(i)%dens))


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
     if(p(j)%id.eq.1)then!fluid particles only for free-slip bc
       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x**2+rij_y**2

        if(rij.le.kf*kf*h*h)then
         call ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)
         !kernal gradient correction
       dwdx_corr(i)=p(i)%mat(1,1)*dwdx+p(i)%mat(1,2)*dwdy
       dwdy_corr(i)=p(i)%mat(2,1)*dwdx+p(i)%mat(2,2)*dwdy
       dwdx_corr(j)=p(j)%mat(1,1)*dwdx+p(j)%mat(1,2)*dwdy
       dwdy_corr(j)=p(j)%mat(2,1)*dwdx+p(j)%mat(2,2)*dwdy
         
      if(visco_type.eq.1)then
           
     !laminar viscosity+sps turbulence
      
       var1=p(i)%vel%x(1)-p(j)%vel%x(1) 
       var2=p(i)%vel%x(2)-p(j)%vel%x(2)
       var3= rij_x*dwdx_corr(i)+ rij_y*dwdy_corr(i)
       var6= rij_x*dwdx_corr(j)+ rij_y*dwdy_corr(j)
       var4=(p(i)%dens+p(j)%dens)
       var5=rij+0.01*h**2
                       
    
      !laminar viscosity 
      f_viscx(i)=f_viscx(i)+((4.0*p(j)%mass*nu*rij_x*var3)/(var4*var5))*var1
      f_viscy(i)=f_viscy(i)+((4.0*p(j)%mass*nu*rij_y*var3)/(var4*var5))*var2
      f_viscx(j)=f_viscx(j)-((4.0*p(j)%mass*nu*rij_x*var6)/(var4*var5))*var1
      f_viscy(j)=f_viscy(j)-((4.0*p(j)%mass*nu*rij_y*var6)/(var4*var5))*var2






      !sps turbulence model
t_forcex(i)=t_forcex(i)+(p(j)%mass*(p(i)%txx/p(i)%dens**2+p(j)%txx/p(j)%dens**2)*dwdx_corr(i))&
                       +(p(j)%mass*(p(i)%txy/p(i)%dens**2+p(j)%txy/p(j)%dens**2)*dwdy_corr(i))

t_forcey(i)=t_forcey(i)+(p(j)%mass*(p(i)%txy/p(i)%dens**2+p(j)%txy/p(j)%dens**2)*dwdx_corr(i))&
                       +(p(j)%mass*(p(i)%tyy/p(i)%dens**2+p(j)%tyy/p(j)%dens**2)*dwdy_corr(i))

t_forcex(j)=t_forcex(j)-(p(i)%mass*(p(i)%txx/p(i)%dens**2+p(j)%txx/p(j)%dens**2)*dwdx_corr(j))&
                       -(p(i)%mass*(p(i)%txy/p(i)%dens**2+p(j)%txy/p(j)%dens**2)*dwdy_corr(j))

t_forcey(j)=t_forcey(j)-(p(i)%mass*(p(i)%txy/p(i)%dens**2+p(j)%txy/p(j)%dens**2)*dwdx_corr(j))&
                       -(p(i)%mass*(p(i)%tyy/p(i)%dens**2+p(j)%tyy/p(j)%dens**2)*dwdy_corr(j))
                       

     end if
  
       if(visco_type.eq.2)then
       
        !artificial viscosity computation  
      var1=(p(i)%vel%x(1)-p(j)%vel%x(1))*rij_x + (p(i)%vel%x(2)-p(j)%vel%x(2))*rij_y    
      var2=rij+(0.01*h)**2
      tij=var1/var2
           
       sum1(i)=sum1(i)+((p(j)%mass/p(j)%dens)*tij*dwdx_corr(i))
       sum2(i)=sum2(i)+((p(j)%mass/p(j)%dens)*tij*dwdy_corr(i))
       sum1(j)=sum1(j)-((p(i)%mass/p(i)%dens)*tij*dwdx_corr(j))
       sum2(j)=sum2(j)-((p(i)%mass/p(i)%dens)*tij*dwdy_corr(j))
      end if

      if(visco_type.eq.3)then !orginal viscosity
                   
         !computing acc due to viscous stress
   
    fxx(i)=fxx(i)+p(i)%mat(2,1)*(p(j)%hxx*wij*(p(j)%mass/p(j)%dens))+p(i)%mat(2,2)*(p(j)%hxx*wij*(-rij_x)*(p(j)%mass/p(j)%dens))&
                 +p(i)%mat(2,3)*(p(j)%hxx*wij*(-rij_y)*(p(j)%mass/p(j)%dens))


    fxx(j)=fxx(j)+p(j)%mat(2,1)*(p(i)%hxx*wij*(p(i)%mass/p(i)%dens))+p(j)%mat(2,2)*(p(i)%hxx*wij*rij_x*(p(i)%mass/p(i)%dens))&
                 +p(j)%mat(2,3)*(p(i)%hxx*wij*rij_y*(p(i)%mass/p(i)%dens))

                                  

    fxy(i)=fxy(i)+p(i)%mat(3,1)*(p(j)%hxy*wij*(p(j)%mass/p(j)%dens))+p(i)%mat(3,2)*(p(j)%hxy*wij*(-rij_x)*(p(j)%mass/p(j)%dens))&
                 +p(i)%mat(3,3)*(p(j)%hxy*wij*(-rij_y)*(p(j)%mass/p(j)%dens))

    fxy(j)=fxy(j)+p(j)%mat(3,1)*(p(i)%hxy*wij*(p(j)%mass/p(j)%dens))+p(i)%mat(3,2)*(p(i)%hxy*wij*rij_x*(p(i)%mass/p(i)%dens))&
                 +p(j)%mat(3,3)*(p(i)%hxy*wij*rij_y*(p(i)%mass/p(i)%dens)) 

                 

    fyy(i)=fyy(i)+p(i)%mat(3,1)*(p(j)%hyy*wij*(p(j)%mass/p(j)%dens))+p(i)%mat(3,2)*(p(j)%hyy*wij*(-rij_x)*(p(j)%mass/p(j)%dens))&
                 +p(i)%mat(3,3)*(p(j)%hyy*wij*(-rij_y)*(p(j)%mass/p(j)%dens))

    fyy(j)=fyy(j)+p(j)%mat(3,1)*(p(i)%hyy*wij*(p(i)%mass/p(i)%dens))+p(j)%mat(3,2)*(p(i)%hyy*wij*rij_x*(p(i)%mass/p(i)%dens))&
                 +p(j)%mat(3,3)*(p(i)%hyy*wij*rij_y*(p(i)%mass/p(i)%dens))



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

                     
if(visco_type.eq.1)then 
!laminar viscosity+sps turbulence
do i=1,ntotal
 p(i)%acc%x(1)=p(i)%acc%x(1)+f_viscx(i)+t_forcex(i)
 p(i)%acc%x(2)=p(i)%acc%x(2)+f_viscy(i)+t_forcey(i)
end do
end if

if(visco_type.eq.2)then 
!artificial viscosity
do i=1,ntotal
p(i)%acc%x(1)=p(i)%acc%x(1)+(alpha*h*c(i)*(rho/p(i)%dens)*sum1(i))
p(i)%acc%x(2)=p(i)%acc%x(2)+(alpha*h*c(i)*(rho/p(i)%dens)*sum2(i))
end do
end if

  if(visco_type.eq.3)then 
!orginal viscosity
do i=1,ntotal
p(i)%acc%x(1)=p(i)%acc%x(1)+(mu/p(i)%dens)*(fxx(i)+fxy(i))
p(i)%acc%x(2)=p(i)%acc%x(2)+(mu/p(i)%dens)*(fxy(i)+fyy(i))
end do
end if   
          
  end subroutine





 

subroutine viscosity_NoSlip(visco_type,kernal_type,p,ntotal,nu,alpha,coeff,h,dp,ep_vel,rho&
                            ,length,height,hswl)

implicit none
type(particle),intent(inout)::p(:)
double precision,intent(in)::h
integer,intent(in)::ntotal
double precision,intent(in)::length,height
double precision,intent(in)::dp
type(vector_2d),intent(in)::ep_vel(:) !extrapolated wall velocity
double precision,intent(in)::nu,hswl !kinematic viscosity
double precision,intent(in)::coeff
integer,intent(in)::visco_type !viscosity type laminar , artificial or orginal
integer,intent(in)::kernal_type

!local variables
integer::i,j,k,sc_index,i1,i2,neigh_cell
double precision::rij_x,rij_y,rij,wij,dwdx,dwdy
double precision,allocatable::dwdx_corr(:),dwdy_corr(:)
integer::cellx,celly,cell_total
double precision::xmax,ymax!max domain size
double precision::mu=1d-3 !dynamic viscosity of water
double precision,intent(in)::rho !reference fluid density
integer,allocatable::box(:),ll(:),np(:)
integer::ix,iy
integer,save::ncellx(4),ncelly(4) !no.of neighbour cells in x and y direction
double precision::var1,var2,var3,var4,var5,var6
double precision::tij
double precision,allocatable::sum1(:),sum2(:)
integer::j1,j2,jj
double precision,allocatable::f_viscx(:),f_viscy(:)
double precision,parameter::g=9.81
double precision,intent(in)::alpha
double precision,allocatable::c(:)
double precision::phi
double precision::kf
double precision,allocatable::t_forcex(:),t_forcey(:)
double precision,allocatable::fxx(:),fyy(:),fxy(:)

  
  

 
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



allocate(sum1(ntotal),sum2(ntotal))
allocate(box(cell_total),ll(ntotal),np(cell_total),f_viscx(ntotal),f_viscy(ntotal))
allocate(c(ntotal))
allocate(dwdx_corr(ntotal),dwdy_corr(ntotal))
allocate(t_forcex(ntotal),t_forcey(ntotal))
allocate(fxx(ntotal),fyy(ntotal),fxy(ntotal))

do i=1,ntotal
if(p(i)%id.eq.1)then!only for fluid particles
 c(i)=coeff*sqrt(g*hswl) 
else
c(i)=0.0
end if
end do

    !initialize cells
       
   box=0  
        
   do i=1,cell_total
   np(i)=0
    end do
    
  do i=1,ntotal
    ll(i)=0
   end do



!initialize wall velocity to extrapolated velocity
  do i=1,ntotal
    if(p(i)%id.gt.1)then
      p(i)%vel%x(1)=ep_vel(i)%x(1)
      p(i)%vel%x(2)=ep_vel(i)%x(2)
      end if
      end do




!initialize laminar-viscosity forces
do i=1,ntotal
 f_viscx(i)=0.0d0
 f_viscy(i)=0.0d0
sum1(i)=0.0
sum2(i)=0.0
dwdx_corr(i)=0.0d0
dwdy_corr(i)=0.0d0
end do

 DO i=1,ntotal
    fxx(i)=0.0d0
    fxy(i)=0.0d0
    fyy(i)=0.0d0
    end do



!initialize turbulence-viscosity forces
do i=1,ntotal
  t_forcex(i)=0.0d0
  t_forcey(i)=0.0d0
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
 
     !computation of acceleration due to viscous interaction using cell-linked list

  !call matrix(p,ntotal,h) 

      do i2=1,celly
          do i1=1,cellx
           sc_index=i1+(i2-1)*cellx
       if(box(sc_index).gt.0)then!there are particles in the cell
         i=box(sc_index)
                         
         do while(i.gt.0)

          if(p(i)%id.eq.1)then !fluid particles only
   
           j=ll(i)
           do while(j.gt.0) !loop in same cell
       !both fluid and boundary particles  
       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x**2+rij_y**2

        if(rij.le.kf*kf*h*h)then
       call ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)
       !kernal gradient correction
       dwdx_corr(i)=p(i)%mat(1,1)*dwdx+p(i)%mat(1,2)*dwdy
       dwdy_corr(i)=p(i)%mat(2,1)*dwdx+p(i)%mat(2,2)*dwdy
       dwdx_corr(j)=p(j)%mat(1,1)*dwdx+p(j)%mat(1,2)*dwdy
       dwdy_corr(j)=p(j)%mat(2,1)*dwdx+p(j)%mat(2,2)*dwdy

     if(visco_type.eq.1)then
              
     !laminar +sps turbulence viscosity
      
       var1=p(i)%vel%x(1)-p(j)%vel%x(1) 
       var2=p(i)%vel%x(2)-p(j)%vel%x(2)
       var3= rij_x*dwdx_corr(i)+ rij_y*dwdy_corr(i)
       var6= rij_x*dwdx_corr(j)+ rij_y*dwdy_corr(j)
       var4=(p(i)%dens+p(j)%dens)
       var5=rij+0.01*h**2
                       
    
     !laminar viscosity  
      f_viscx(i)=f_viscx(i)+((4.0*p(j)%mass*nu*rij_x*var3)/(var4*var5))*var1
      f_viscy(i)=f_viscy(i)+((4.0*p(j)%mass*nu*rij_y*var3)/(var4*var5))*var2
      f_viscx(j)=f_viscx(j)-((4.0*p(j)%mass*nu*rij_x*var6)/(var4*var5))*var1
      f_viscy(j)=f_viscy(j)-((4.0*p(j)%mass*nu*rij_y*var6)/(var4*var5))*var2

      !sps turbulence model
t_forcex(i)=t_forcex(i)+(p(j)%mass*(p(i)%txx/p(i)%dens**2+p(j)%txx/p(j)%dens**2)*dwdx_corr(i))&
                       +(p(j)%mass*(p(i)%txy/p(i)%dens**2+p(j)%txy/p(j)%dens**2)*dwdy_corr(i))

t_forcey(i)=t_forcey(i)+(p(j)%mass*(p(i)%txy/p(i)%dens**2+p(j)%txy/p(j)%dens**2)*dwdx_corr(i))&
                       +(p(j)%mass*(p(i)%tyy/p(i)%dens**2+p(j)%tyy/p(j)%dens**2)*dwdy_corr(i))

t_forcex(j)=t_forcex(j)-(p(i)%mass*(p(i)%txx/p(i)%dens**2+p(j)%txx/p(j)%dens**2)*dwdx_corr(j))&
                       -(p(i)%mass*(p(i)%txy/p(i)%dens**2+p(j)%txy/p(j)%dens**2)*dwdy_corr(j))

t_forcey(j)=t_forcey(j)-(p(i)%mass*(p(i)%txy/p(i)%dens**2+p(j)%txy/p(j)%dens**2)*dwdx_corr(j))&
                       -(p(i)%mass*(p(i)%tyy/p(i)%dens**2+p(j)%tyy/p(j)%dens**2)*dwdy_corr(j))
      

    end if

   if(visco_type.eq.2)then
    !artificial viscosity
      var1=(p(i)%vel%x(1)-p(j)%vel%x(1))*rij_x + (p(i)%vel%x(2)-p(j)%vel%x(2))*rij_y    
      var2=rij+(0.01*h)**2
      tij=var1/var2
         
        !artificial viscosity computation  
       sum1(i)=sum1(i)+((p(j)%mass/p(j)%dens)*tij*dwdx_corr(i))
       sum2(i)=sum2(i)+((p(j)%mass/p(j)%dens)*tij*dwdy_corr(i))
       sum1(j)=sum1(j)-((p(i)%mass/p(i)%dens)*tij*dwdx_corr(j))
       sum2(j)=sum2(j)-((p(i)%mass/p(i)%dens)*tij*dwdy_corr(j))
        end if
           
               if(visco_type.eq.3)then !orginal viscosity
                   
         !computing acc due to viscous stress
   
    fxx(i)=fxx(i)+p(i)%mat(2,1)*(p(j)%hxx*wij*(p(j)%mass/p(j)%dens))+p(i)%mat(2,2)*(p(j)%hxx*wij*(-rij_x)*(p(j)%mass/p(j)%dens))&
                 +p(i)%mat(2,3)*(p(j)%hxx*wij*(-rij_y)*(p(j)%mass/p(j)%dens))


    fxx(j)=fxx(j)+p(j)%mat(2,1)*(p(i)%hxx*wij*(p(i)%mass/p(i)%dens))+p(j)%mat(2,2)*(p(i)%hxx*wij*rij_x*(p(i)%mass/p(i)%dens))&
                 +p(j)%mat(2,3)*(p(i)%hxx*wij*rij_y*(p(i)%mass/p(i)%dens))

                                  

    fxy(i)=fxy(i)+p(i)%mat(3,1)*(p(j)%hxy*wij*(p(j)%mass/p(j)%dens))+p(i)%mat(3,2)*(p(j)%hxy*wij*(-rij_x)*(p(j)%mass/p(j)%dens))&
                 +p(i)%mat(3,3)*(p(j)%hxy*wij*(-rij_y)*(p(j)%mass/p(j)%dens))

    fxy(j)=fxy(j)+p(j)%mat(3,1)*(p(i)%hxy*wij*(p(j)%mass/p(j)%dens))+p(i)%mat(3,2)*(p(i)%hxy*wij*rij_x*(p(i)%mass/p(i)%dens))&
                 +p(j)%mat(3,3)*(p(i)%hxy*wij*rij_y*(p(i)%mass/p(i)%dens)) 

                 

    fyy(i)=fyy(i)+p(i)%mat(3,1)*(p(j)%hyy*wij*(p(j)%mass/p(j)%dens))+p(i)%mat(3,2)*(p(j)%hyy*wij*(-rij_x)*(p(j)%mass/p(j)%dens))&
                 +p(i)%mat(3,3)*(p(j)%hyy*wij*(-rij_y)*(p(j)%mass/p(j)%dens))

    fyy(j)=fyy(j)+p(j)%mat(3,1)*(p(i)%hyy*wij*(p(i)%mass/p(i)%dens))+p(j)%mat(3,2)*(p(i)%hyy*wij*rij_x*(p(i)%mass/p(i)%dens))&
                 +p(j)%mat(3,3)*(p(i)%hyy*wij*rij_y*(p(i)%mass/p(i)%dens))


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
     
       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x**2+rij_y**2

        if(rij.le.kf*kf*h*h)then
         call ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)
         !kernal gradient correction
       dwdx_corr(i)=p(i)%mat(1,1)*dwdx+p(i)%mat(1,2)*dwdy
       dwdy_corr(i)=p(i)%mat(2,1)*dwdx+p(i)%mat(2,2)*dwdy
       dwdx_corr(j)=p(j)%mat(1,1)*dwdx+p(j)%mat(1,2)*dwdy
       dwdy_corr(j)=p(j)%mat(2,1)*dwdx+p(j)%mat(2,2)*dwdy
         
      if(visco_type.eq.1)then
           
     !laminar viscosity+sps turbulence
      
       var1=p(i)%vel%x(1)-p(j)%vel%x(1) 
       var2=p(i)%vel%x(2)-p(j)%vel%x(2)
       var3= rij_x*dwdx_corr(i)+ rij_y*dwdy_corr(i)
       var6= rij_x*dwdx_corr(j)+ rij_y*dwdy_corr(j)
       var4=(p(i)%dens+p(j)%dens)
       var5=rij+0.01*h**2
                       
    
      !laminar viscosity 
      f_viscx(i)=f_viscx(i)+((4.0*p(j)%mass*nu*rij_x*var3)/(var4*var5))*var1
      f_viscy(i)=f_viscy(i)+((4.0*p(j)%mass*nu*rij_y*var3)/(var4*var5))*var2
      f_viscx(j)=f_viscx(j)-((4.0*p(j)%mass*nu*rij_x*var6)/(var4*var5))*var1
      f_viscy(j)=f_viscy(j)-((4.0*p(j)%mass*nu*rij_y*var6)/(var4*var5))*var2






      !sps turbulence model
t_forcex(i)=t_forcex(i)+(p(j)%mass*(p(i)%txx/p(i)%dens**2+p(j)%txx/p(j)%dens**2)*dwdx_corr(i))&
                       +(p(j)%mass*(p(i)%txy/p(i)%dens**2+p(j)%txy/p(j)%dens**2)*dwdy_corr(i))

t_forcey(i)=t_forcey(i)+(p(j)%mass*(p(i)%txy/p(i)%dens**2+p(j)%txy/p(j)%dens**2)*dwdx_corr(i))&
                       +(p(j)%mass*(p(i)%tyy/p(i)%dens**2+p(j)%tyy/p(j)%dens**2)*dwdy_corr(i))

t_forcex(j)=t_forcex(j)-(p(i)%mass*(p(i)%txx/p(i)%dens**2+p(j)%txx/p(j)%dens**2)*dwdx_corr(j))&
                       -(p(i)%mass*(p(i)%txy/p(i)%dens**2+p(j)%txy/p(j)%dens**2)*dwdy_corr(j))

t_forcey(j)=t_forcey(j)-(p(i)%mass*(p(i)%txy/p(i)%dens**2+p(j)%txy/p(j)%dens**2)*dwdx_corr(j))&
                       -(p(i)%mass*(p(i)%tyy/p(i)%dens**2+p(j)%tyy/p(j)%dens**2)*dwdy_corr(j))
                       

     end if
  
       if(visco_type.eq.2)then
       
        !artificial viscosity computation  
      var1=(p(i)%vel%x(1)-p(j)%vel%x(1))*rij_x + (p(i)%vel%x(2)-p(j)%vel%x(2))*rij_y    
      var2=rij+(0.01*h)**2
      tij=var1/var2
           
       sum1(i)=sum1(i)+((p(j)%mass/p(j)%dens)*tij*dwdx_corr(i))
       sum2(i)=sum2(i)+((p(j)%mass/p(j)%dens)*tij*dwdy_corr(i))
       sum1(j)=sum1(j)-((p(i)%mass/p(i)%dens)*tij*dwdx_corr(j))
       sum2(j)=sum2(j)-((p(i)%mass/p(i)%dens)*tij*dwdy_corr(j))
      end if

        if(visco_type.eq.3)then !orginal viscosity
                   
         !computing acc due to viscous stress
   
    fxx(i)=fxx(i)+p(i)%mat(2,1)*(p(j)%hxx*wij*(p(j)%mass/p(j)%dens))+p(i)%mat(2,2)*(p(j)%hxx*wij*(-rij_x)*(p(j)%mass/p(j)%dens))&
                 +p(i)%mat(2,3)*(p(j)%hxx*wij*(-rij_y)*(p(j)%mass/p(j)%dens))


    fxx(j)=fxx(j)+p(j)%mat(2,1)*(p(i)%hxx*wij*(p(i)%mass/p(i)%dens))+p(j)%mat(2,2)*(p(i)%hxx*wij*rij_x*(p(i)%mass/p(i)%dens))&
                 +p(j)%mat(2,3)*(p(i)%hxx*wij*rij_y*(p(i)%mass/p(i)%dens))

                                  

    fxy(i)=fxy(i)+p(i)%mat(3,1)*(p(j)%hxy*wij*(p(j)%mass/p(j)%dens))+p(i)%mat(3,2)*(p(j)%hxy*wij*(-rij_x)*(p(j)%mass/p(j)%dens))&
                 +p(i)%mat(3,3)*(p(j)%hxy*wij*(-rij_y)*(p(j)%mass/p(j)%dens))

    fxy(j)=fxy(j)+p(j)%mat(3,1)*(p(i)%hxy*wij*(p(j)%mass/p(j)%dens))+p(i)%mat(3,2)*(p(i)%hxy*wij*rij_x*(p(i)%mass/p(i)%dens))&
                 +p(j)%mat(3,3)*(p(i)%hxy*wij*rij_y*(p(i)%mass/p(i)%dens)) 

                 

    fyy(i)=fyy(i)+p(i)%mat(3,1)*(p(j)%hyy*wij*(p(j)%mass/p(j)%dens))+p(i)%mat(3,2)*(p(j)%hyy*wij*(-rij_x)*(p(j)%mass/p(j)%dens))&
                 +p(i)%mat(3,3)*(p(j)%hyy*wij*(-rij_y)*(p(j)%mass/p(j)%dens))

    fyy(j)=fyy(j)+p(j)%mat(3,1)*(p(i)%hyy*wij*(p(i)%mass/p(i)%dens))+p(j)%mat(3,2)*(p(i)%hyy*wij*rij_x*(p(i)%mass/p(i)%dens))&
                 +p(j)%mat(3,3)*(p(i)%hyy*wij*rij_y*(p(i)%mass/p(i)%dens))

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

                     
if(visco_type.eq.1)then 
!laminar viscosity+sps turbulence
do i=1,ntotal
 p(i)%acc%x(1)=p(i)%acc%x(1)+f_viscx(i)+t_forcex(i)
 p(i)%acc%x(2)=p(i)%acc%x(2)+f_viscy(i)+t_forcey(i)
end do
end if

if(visco_type.eq.2)then 
!artificial viscosity
do i=1,ntotal
p(i)%acc%x(1)=p(i)%acc%x(1)+(alpha*h*c(i)*(rho/p(i)%dens)*sum1(i))
p(i)%acc%x(2)=p(i)%acc%x(2)+(alpha*h*c(i)*(rho/p(i)%dens)*sum2(i))
end do
end if

 !orginal viscosity formulation
 if(visco_type.eq.3)then 
  do i=1,ntotal
  p(i)%acc%x(1)=p(i)%acc%x(1)+(mu/p(i)%dens)*(fxx(i)+fxy(i))
  p(i)%acc%x(2)=p(i)%acc%x(2)+(mu/p(i)%dens)*(fxy(i)+fyy(i))
  end do
    end if   


 end subroutine

 

   
   end module

   
    
    





