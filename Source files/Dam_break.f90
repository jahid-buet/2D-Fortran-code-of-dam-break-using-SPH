
program main
use dam
use particles
use ini
use properties
use kernal
use density
use pressure
use Compute_press_grad
use CompStress
use CompStrainRate
use CompGravity
use kernal_gradient_correction
use ptcl_matrix
use step
implicit none

double precision::x !length
double precision::y !width
double precision::z  !height
double precision::dp
double precision::ini_posx,ini_posy
double precision::fl_length!fluid length
double precision::fl_height !fluid depth
double precision::dp_z
double precision::dt !time step
double precision::coeffs !coeff of sound
double precision::nu!kinematic viscosity of water
double precision::a!artificial viscosity of fluid
double precision::rho0 !density of water
double precision::rhob!density of boundary
double precision::h,time,time1,time2
integer::nx,nz,ntotal,ntotal_f,nx_fl,nz_fl,i,j
type(particle),allocatable::p(:)
double precision,allocatable::drodt(:),drodt_old(:),dens_old(:)
integer::rank,ierr,ierror
integer::nt,t,out_nt,n1,run
double precision::out_time
double precision,allocatable::vel_halfx(:),vel_halfy(:),pos_half_x(:),pos_half_y(:),acc_oldx(:),acc_oldy(:)
double precision,allocatable::vel_oldx(:),vel_oldy(:),pos_oldx(:),pos_oldy(:),corr_vx(:),corr_vy(:)
type(vector_2d),allocatable::pres_vel(:)!prescribed wall velocity
type(vector_2d),allocatable::extra_vel(:)!extrapolated wall velocity
double precision,allocatable::v1(:)
double precision::simtime!simulation  time
double precision::var1,dp2,eta
integer::bound_type    !type 1=slip-boundary condition
                       !type 2=no-slip boundary condition


integer::type              !type 1=laminar +sps turbulnce viscosity formulation                                     
                           !type 2=artificial viscosity formulation
                           !type 3=orginal viscosity formulation

integer::k_type     !type 1=quintic spline kernal
                      !type 2=wendland kernal

integer,parameter::time_integrator=1   !1=predictor-corrector time integration
                                      !2=velocity-verlet time integration
character::err_msg
integer::t1,t2,clock_rate,clock_max
real::a1,b1,c1
integer::day,hr,mint
real::sec
integer::status=0

!do while(run.lt.2)

write(*,140)"This is a simple code of  2d-dam break problem using SPH numerical method developed by Md.Jahid Hasan.Water resources engineering,M.S.C,BUET id-0419162024"
140 format(/,2x,a)
write(*,141)"##########################################################################################################################################################"
141 format(2x,a)
write(*,200)'Note:If you encounter any bug/error feel free to notify me or you can modify the source files.'
200 format(2x,a123)
write(*,201)"***********************************************************************************************"
201 format(2x,a123)
write(*,142)"All values should be given in SI units"
142 format(2x,a80)
write(*,143)"*************************************************************"
143 format(15x,a80)
write(*,160)"spacing between particles(dp)=?"
160 format(/,1x,a32)
write(*,166)"Note:dp should be choosen in such a way that ratio of length of wall in x or y direction and dp is an integer"
166 format(1x,a110)
read(*,*) dp
write(*,144)"Its assuming that all particles in domain must have positive coordinate."
144 format(1x,a82)
write(*,170)"***********************************************************************"
170 format(1x,a80)
write(*,171)"length of confined-wall in x-direction=?"
171 format(/,1x,a)
read(*,*) x
write(*,*)"length of confined-wall in y-direction=?"
read(*,*) z
write(*,*)"density of fluid particles=?"
read(*,*) rho0
write(*,133)"density of wall/boundary particles=?"
133 format(/,1x,a)
write(*,145)"Note:its suggessed that density of wall particles should be equal to that of fluid particles."
145 format(1x,a93)
write(*,146)"higher density cause higher repulsive force excerted by wall to fluid partlces which results in gap betn wall and fluid."
146 format(1x,a120)
read(*,*) rhob
write(*,134)"initial x-position of fluid column=?"
134 format(/,1x,a)
read(*,*) ini_posx
write(*,135)"initial y-position of fluid column=?"
135 format(/,1x,a)
read(*,*) ini_posy

 if(ini_posx.ge.x.or.ini_posy.ge.z)then
  write(*,*)"Error:fluid must contains within domain"  
  pause("press enter to exit....")
  stop
  end if


write(*,136)"width of fluid column=?"
136 format(/,1x,a)
read(*,*) fl_length
write(*,137)"height of fluid column=?"
137 format(/,1x,a)
read(*,*) fl_height
write(*,154)"type of kernal function that will be used =?(1=quintic spline kernal,2=wendland kernal)"
154 format(/,1x,a)
read(*,*) k_type
write(*,155)"coefficient to calculate artificial sound speed in fluid=?(10 to 20 gives good approximation of sound)"
155 format(/,1x,a)
!write(*,158)"Note:Higher value of coefficent results in small time-step size which further increases no.of time-step."
!158 format(1x,a)
read(*,*) coeffs
write(*,156)"which type of viscosity will be used=?(1=laminar+sps turbulence,2=artificial viscosity,3=orginal viscosity)"
156 format(/,1x,a)
read(*,*) type
if(type.eq.1.or.type.eq.3)then
  write(*,*)"kinematic viscosity of fluid=?"
  read(*,*) nu
elseif(type.eq.2)then
 write(*,157)"artificial viscosity coefficient(alpha)=?(0.01 to 0.04 give good approximation of fluid viscosity)"
 157 format(1x,a98)
 read(*,*) a
 end if
 
 write(*,169)"Note:wall boundary is modelled by using dummy particles(paper:A generalized wall boundary condition for smoothed particle hydrodynamics by Adami 2012)"
 169 format(/,1x,a)
write(*,138)"type of boundary condion=?(1-slip,2-no-slip)"
138 format(/,1x,a)
write(*,203)"Note:Its better to use slip bc when using artificial viscosity,otherwise particles sticking to wall/boundary can occur."
203 format(1x,a)
write(*,204)"In case of orginal viscosity as kgf-sph formulation is used,no-slip boundary condition causes program to crash after some timesteps(may be its related to  matrix formulation for boundary particles,and this problm is still under investigation)"
write(*,205)"so its advisable to use slip boundary condition when using orginal viscosity formulation"
 204 format(1x,a)
 205 format(1x,a)
read(*,*) bound_type


if(k_type.eq.1)then
h=1.0*dp
else
  h=1.5*dp
  end if

dp2=dp !cutoff radius for boundary particle creation

nx=nint(x/dp2)
nz=nint(z/dp2)
nx_fl=nint(fl_length/dp)
nz_fl=nint(fl_height/dp)
ntotal_f=nx_fl*nz_fl
ntotal=(4*nx+8*nz-22)+ntotal_f
!write(*,*)"ntotal=",ntotal
allocate(p(ntotal))
allocate(vel_halfx(ntotal),vel_halfy(ntotal),pos_half_x(ntotal),pos_half_y(ntotal),stat=status,errmsg=err_msg)
allocate(vel_oldx(ntotal),vel_oldy(ntotal),pos_oldx(ntotal),pos_oldy(ntotal),corr_vx(ntotal),corr_vy(ntotal),stat=status,errmsg=err_msg)
allocate(drodt(ntotal),drodt_old(ntotal),dens_old(ntotal),acc_oldx(ntotal),acc_oldy(ntotal),stat=status,errmsg=err_msg)
allocate(pres_vel(ntotal),extra_vel(ntotal),v1(ntotal))



write(*,147) "Generating particles coordinate(fluid,wall boundary.etc)...................................."
147 format(1x,a100)


call system_clock(t1,clock_rate,clock_max)
!creation of problem domain
 call create_domain(p,ntotal,x,ini_posx,fl_length,z,ini_posy,fl_height,dp,dp2)

write(*,148)"done"
148 format(50x,a5)
write(*,180)'*****************************************************'
180 format(1x,a80)
write(*,187)"total no. of fluid particles=",ntotal_f
187 format(/,1x,a,i7)
write(*,188)"total no. of boundary particles=",(ntotal-ntotal_f)
188 format(1x,a,i7)
write(*,191)"total no. of particles=",ntotal
191 format(1x,a,i7)


write(*,150)"Assigning properties(mass,density) and initializing variables(velocity,pressure,etc)to all particles in the domain......................"
150 format(/,2x,a120)

!assign properties(mass,density) to each of particles in domain
 call PROP(p,ntotal,rho0,rhob,dp,dp2)
 
 !initialize velocity,acceleration for all particles in the domain
 call initialize(p,ntotal)
 
write(*,182)"writing out initial particles position,density and pressures to a text file..............................."
182 format(/,1x,a80)

!save geometry to txt file
open(unit=10,file='coordinate_ini.txt',action='write'&
   ,status='replace',iostat=ierror)
     
    write(10,102) 'coordx','coordy','coordz',"density","pressure"
    102 format(1x,a12,2x,a12,2x,a12,2x,a12,2x,a12)
   do i=1,ntotal
    write(10,105) p(i)%coord%x(1),p(i)%coord%x(2),p(i)%coord%x(3),p(i)%dens,p(i)%press
      
    105 format(1x,f12.5,2x,f12.5,2x,f12.5,2x,f12.5,2x,f12.5)
  end do


write(*,151)"done"
151 format(50x,a5)
write(*,183)"**********************************************************"
183 format(1x,a80)

!time integration should be called here
write(*,*)"Simulation time(s)=?"
read(*,*) simtime
write(*,206)"At which time(s)output results(coordinates,density and pressure) will be saved=?"
206 format(1x,a)

read(*,*) out_time
 write(*,*)"Timestep size(dt=?)""(1e-5 to 5e-5 hopefully give stable simulation)"

  read(*,*) dt
  !call timestep routine
 
  !call timestep(1,p,ntotal,h,fl_height,coeffs,nu,dt,v1) !still buggy :()
  write(*,130)"Initial timestep size=",dt,"sec"
  130 format(1x,a28,ES12.4,a5)
  
  write(*,*)"*****************************************************"
  nt=int(simtime/dt)
  write(*,152)"Total no. of timestep=",nt  
  152 format(1x,a30,1x,i8)
  write(*,*)"****************************************************"
  out_nt=int(out_time/dt)

 
!For debugging purpose

!****************************************************************************
!$$$$$$  do i=1,5
!$$$$$$  call matrix(p,ntotal,h)
!$$$$$$  !call Strain_Rate(p,ntotal,h)
!$$$$$$  call gravity(p,ntotal)
!$$$$$$  !call timestep(i,p,ntotal,coeffs,h,fl_height,nu,dt)
!$$$$$$  call Compute_density_fluid(p,ntotal,h,drodt)
!$$$$$$  !call compute_press(p,ntotal,h,coeffs,fl_height)
!$$$$$$   !call Press_Gradient(p,ntotal,h)
!$$$$$$   !write(*,*)"dt=",dt
!$$$$$$   
!$$$$$$  
!$$$$$$ end do
!$$$$$$   do i=1,ntotal
!$$$$$$   write(*,*)"p%mat(1,1)=",p(i)%mat(1,1)
!$$$$$$   end do
!$$$$$$             do i=1,ntotal
!$$$$$$             write(*,*)"p(i)%drodt=",drodt(i)
!$$$$$$             end do
   !write(*,*)"p%dens=",p%dens
   !write(*,*)"p%acc=",p(:)%acc%x(3)
!****************************************************************************


  !initialize variables at time=0
  do i=1,ntotal
   if(p(i)%id.eq.1)then
   vel_oldx(i)=0.0d0
   vel_oldy(i)=0.0d0
   pos_oldx(i)=p(i)%coord%x(1)
   pos_oldy(i)=p(i)%coord%x(2)
   drodt(i)=0.0d0
   drodt_old(i)=0.0
   dens_old(i)=p(i)%dens!initialize old_density of fluid particles
   corr_vx(i)=0.0
   corr_vy(i)=0.0 
   acc_oldx(i)=0.0d0
   acc_oldy(i)=0.0d0
   vel_halfx(i)=0.0d0
   vel_halfy(i)=0.0d0
    end if
   end do
!initialize prescribed velocity
  do i=1,ntotal
   pres_vel(i)%x(1)=0.0d0
   pres_vel(i)%x(2)=0.0d0
   end do
write(*,*)"Starting time integration............................"
write(*,153)"predictor-corrector time integration is used."
153 format(3x,a50)
write(*,*)"****************************************************************"





  write(*,110)"Timestep" ,"Time(sec)"             
  110 format(4x,a12,2x,a12)
  write(*,111)"========","=========="
  111 format(4x,a12,2x,a12)
  time=0.0
  time1=0.0
  time2=0.0
  n1=0

!start time integration
  do t=1,nt
    if(t.eq.1)then
   time=time+dt
   end if    
     write(*,120) t,time
     120 format(1x,i12,4x,f12.6)
 

 if(time_integrator.eq.1)then                           
  !predictor-corrector time integration
 
!predictor step

!compute acceleration of fluid particles at n timestep
  if(type.eq.3)then
  call matrix(type,k_type,x,z,p,ntotal,h)
  end if
  call compute_press_fluid(p,ntotal,rho0,coeffs,fl_height) !compute pressure between fluid  particles
  call Compute_press_and_vel_bound(p,ntotal,h,dp,pres_vel,extra_vel,x,z,k_type) !compute pressure and velocity of boundary particles by extrapolation
  call Compute_density_bound(p,ntotal,rho0,coeffs,fl_height)!compute density of boundary particles
  call KGC(p,ntotal,h,dp,x,z,k_type) !kernal gradient correction
  call gravity(time,p,ntotal)
  call Press_Gradient(p,ntotal,h,dp,x,z,k_type)
  
     if(bound_type.eq.1)then
       if(type.eq.1.or.type.eq.3)then
         call  Strain_Rate_Slip(k_type,type,p,ntotal,h,dp,x,z)
      end if
       end if
    
       if(bound_type.gt.1)then
         if(type.eq.1.or.type.eq.3)then          
       call Strain_Rate_NoSlip(k_type,type,p,ntotal,h,dp,extra_vel,x,z)
         end if
            end if
     
   if(bound_type.eq.1)then
     if(type.eq.1.or.type.eq.3)then
   call viscosity_Slip(type,k_type,p,ntotal,nu,a,coeffs,h,dp,rho0,x,z,fl_height)
   end if
    end if
    
   if(bound_type.gt.1)then
     if(type.eq.1.or.type.eq.3)then
     call viscosity_NoSlip(type,k_type,p,ntotal,nu,a,coeffs,h,dp,extra_vel,rho0,x,z,fl_height)
     end if
      end if

 
  
 !compute drodt at n time level    
  call Compute_density_fluid(p,ntotal,h,dp,pres_vel,x,z,k_type,coeffs,fl_height,drodt,v1)
  
  !call timestep(t,p,ntotal,h,fl_height,coeffs,nu,dt,v1) 
  time=time+dt

     

!update density at half time step(n+1/2)
   do i=1,ntotal
     if(p(i)%id.eq.1)then ! density of fluid only
        p(i)%dens=dens_old(i)+(dt/2.0)*drodt(i)
     end if
   end do


!calculate velocity at half timestep(n+1/2)
  do i=1,ntotal
       if(p(i)%id.eq.1)then !onlyfluid particles
   
    p(i)%vel%x(1)=vel_oldx(i)+(dt/2.0)*p(i)%acc%x(1)
    p(i)%vel%x(2)=vel_oldy(i)+(dt/2.0)*p(i)%acc%x(2)                      
       end if
    end do


 
!update posotion at half timestep(n+1/2)
 do i=1,ntotal
   if(p(i)%id.eq.1)then !onlyfluid particles
     p(i)%coord%x(1)=pos_oldx(i)+(dt/2.0)*vel_oldx(i)
     p(i)%coord%x(2)=pos_oldy(i)+(dt/2.0)*vel_oldy(i)
     !p(i)%coord%x(1)=pos_oldx(i)+(dt/2.0)*corr_vx(i)
     !p(i)%coord%x(2)=pos_oldy(i)+(dt/2.0)*corr_vy(i)
     
  end if
    end do


!resetting particles acc to zero
   do i=1,ntotal
     p(i)%acc%x(1)=0.0
     p(i)%acc%x(2)=0.0
    end do


!corrector step
!compute acceleration of fluid particles at half timestep(n+1/2)
  if(type.eq.3)then
  call matrix(type,k_type,x,z,p,ntotal,h)
  end if
  call compute_press_fluid(p,ntotal,rho0,coeffs,fl_height) !compute pressure between fluid  particles
  call Compute_press_and_vel_bound(p,ntotal,h,dp,pres_vel,extra_vel,x,z,k_type) !compute pressure and velocity of boundary particles by extrapolation
  call Compute_density_bound(p,ntotal,rho0,coeffs,fl_height)!compute density of boundary particles
  call KGC(p,ntotal,h,dp,x,z,k_type) !kernal gradient correction
  call gravity(time,p,ntotal)
  call Press_Gradient(p,ntotal,h,dp,x,z,k_type)
  
    if(bound_type.eq.1)then
       if(type.eq.1.or.type.eq.3)then
         call  Strain_Rate_Slip(k_type,type,p,ntotal,h,dp,x,z)
      end if
       end if
    
       if(bound_type.gt.1)then
         if(type.eq.1.or.type.eq.3)then          
       call Strain_Rate_NoSlip(k_type,type,p,ntotal,h,dp,extra_vel,x,z)
         end if
            end if
     
   if(bound_type.eq.1)then
     if(type.eq.1.or.type.eq.3)then
   call viscosity_Slip(type,k_type,p,ntotal,nu,a,coeffs,h,dp,rho0,x,z,fl_height)
   end if
    end if
    
   if(bound_type.gt.1)then
     if(type.eq.1.or.type.eq.3)then
     call viscosity_NoSlip(type,k_type,p,ntotal,nu,a,coeffs,h,dp,extra_vel,rho0,x,z,fl_height)
     end if
      end if


   

!compute corrected  drodt at half timestep(n+1/2)
call Compute_density_fluid(p,ntotal,h,dp,pres_vel,x,z,k_type,coeffs,fl_height,drodt,v1)


  
  
  
    
!corrrect density of fluid particles at n+1/2 time level
 do i=1,ntotal
 if(p(i)%id.eq.1)then 
   p(i)%dens=dens_old(i)+(dt/2.0)*drodt(i)
 end if
  end do
 
!correct velocity at n+1/2 timestep
  do i=1,ntotal
   if(p(i)%id.eq.1)then !onlyfluid particles
   
    p(i)%vel%x(1)=vel_oldx(i)+(dt/2.0)*p(i)%acc%x(1)
    p(i)%vel%x(2)=vel_oldy(i)+(dt/2.0)*p(i)%acc%x(2)                      
       end if
    end do



  !XSph correction
 !call xsph_corr(p,ntotal,h,corr_vx,corr_vy)
 
!advance posotion at n+1/2 timestep
do i=1,ntotal
   if(p(i)%id.eq.1)then !onlyfluid particles
    p(i)%coord%x(1)= pos_oldx(i)+(dt/2.0)*p(i)%vel%x(1)
    p(i)%coord%x(2)= pos_oldy(i)+(dt/2.0)*p(i)%vel%x(2)

     !p(i)%coord%x(1)= pos_oldx(i)+(dt/2.0)*corr_vx(i)
     !p(i)%coord%x(2)= pos_oldy(i)+(dt/2.0)*corr_vy(i)

  end if
    end do


!advance density,position and velocity for next time step(n+1)

do i=1,ntotal
 if(p(i)%id.eq.1)then !onlyfluid particles
p(i)%dens=2.0*p(i)%dens-dens_old(i)
p(i)%vel%x(1)=2.0*p(i)%vel%x(1)-vel_oldx(i)
p(i)%vel%x(2)=2.0*p(i)%vel%x(2)-vel_oldy(i)
p(i)%coord%x(1)=2.0*p(i)%coord%x(1)-pos_oldx(i)
p(i)%coord%x(2)=2.0*p(i)%coord%x(2)-pos_oldy(i)
end if
end do

!checking if particles go out of domain if then these particles will be excluded from domain
do i=1,ntotal
if(p(i)%id.eq.1)then
  if(p(i)%coord%x(2).le.2.0*dp)then
     n1=n1+1    
       if(n1.lt.7)then   
     write(*,290)"warning*******"
  290 format(1x,a)
  write(*,291)"particle=",i,"tries to escape minimum domain limit in y direction which is=",0
  291 format(1x,a8,i6,a,f5.3)
  write(*,292)"particles=",i,"is excluded from domain"
  292 format(1x,a8,i6,a25)
   p(i)%coord%x(1)=0.0
   p(i)%coord%x(2)=0.0    
   p(i)%id=-1
  end if
    end if
  if(p(i)%coord%x(1).le.2.0*dp)then 
    p(i)%coord%x(1)=0.0
    p(i)%coord%x(2)=0.0  
   p(i)%id=-1
  end if
  if(p(i)%coord%x(1).gt.(x-3.0*dp))then
   p(i)%coord%x(1)=0.0
   p(i)%coord%x(2)=0.0 
   p(i)%id=-1
  end if
 
   end if
   
  end do


    if(n1.ge.7)then
      write(*,199)"Error:Too much particles try to escape the domain limit in minmum y coordinate"
      199 format(1x,a)
      write(*,190)"Try to reduce timestep or check if initial geometry is ok or not."
      190 format(1x,a)
      write(*,193)"Exiting program....................."
      193 format(1x,a)
      write(*,194)"-_-"
      194 format(1x,a)
      pause(" ")
      stop      
       end if 








!Dynamically changing the cell size in upper z limit
  z=maxval(p(:)%coord%x(2))
    




  
           

!XSph correction
   !call xsph_corr(p,ntotal,h,corr_vx,corr_vy)
!extrapolate correct velocity of boundary particles from fluid  for no-slip boundary condition at n+1 timestep

   !call No_slip_boundary_condition(p,ntotal,h)


 do i=1,ntotal
  if(p(i)%id.eq.1)then
   vel_oldx(i)=p(i)%vel%x(1)
   vel_oldy(i)=p(i)%vel%x(2)
   pos_oldx(i)=p(i)%coord%x(1)
   pos_oldy(i)=p(i)%coord%x(2)
   dens_old(i)=p(i)%dens
  end if
!dens_old(i)=p(i)%dens
 end do

!resetting particles acc to zero
   do i=1,ntotal
     p(i)%acc%x(1)=0.0
     p(i)%acc%x(2)=0.0
    end do


    if(mod(t,out_nt).eq.0)then
      write(*,*)"writing out final position,density and pressure of particles in a text file......."
      open(unit=11,file='coordinate_final.txt',action='write'&
   ,status='replace',iostat=ierror)
     
    write(11,101) 'coordx','coordy','coordz','density','pressure'
    101 format(1x,a12,2x,a12,2x,a12,2x,a12,2x,a12)
   do i=1,ntotal
    write(11,100) p(i)%coord%x(1),p(i)%coord%x(2),p(i)%coord%x(3),p(i)%dens,p(i)%press
      
    100 format(1x,f12.5,2x,f12.5,2x,f12.5,2x,f12.5,2x,f12.5)
  end do
   
 write(*,161)"done"
 161 format(15x,a)
 write(*,162)"**********************"
 162 format(6x,a)
   end if

 !time=time1+time2


      end if










  !velocity verlet algorithm(implementation of this method is buggy and still under investigation....)
  
  if(time_integrator.gt.1)then !velocity verlet algorithm
    !compute velocity and position at half timestep
    do i=1,ntotal
    if(p(i)%id.eq.1)then
      vel_halfx(i)=p(i)%vel%x(1)+(dt/2.0)*acc_oldx(i)
      vel_halfy(i)=p(i)%vel%x(2)+(dt/2.0)*acc_oldy(i)
      
     
      p(i)%coord%x(1)=p(i)%coord%x(1)+(dt/2.0)*vel_halfx(i)
      p(i)%coord%x(2)=p(i)%coord%x(2)+(dt/2.0)*vel_halfy(i)
     end if
       end do

   !compute dro/dt at n timestep
    call Compute_density_fluid(p,ntotal,h,dp,pres_vel,x,z,k_type,coeffs,fl_height,drodt,v1)

     !compute density of fluid  at half timestep
     do i=1,ntotal
       if(p(i)%id.eq.1)then
       p(i)%dens=dens_old(i)+(dt/2.0)*drodt(i)
       end if
       end do
 
   !compute dro/dt at halftime step(n+1/2)
      call Compute_density_fluid(p,ntotal,h,dp,pres_vel,x,z,k_type,coeffs,fl_height,drodt,v1)
      
       !compute density of fluid  at next timestep(n+1)
         do i=1,ntotal
       if(p(i)%id.eq.1)then
       p(i)%dens=dens_old(i)+(dt/2.0)*drodt(i)
       end if
       end do

    !compute position at next timestep(n+1)
     do i=1,ntotal
    if(p(i)%id.eq.1)then           
      p(i)%coord%x(1)=p(i)%coord%x(1)+(dt/2.0)*p(i)%vel%x(1)
      p(i)%coord%x(2)=p(i)%coord%x(2)+(dt/2.0)*p(i)%vel%x(2)
     end if
       end do

  !compute velocity at next timestep(n+1)
      
        do i=1,ntotal
         if(p(i)%id.eq.1)then
      p(i)%vel%x(1)=vel_oldx(i) +dt*acc_oldx(i)
      p(i)%vel%x(2)=vel_oldy(i) +dt*acc_oldy(i)
      
          end if
            end do




   !compute acceleration of fluid particles for next time step(n+1)
  call compute_press_fluid(p,ntotal,rho0,coeffs,fl_height) !compute pressure between fluid  particles
  call Compute_press_and_vel_bound(p,ntotal,h,dp,pres_vel,extra_vel,x,z,k_type) !compute pressures of boundary particles
  call Compute_density_bound(p,ntotal,rho0,coeffs,fl_height)!compute density of boundary particles
  call KGC(p,ntotal,h,dp,x,z,k_type) !kernal gradient correction
  call gravity(time,p,ntotal)
  call Press_Gradient(p,ntotal,h,dp,x,z,k_type)
  
     if(bound_type.eq.1)then
       if(type.eq.1.or.type.eq.3)then
         call  Strain_Rate_Slip(k_type,type,p,ntotal,h,dp,x,z)
      end if
       end if
    
       if(bound_type.gt.1)then
         if(type.eq.1.or.type.eq.3)then          
       call Strain_Rate_NoSlip(k_type,type,p,ntotal,h,dp,extra_vel,x,z)
         end if
            end if
     
   if(bound_type.eq.1)then
     if(type.eq.1.or.type.eq.3)then
   call viscosity_Slip(type,k_type,p,ntotal,nu,a,coeffs,h,dp,rho0,x,z,fl_height)
   end if
    end if
    
   if(bound_type.gt.1)then
     if(type.eq.1.or.type.eq.3)then
     call viscosity_NoSlip(type,k_type,p,ntotal,nu,a,coeffs,h,dp,extra_vel,rho0,x,z,fl_height)
     end if
      end if

     !correct velocity at next timestep(n+1)
      
         do i=1,ntotal
       if(p(i)%id.eq.1)then
        p(i)%vel%x(1)=vel_halfx(i)+(dt/2.0)*p(i)%acc%x(1)
        p(i)%vel%x(2)=vel_halfy(i)+(dt/2.0)*p(i)%acc%x(2)                 
       end if
         end do


 !checking if particles go out of domain if then resetting velocities of these particles to zero
    do i=1,ntotal
   if(p(i)%id.eq.1)then
  if(p(i)%coord%x(2).lt.2.0*dp)then
   p(i)%vel%x(1)=0.0d0
   p(i)%vel%x(2)=0.0
   p(i)%press=0.0
    end if
  if(p(i)%coord%x(1).gt.(x-2.0*dp))then
   p(i)%vel%x(1)=0.0d0
   p(i)%vel%x(2)=0.0
   p(i)%press=0.0
   end if
   end if
  end do




    do i=1,ntotal
    if(p(i)%id.eq.1)then
      vel_oldx(i)=p(i)%vel%x(1)
      vel_oldy(i)=p(i)%vel%x(2)
      acc_oldx(i)=p(i)%acc%x(1)
      acc_oldy(i)=p(i)%acc%x(2)   
      dens_old(i)=p(i)%dens
     end if
      end do

   !resetting particles acc to zero
    do i=1,ntotal
     p(i)%acc%x(1)=0.0
     p(i)%acc%x(2)=0.0
    end do

  
        end if
                               
        
           end do

   
   

   
  call system_clock(t2,clock_rate,clock_max)
  var1=(real(t2-t1)/real(clock_rate))
  if(var1.le.60.0d0)then
  write(*,113)var1,"sec"
  113 format("time need to finish the simulation =",1x,f8.5,a4)
  end if
  if(var1.gt.60.0d0 .and. var1.lt.3600.0d0)then
    a1=dble(var1/60.0d0)
    mint=int(var1/60.0d0)
    sec=dble((a1-mint)*60.0)
    write(*,114)"Time need to finish the simulation =",mint,"min",sec,"sec"
    114 FORMAT(1x,a,i3,a5,f6.2,a4)
    end if
   if(var1.gt.3600.0d0 .and. var1.le.86400.0d0)then
     a1=dble(var1/3600.0d0)
     hr=int(var1/3600.0d0)
     b1=dble((a1-hr)*60.0)
     mint=int((a1-hr)*60.0)
     sec=dble((b1-mint)*60.0)
     write(*,115)"Time needs to finish the simulation =",hr,"hr",mint,"min",sec,"sec"
     115 format(1x,a,i3,a4,i3,a4,f6.2,a4)
    end if
   if(var1.gt.86400.0d0)then
     a1=dble(var1/86400.0d0)
     day=int(var1/86400.0d0)
     b1=dble((a1-day)*24)
     hr=int((a1-day)*24)
     mint=int((b1-hr)*60.0)
     c1=dble((b1-hr)*60.0)
     sec=dble((c1-mint)*60.0)     
     write(*,116)"Time need to finish the simulation =",day,"day",hr,"hr",mint,"min",sec,"sec"                 
     116 format(1x,a,i2,a5,i3,a5,i3,a5,f6.2,a4)
    end if   
 
 pause("Press any key button to exit")
 
end program
