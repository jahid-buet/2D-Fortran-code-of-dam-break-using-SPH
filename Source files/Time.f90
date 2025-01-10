module step
contains
subroutine timestep(curtime,p,ntotal,h,hswl,coefsound,nu,dt,v1)
use particles
type(particle),intent(in)::p(:)
double precision,intent(in)::v1(:)
double precision,intent(out)::dt
double precision,intent(in)::h,nu,coefsound,hswl
integer,intent(in)::ntotal,curtime
!local variable
double precision::dt_force,dt_visc,cs,dt_cfl
double precision,allocatable::dt_cv(:)      
double precision,allocatable::var(:),velocity(:)
double precision::vmax

double precision::g=9.81
double precision::f_mag
integer::i

allocate(var(ntotal),dt_cv(ntotal),velocity(ntotal))

var=0.0
cs=coefsound*sqrt(g*hswl)
vmax=sqrt(g*hswl)
  
     if(curtime.eq.1)then
      dt=0.15*(h/(cs+vmax))
     end if


  if(curtime.gt.1)then
!$$$$$$     do i=1,ntotal   
!$$$$$$     if(p(i)%id.eq.1)then
!$$$$$$     dt_cv(i)=h/(cs+v1(i))    
!$$$$$$      if (dt_cv(i).gt.0.0)then
!$$$$$$       dt_cfl=dt_cv(i)
!$$$$$$       dt_cfl=min(dt_cv(i),dt_cfl)     
!$$$$$$    
!$$$$$$    end if
!$$$$$$     end if
!$$$$$$      end do
  
  dt_cfl=0.15*(h/(cs+vmax))  
  
 
   do i=1,ntotal
    if(p(i)%id.eq.1)then
   f_mag=sqrt((p(i)%acc%x(1))**2+(p(i)%acc%x(2))**2)
    var(i)=0.15*sqrt(h/(f_mag))
    if (var(i).gt.0.0)then
      dt_force=var(i)
      dt_force=min(var(i),dt_force)     
   
   end if
  end if
  end do


!$$$$$$    do i=1,ntotal
!$$$$$$     if(p(i)%id.eq.1)then
!$$$$$$    f_mag=sqrt((p(i)%acc%x(1))**2+(p(i)%acc%x(2))**2)
!$$$$$$    var(i)=0.25*sqrt(h/(f_mag))
!$$$$$$    end if
!$$$$$$    end do
!$$$$$$   dt_force=maxval(var)
  
 
!$$$$$$  do i=1,ntotal
!$$$$$$    if(p(i)%id.eq.1)then
!$$$$$$    velocity(i)=sqrt(p(i)%vel%x(1)**2+p(i)%vel%x(2)**2)
!$$$$$$   end if
!$$$$$$    end do
!$$$$$$   vmax=maxval(velocity)
!$$$$$$   dt_cfl=0.25*(h/(cs+vmax)) 

  dt_visc=0.125*(h**2/nu)

  dt=min(dt_force,dt_cfl,dt_visc)

  end if

  end subroutine
   end module
