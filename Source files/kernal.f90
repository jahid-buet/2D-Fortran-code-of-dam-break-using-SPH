module kernal
contains
subroutine ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)
implicit none
double precision,intent(in)::h
double precision,intent(in)::rij,rij_x,rij_y !distance between two particles
double precision,intent(out)::wij !kernal
double precision,intent(out)::dwdx,dwdy!kernal derivative
double precision::q,fac,dwdq,tmp,h1,facd,qq1,qq2,qqs,tmp1
integer,intent(in)::kernal_type
if(kernal_type.eq.1)then
 !quintic spline kernal
    h1=1./h
    q= h1*dsqrt(rij)
    fac=7.0/(478.0*3.1416*h*h)
    !quintic kernal
     if(q.gt.0.0 .and. q.le.1.0)then
   wij=fac*((3.-q)**5 -6.*(2.-q)**5 +15.*(1.-q)**5)

     !kernal gradient
    dwdq=fac*((-5.0*(3.0-q)**4) +(30.0*(2.0-q)**4) -(75.0*(1.0-q)**4))
    tmp=dwdq*h1 /dsqrt(rij)
    dwdx=tmp*rij_x
    dwdy=tmp*rij_y
    !dwdx=dwdq/h*(rij_x/dsqrt(rij))
    !dwdy=dwdq/h*(rij_y/dsqrt(rij))

  
     elseif(q.gt.1.0 .and.q.le.2.0)then
        wij=fac*((3.-q)**5 -(6.*(2.-q)**5))

       !kernal gradient
        dwdq=fac*((-5.0*(3.0-q)**4) + (30.0*(2.0-q)**4))
        tmp=dwdq*h1 /dsqrt(rij)
        dwdx=tmp*rij_x
        dwdy=tmp*rij_y
       !dwdx=dwdq/h*(rij_x/dsqrt(rij))
       !dwdy=dwdq/h*(rij_y/dsqrt(rij))

     elseif(q.gt.2.0 .and. q.le.3.0)then
       wij=fac*((3.0-q)**5)

       !kernal gradient
        dwdq=fac*(-5.0*(3.0-q)**4)
        tmp=dwdq*h1 /dsqrt(rij)
        dwdx=tmp*rij_x
        dwdy=tmp*rij_y
        
       !dwdx=dwdq/h*(rij_x/dsqrt(rij))
       !dwdy=dwdq/h*(rij_y/dsqrt(rij))
       
     else
      wij=0.0
      dwdx=0.0
      dwdy=0.0
   end if

 elseif(kernal_type.eq.2)then
!wendland spline kernal
    h1=1./h
    q= h1*dsqrt(rij)
    fac=0.557*h1*h1
    tmp1=1.0-0.5*q
     if(q.ge.0.0 .and.q.le.2.0)then     
      wij=fac*tmp1*tmp1*tmp1*tmp1*(2.0*q+1.0)
      dwdq=-5.0*q*tmp1*tmp1*tmp1
      tmp=dwdq*h1 /dsqrt(rij)      
      dwdx=fac*rij_x*tmp
      dwdy=fac*rij_y*tmp
 else
 wij=0.0
 dwdx=0.0
 dwdy=0.0
 end if

  end if
 !hyperbolic-shaped kernal
!$$$$$$     h1=1./h
!$$$$$$     q= h1*dsqrt(rij)
!$$$$$$     fac=0.1061*h1*h1
!$$$$$$     if(q.ge.0.0 .and.q.le.1.0)then 
!$$$$$$       wij=fac*(q*q*q -6.0*q +6.0d0)
!$$$$$$       dwdq=3.0d0*q*q -6.0d0
!$$$$$$       tmp=dwdq*h1 /dsqrt(rij)
!$$$$$$       dwdx=fac*rij_x*tmp
!$$$$$$       dwdy=fac*rij_y*tmp
!$$$$$$      elseif(q.gt.1.0 .and.q.le.2.0)then
!$$$$$$       wij=fac*((2.0-q)*(2.0-q)*(2.0-q))
!$$$$$$       dwdq=-3.0*(2.0-q)*(2.0-q)
!$$$$$$       tmp=dwdq*h1 /dsqrt(rij)
!$$$$$$       dwdx=fac*rij_x*tmp
!$$$$$$       dwdy=fac*rij_y*tmp
!$$$$$$       else
!$$$$$$        wij=0.0
!$$$$$$        dwdx=0.0
!$$$$$$        dwdy=0.0
!$$$$$$        end if 

    end subroutine


   
      end module
 