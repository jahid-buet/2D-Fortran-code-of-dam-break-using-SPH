module CompGravity
contains
 subroutine gravity(time,p,ntotal)
 use particles
 implicit none
 type(particle),intent(inout)::p(:)
 integer,intent(in)::ntotal
 integer::i
 double precision,parameter::g=9.81
 double precision::fac
 double precision::tdamp=0.050d0 !assuming damping time duration 0.05 s
 double precision,intent(in)::time
 double precision::pi=3.1416

   !if(time<tdamp)then
     !fac=0.5*(sin((-0.5 +time/tdamp)*pi)+1.0)
    !elseif(time>tdamp)then
    fac=1.0
    !end if

    do i=1,ntotal
     if(p(i)%id.lt.2)then !fluids particles
      p(i)%acc%x(1)=p(i)%acc%x(1)+0.0d0           !note:here damping only for external forces such as gravity 
      p(i)%acc%x(2)=p(i)%acc%x(2)-g*fac
     end if
      end do

      end subroutine
       end module


