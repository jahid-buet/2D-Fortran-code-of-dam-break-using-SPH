module ini

contains
subroutine initialize(p,ntotal)
use particles
type(particle),intent(inout)::p(:)
integer,intent(in)::ntotal
integer::i
 
!initialize velocity,pressure and acceleration of each particlesin the domain
 do i=1,ntotal
   p(i)%vel%x(1)=0.0d0
   p(i)%vel%x(2)=0.0d0
   p(i)%acc%x(1)=0.0d0
   p(i)%acc%x(2)=0.0d0
   p(i)%press=0.0d0
  end do

 end subroutine
  end module
 

