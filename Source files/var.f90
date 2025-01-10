module properties
contains
subroutine PROP(p,ntotal,rho,rho_boundary,dp,dp1)
use particles
implicit none
type(particle),intent(out)::p(:)
integer,intent(in)::ntotal
double precision,intent(in)::dp,dp1
double precision,intent(in)::rho,rho_boundary
integer::i
do i=1,ntotal
  if(p(i)%id.eq.1)then
  p(i)%dens=rho
  p(i)%mass=p(i)%dens*dp*dp
 else
  p(i)%dens=rho_boundary 
  p(i)%mass=p(i)%dens*dp1*dp1
   end if
 
end do
 end subroutine
  end module
