module particles
use vectors

type:: particle !definition of  particles
type(vector_3d)::coord
type(vector_3d)::vel
type(vector_3d)::acc
double precision::press
double precision::mass
double precision::dens
double precision::mat(3,3)!matrix for kernal gradient correction
double precision::txx
double precision::txy
double precision::tyy!stress rate for particles
double precision::hxx
double precision::hxy  !strain rate for particles
double precision::hyy
integer::id
 end type

!definition of boundary
 !type::boundary
 !type(vector_3d)::coord,vel,acc
 !double precision::h,press,mass,dens,vol
 !!end type
 !assignment of double operator only in boundary particles

  !interface assignment(=)
  !module procedure scalar_to_boundary
  !end interface
  !contains
   !subroutine scalar_to_boundary(a,b)
   !type(boundary),intent(out)::a
   !double precision,intent(in)::b
   !a%coord%x(1:3)=b
   !end subroutine
   
  end module





