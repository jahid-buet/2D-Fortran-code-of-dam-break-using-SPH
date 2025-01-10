module vectors
private
public::operator(-),operator(+),operator(*),operator(.cross.),operator(.dot.),assignment(=),operator(.root.),operator(/)
  type,public:: vector_3d
  double precision::x(3)
  end type

  type,public:: vector_2d
  double precision::x(2)
  end type



interface operator(-)
module procedure vector_3d_minus
module procedure vector_2d_minus
end interface

interface operator(+)
module procedure vector_3d_plus
module procedure vector_2d_plus
module procedure vector_plus_vector_multiply_scalar
end interface

interface operator(*)
module procedure vector_3d_multiply_scalar_3d
module procedure vector_2d_multiply_scalar_2d
!$$$$$$ module procedure scalar_2d_multiply_scalar_2d
end interface

interface operator(.cross.)
module procedure cross_product
end interface

interface operator(.dot.)
module procedure dot_product_3d
module procedure dot_product_2d
end interface

interface assignment(=)
module procedure vector_3d_to_scalar_3d
 module procedure scalar_to_vector
module procedure vector_3d_to_vector_3d
module procedure vector_2d_to_scalar_2d
module procedure vector_2d_to_vector_2d
end interface

interface operator(.root.)
module procedure scalar_root
end interface

interface operator(/)
module procedure division_by_scalar_3d
module procedure division_by_scalar_2d
end interface

contains

subroutine vector_3d_to_scalar_3d(a,b)
double precision,intent(out)::a(3)
type(vector_3d),intent(in)::b
a(1)=b%x(1)
a(2)=b%x(2)
a(3)=b%x(3)
end subroutine

 subroutine vector_2d_to_scalar_2d(a,b)
 double precision,intent(out)::a(2)
 type(vector_2d),intent(in)::b
 a(1)=b%x(1)
 a(2)=b%x(2)
 end subroutine
 

 subroutine scalar_to_vector(a,b)
 type(vector_3d),intent(out)::a
 double precision,intent(in)::b
 a%x(1)=b
 a%x(2)=b
 a%x(3)=b
end subroutine

subroutine vector_3d_to_vector_3d(a,b) !vector to vector assignment
type(vector_3d),intent(out)::a
type(vector_3d),intent(in)::b
a%x(1)=b%x(1)
a%x(2)=b%x(2)
a%x(3)=b%x(3)
end subroutine

  subroutine  vector_2d_to_vector_2d(a,b) !vector to vector assignment
  type(vector_2d),intent(out)::a
  type(vector_2d),intent(in)::b
   a%x(1)=b%x(1)
   a%x(2)=b%x(2)
  end subroutine

function vector_3d_minus(a,b)
type(vector_3d),intent(in)::a
type(vector_3d),intent(in)::b
type(vector_3d)::vector_3d_minus
vector_3d_minus%x(1)=a%x(1)-b%x(1)
vector_3d_minus%x(2)=a%x(2)-b%x(2)
vector_3d_minus%x(3)=a%x(3)-b%x(3)
end function
  
   function vector_2d_minus(a,b)
   type(vector_2d),intent(in)::a
   type(vector_2d),intent(in)::b
   type(vector_2d)::vector_2d_minus
   vector_2d_minus%x(1)=a%x(1)-b%x(1)
   vector_2d_minus%x(2)=a%x(2)-b%x(2)
   end function
   
function vector_3d_plus(a,b)
type(vector_3d),intent(in)::a
type(vector_3d),intent(in)::b
type(vector_3d)::vector_3d_plus
vector_3d_plus%x(1)=a%x(1)+b%x(1)
vector_3d_plus%x(2)=a%x(2)+b%x(2)
vector_3d_plus%x(3)=a%x(3)+b%x(3)
end function

     function vector_2d_plus(a,b)
     type(vector_2d),intent(in)::a
     type(vector_2d),intent(in)::b
     type(vector_2d)::vector_2d_plus
     vector_2d_plus%x(1)=a%x(1)+b%x(1)
     vector_2d_plus%x(2)=a%x(2)+b%x(2)
     end function

  function vector_plus_vector_multiply_scalar(a,b)
  double precision,intent(in)::b(4)
  type(vector_2d),intent(in)::a
  type(vector_2d)::vector_plus_vector_multiply_scalar
  
  integer::i
  vector_plus_vector_multiply_scalar%x(1)=0
  vector_plus_vector_multiply_scalar%x(2)=0
  do i=1,4 !no.of sum term
 vector_plus_vector_multiply_scalar%x(1)=vector_plus_vector_multiply_scalar%x(1)+b(i)*a%x(1)
 vector_plus_vector_multiply_scalar%x(2)=vector_plus_vector_multiply_scalar%x(2)+b(i)*a%x(2)
 end do 
!$$$$$$ vector_plus_vector_multiply_scalar%x(2)=vector_plus_vector_multiply_scalar%x(2)+b*a%x(2)
!$$$$$$ vector_plus_vector_multiply_scalar%x(3)=vector_plus_vector_multiply_scalar%x(3)+b*a%x(3)
 end function
 
 function vector_3d_multiply_scalar_3d(a,b)
 type(vector_3d),intent(in)::a
 double precision,intent(in)::b
 type(vector_3d)::vector_3d_multiply_scalar_3d
 vector_3d_multiply_scalar_3d%x(1)=b*a%x(1)
 vector_3d_multiply_scalar_3d%x(2)=b*a%x(2)
 vector_3d_multiply_scalar_3d%x(3)=b*a%x(3)

 end function

    function vector_2d_multiply_scalar_2d(a,b)
    type(vector_2d),intent(in)::a
    double precision,intent(in)::b
    type(vector_2d)::vector_2d_multiply_scalar_2d
    vector_2d_multiply_scalar_2d%x(1:2)=b*a%x(1:2)
   
    end function
    
 

!$$$$$$    function scalar_2d_multiply_scalar_2d(a,b)
!$$$$$$   double precision,intent(in)::a
!$$$$$$   double precision,intent(in)::b
!$$$$$$   double precision::scalar_2d_multiply_scalar_2d
!$$$$$$    scalar_2d_multiply_scalar_2d=a*b
!$$$$$$   end function

function cross_product(a,b)
type(vector_3d),intent(in)::a
type(vector_3d),intent(in)::b
type(vector_3d)::cross_product
cross_product%x(1)=a%x(2)*b%x(3)-b%x(2)*a%x(3)
cross_product%x(2)=b%x(1)*a%x(3)-a%x(1)*b%x(3)
cross_product%x(3)=a%x(1)*b%x(2)-b%x(1)*a%x(2)
end function

function dot_product_3d(a,b)
type(vector_3d),intent(in)::a
type(vector_3d),intent(in)::b
double precision::dot_product_3d
dot_product_3d=a%x(1)*b%x(1)+a%x(2)*b%x(2)+a%x(3)*b%x(3)
end function

function dot_product_2d(a,b)
type(vector_2d),intent(in)::a
type(vector_2d),intent(in)::b
double precision::dot_product_2d
dot_product_2d=a%x(1)*b%x(1)+a%x(2)*b%x(2)
end function


 function scalar_root(a)
 double precision,intent(in)::a
 double precision::scalar_root
 scalar_root=sqrt(a)
 end function
 
 function division_by_scalar_3d(a,b)
 type(vector_3d),intent(in)::a
 double precision,intent(in)::b
 type(vector_3d)::division_by_scalar_3d
 
 division_by_scalar_3d%x(1)=(a%x(1))/b
 division_by_scalar_3d%x(2)=(a%x(2))/b
 division_by_scalar_3d%x(3)=(a%x(3))/b
 
 end function
  
   function division_by_scalar_2d(a,b)
   type(vector_2d),intent(in)::a
   double precision,intent(in)::b
   type(vector_2d)::division_by_scalar_2d
   division_by_scalar_2d%x(1)=(a%x(1))/b
   division_by_scalar_2d%x(2)=(a%x(2))/b
   end function
 
end module

