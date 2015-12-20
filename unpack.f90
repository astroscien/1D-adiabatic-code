
module unpack

contains

subroutine unpack_soln(x,xi,psi,phi,gprime)
use defs, only : dp
use grid
implicit none
real (dp), intent (in) :: x(:)
real (dp), intent (out) :: xi(:),psi(:),phi(:),gprime(:)
integer :: i

do i=1,g%n
 xi(i)=x(4*i-3)
 psi(i)=x(4*i-2)
 phi(i)=x(4*i-1)
 gprime(i)=x(4*i)
end do

end subroutine unpack_soln


end module unpack
