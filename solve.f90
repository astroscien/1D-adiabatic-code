
module solve
use defs, only : dp
use grid
implicit none
real (dp), allocatable :: xi(:),psi(:),phi(:),gprime(:)

contains

subroutine solve_axb(omegasq)
use defs, only : dp
use grid
use matrix
use unpack
implicit none
real (dp), intent (in) :: omegasq
integer :: nvar,bandwidth
real (dp), allocatable :: al(:,:)
integer, allocatable :: indx(:)
real (dp) :: detsign
integer :: i

if (allocated(xi)) deallocate(xi)
if (allocated(psi)) deallocate(psi)
if (allocated(phi)) deallocate(phi)
if (allocated(gprime)) deallocate(gprime)
allocate( xi(g%n), psi(g%n), phi(g%n), gprime(g%n))

nvar=4*g%n
bandwidth=nleft+nright+1
allocate( al(nvar,nleft), indx(nvar) )
call make_matrix(omegasq)


call bandec(a,nvar,nleft,nright,nvar,bandwidth,al,nleft,indx,detsign)
x=b


call banbks(a,nvar,nleft,nright,nvar,bandwidth,al,nleft,indx,x)
deallocate( al , indx )


call unpack_soln(x,xi,psi,phi,gprime)

end subroutine solve_axb


end module solve
