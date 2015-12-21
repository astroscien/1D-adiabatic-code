
module matrix
use defs, only : dp
implicit none
integer, parameter :: ixi = 1, ipsi = 2, iphi = 3, igprime = 4
integer :: nleft,nright
real (dp), allocatable :: a(:,:),x(:),b(:)

contains


subroutine make_matrix(omegasq)
use defs
use grid
use forcing
implicit none
real (dp) :: omegasq
integer ::  i,eqnnum,varnum,ndiag

a=0.d0
b=0.d0
x=0.d0

ndiag=nleft+1

! at r=1e-3 boundary condition
eqnnum=1
varnum=iphi
a(eqnnum,varnum-eqnnum+ndiag) = f%ell/g%r(1)	
varnum=igprime
a(eqnnum,varnum-eqnnum+ndiag) = -1.d0

b(eqnnum)=0.d0


eqnnum=2
varnum=ixi
a(eqnnum,varnum-eqnnum+ndiag) = omegasq
varnum=ipsi
a(eqnnum,varnum-eqnnum+ndiag) = -f%ell/g%r(1)

b(eqnnum)=0.d0

! equations differenced at n-1 points
do i=1,g%n-1


  eqnnum=4*i-1
  varnum=ixi + 4*(i-1)					! 4i-3
  a(eqnnum,varnum-eqnnum+ndiag) = 0.25d0*(g%Nsq(i+1)+g%Nsq(i))-0.5d0*omegasq
  varnum=ixi + 4*(i)					! 4i+1
  a(eqnnum,varnum-eqnnum+ndiag) = 0.25d0*(g%Nsq(i+1)+g%Nsq(i))-0.5d0*omegasq
  varnum=ipsi + 4*(i-1)					! 4i-2
  a(eqnnum,varnum-eqnnum+ndiag) = - 1.d0/g%delta(i) - 0.5d0*(g%Nsq(i+1)+g%Nsq(i))/(g%g0(i+1)+g%g0(i))
  varnum=ipsi + 4*(i)					! 4i+2
  a(eqnnum,varnum-eqnnum+ndiag) = 1.d0/g%delta(i) - 0.5d0*(g%Nsq(i+1)+g%Nsq(i))/(g%g0(i+1)+g%g0(i))
  varnum=iphi + 4*(i-1)					! 4i-1
  a(eqnnum,varnum-eqnnum+ndiag) = 0.5d0*(g%Nsq(i+1)+g%Nsq(i))/(g%g0(i+1)+g%g0(i))
  varnum=iphi + 4*(i)					! 4i+3
  a(eqnnum,varnum-eqnnum+ndiag) = 0.5d0*(g%Nsq(i+1)+g%Nsq(i))/(g%g0(i+1)+g%g0(i))

  b(eqnnum) = -(g%Nsq(i+1)+g%Nsq(i))/(g%g0(i+1)+g%g0(i))*0.5d0*(f%Ur(i+1)+f%Ur(i))


  eqnnum=4*i
  varnum=ixi + 4*(i-1)					! 4i-3
  a(eqnnum,varnum-eqnnum+ndiag) = - 1.d0/g%delta(i) + 0.5d0*( 4.d0/(g%r(i+1)+g%r(i))-(g%g0(i+1)+g%g0(i))/(g%csq(i+1)+g%csq(i)) )
  varnum=ixi + 4*(i)					! 4i+1
  a(eqnnum,varnum-eqnnum+ndiag) = 1.d0/g%delta(i) + 0.5d0*( 4.d0/(g%r(i+1)+g%r(i))-(g%g0(i+1)+g%g0(i))/(g%csq(i+1)+g%csq(i)) )
  varnum=ipsi + 4*(i-1)					! 4i-2
  a(eqnnum,varnum-eqnnum+ndiag) = 0.5d0*(2.d0/(g%csq(i+1)+g%csq(i))-0.5d0*(g%khsq(i+1)+g%khsq(i))/omegasq )
  varnum=ipsi + 4*(i)					! 4i+2
  a(eqnnum,varnum-eqnnum+ndiag) = 0.5d0*(2.d0/(g%csq(i+1)+g%csq(i))-0.5d0*(g%khsq(i+1)+g%khsq(i))/omegasq )
  varnum=iphi + 4*(i-1)					! 4i-1
  a(eqnnum,varnum-eqnnum+ndiag) = -1.d0/(g%csq(i+1)+g%csq(i))
  varnum=iphi + 4*(i)					! 4i+3
  a(eqnnum,varnum-eqnnum+ndiag) = -1.d0/(g%csq(i+1)+g%csq(i))

  b(eqnnum) = (f%Ur(i+1)+f%Ur(i))/(g%csq(i+1)+g%csq(i))


  eqnnum=4*i+1
  varnum=iphi + 4*(i-1)					! 4i-1
  a(eqnnum,varnum-eqnnum+ndiag) = - 1.d0/g%delta(i)
  varnum=iphi + 4*(i)					! 4i+3
  a(eqnnum,varnum-eqnnum+ndiag) = 1.d0/g%delta(i)
  varnum=igprime + 4*(i-1)					! 4i
  a(eqnnum,varnum-eqnnum+ndiag) = -0.5d0
  varnum=igprime + 4*(i)					! 4i+4
  a(eqnnum,varnum-eqnnum+ndiag) = -0.5d0

  b(eqnnum) = 0.d0


  eqnnum=4*i+2
  varnum=ixi + 4*(i-1)					! 4i-3
  a(eqnnum,varnum-eqnnum+ndiag) = pi*cgrav*(g%rho(i+1)+g%rho(i))*(g%Nsq(i+1)+g%Nsq(i))/(g%g0(i+1)+g%g0(i))
  varnum=ixi + 4*(i)					! 4i+1
  a(eqnnum,varnum-eqnnum+ndiag) = pi*cgrav*(g%rho(i+1)+g%rho(i))*(g%Nsq(i+1)+g%Nsq(i))/(g%g0(i+1)+g%g0(i))
  varnum=ipsi + 4*(i-1)					! 4i-2
  a(eqnnum,varnum-eqnnum+ndiag) = 2.d0*pi*cgrav*(g%rho(i+1)+g%rho(i))/(g%csq(i+1)+g%csq(i))
  varnum=ipsi + 4*(i)					! 4i+2
  a(eqnnum,varnum-eqnnum+ndiag) = 2.d0*pi*cgrav*(g%rho(i+1)+g%rho(i))/(g%csq(i+1)+g%csq(i))
  varnum=iphi + 4*(i-1)					! 4i-1
  a(eqnnum,varnum-eqnnum+ndiag) = 0.25d0*(g%khsq(i+1)+g%khsq(i)) - 2.d0*pi*cgrav*(g%rho(i+1)+g%rho(i))/(g%csq(i+1)+g%csq(i))
  varnum=iphi + 4*(i)					! 4i+3
  a(eqnnum,varnum-eqnnum+ndiag) = 0.25d0*(g%khsq(i+1)+g%khsq(i)) - 2.d0*pi*cgrav*(g%rho(i+1)+g%rho(i))/(g%csq(i+1)+g%csq(i))
  varnum=igprime + 4*(i-1)					! 4i
  a(eqnnum,varnum-eqnnum+ndiag) = 1.d0/g%delta(i) - 2.d0/(g%r(i+1)+g%r(i))
  varnum=igprime + 4*(i)					! 4i+4
  a(eqnnum,varnum-eqnnum+ndiag) = -1.d0/g%delta(i) - 2.d0/(g%r(i+1)+g%r(i))

  b(eqnnum) = 2.d0*pi*cgrav*(g%rho(i+1)+g%rho(i))/(g%csq(i+1)+g%csq(i))*(f%Ur(i+1)+f%Ur(i))

end do

! at r=R boundary condition
eqnnum=4*g%n-1
varnum= -1 + 4*g%n					! 4N-1
a(eqnnum,varnum-eqnnum+ndiag) = (f%ell+1.d0)/g%r(g%n)
varnum= 4*g%n						! 4N
a(eqnnum,varnum-eqnnum+ndiag) = 1.d0

b(eqnnum)=0.d0


eqnnum=4*g%n
varnum= -3 + 4*g%n					! 4N-3
a(eqnnum,varnum-eqnnum+ndiag) = g%g0(g%n)
varnum= -2 + 4*g%n					! 4N-2
a(eqnnum,varnum-eqnnum+ndiag) = -1.d0
varnum= -1 + 4*g%n					! 4N-1
a(eqnnum,varnum-eqnnum+ndiag) = 1.d0

b(eqnnum)= - f%Ur(g%n)

!do i = 1,4*g%n
!print*,a(i,:)
!enddo

end subroutine make_matrix



subroutine setup_matrix
use grid
implicit none
integer :: nvar

nvar=4*g%n

! keep track of these as you make the matrix
nleft = 7
nright = 7

if (allocated(a)) deallocate(a)
if (allocated(x)) deallocate(x)
if (allocated(b)) deallocate(b)
allocate( a(nvar,nleft+nright+1),x(nvar),b(nvar) )

end subroutine setup_matrix

end module matrix
