
module grid
use defs, only : dp
implicit none
type gridtype
  integer :: n
!  real (dp) :: delta
  real (dp), allocatable :: delta(:)
  real (dp), allocatable :: r(:),g0(:),Nsq(:),csq(:),rho(:),khsq(:),Ur(:)
end type gridtype
type (gridtype) :: g

contains

subroutine grid_setup
implicit none
integer :: n,i,iostatus
real (dp) :: mytemp


!namelist /grid/ n
!open(11,file='input')
!read (11,nml=grid)
!close(11)
!g%n=n

!g%delta = 1.d0/(g%n-1.d0)

open(11,file='poly.data')

i=0 
do  ! how many points in file
  i=i+1
  read(11,*,iostat=iostatus) mytemp 
  if (iostatus > 0 ) then
     print *,"PROBLEM WITH READ"
     stop
  else if (iostatus < 0) then
     exit
  endif 
enddo

g%n=i-1
if (allocated(g%r)) deallocate(g%r)
if (allocated(g%g0)) deallocate(g%g0)
if (allocated(g%Nsq)) deallocate(g%Nsq)
if (allocated(g%csq)) deallocate(g%csq)
if (allocated(g%rho)) deallocate(g%rho)
if (allocated(g%khsq)) deallocate(g%khsq)
if (allocated(g%Ur)) deallocate(g%Ur)

allocate( g%r(g%n),g%g0(g%n),g%Nsq(g%n),g%csq(g%n),g%rho(g%n),g%khsq(g%n),g%Ur(g%n)) ! allocate variables
allocate( g%delta(g%n-1) )

rewind(11)

do i=1,g%n
   read(11,*)g%r(i),g%g0(i),g%Nsq(i),g%csq(i),g%rho(i),g%khsq(i),g%Ur(i)
!   print*,
enddo

close(11)

do i=1,g%n-1
   g%delta(i) = g%r(i+1)-g%r(i)
end do

!do i=1,g%n
! g%r(i) = g%delta * (i-1.d0)
!end do

end subroutine grid_setup

end module grid
