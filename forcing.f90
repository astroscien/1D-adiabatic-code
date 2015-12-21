
module forcing
use defs, only : dp
implicit none
type forcingtype
	real (dp), allocatable :: Ur(:)
end type forcingtype
type (forcingtype) :: f

contains


subroutine setup_forcing
use defs
use grid
implicit none
integer :: i,j

if (allocated(f%Ur)) deallocate(f%Ur)
allocate(f%Ur(g%n))

f%Ur(:)=0.d0
f%ell=2

do i=1,g%n
   f%Ur(i)=-cgrav*Msun/Rsun*(g%r(i)/Rsun)**f%ell
end do

end subroutine setup_forcing

end module forcing
