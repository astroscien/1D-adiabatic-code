
!	this program finds the oscillation modes for 1-d cartesian sound waves.

! compile: gfortran defs.f90 bandec.f banbks.f grid.f90 forcing.f90 matrix.f90 unpack.f90 solve.f90 exact.f90 sweep.f90 main.f90

program sho
use grid
use forcing
use matrix
use sweep
implicit none
integer :: i

 call grid_setup
 call setup_forcing
 call setup_matrix
 call sweep_omega

end program
