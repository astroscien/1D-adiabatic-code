
module sweep

contains


subroutine sweep_omega
use defs
use solve
!use exact
implicit none
real (dp) :: omegamin,omegamax
integer :: nomega,i
real (dp), dimension (:),allocatable :: xi_r,omega

omegamin=2e-4
omegamax=1e-2
nomega=10000

if (allocated(xi_r)) deallocate(xi_r)
if (allocated(omega)) deallocate(omega)
allocate( xi_r(nomega),omega(nomega) )

! print out size of response versus forcing frequency
open (13,file='output.data')
do i=1,nomega
  omega(i)=omegamin + (i-1.d0)*(omegamax-omegamin)/(nomega-1.d0)
  call solve_axb(omega(i)**2)
  write(13,*) omega(i)/pi,omega(i)**2*(Rsun**3)/(cgrav*Msun),xi(g%n)!sqrt(sum(xi**2))
  xi_r(i)=xi(g%n)
end do

 call find_max(xi_r,omega,nomega)

 close(13)

end subroutine sweep_omega


subroutine find_max(xi_r,omega,nomega)
use defs
use solve
implicit none
real (dp), dimension (:), intent (in) :: xi_r,omega
integer :: i
integer, intent (in) :: nomega
real (dp) :: diff1,diff2,new_omega,left,right,mid,fl,fm,fr

do i=2,nomega-1
   fl = abs(xi_r(i-1))
   fm = abs(xi_r(i))
   fr = abs(xi_r(i+1))
   diff1 = fm-fl
   diff2 = fr-fm

   if ((diff1 .gt. 0).and.(diff2 .lt. 0)) then
      left = abs(omega(i-1))
      right = abs(omega(i+1))
      mid = abs(omega(i))

      if (fl.lt.fr) then
         left = mid
         fl = fm
         do while ( abs(left-right)>(1.d-8)*(left+right) )
            mid = 0.5d0*(left+right)
	    if (mid==left) then
	       exit
	    elseif (mid==right) then
	       exit
	    else
	    endif

            call solve_axb(mid**2)
            fm = xi(g%n)
            if (fl.lt.fr) then
               left = mid
               fl = fm
            else
               right = mid
               fr = fm
            endif
         enddo
      else
         right = mid
         fr = fm
         do while ( abs(left-right)>(1.d-8)*(left+right) )
            mid = 0.5d0*(left+right)
	    if (mid==left) then
	       exit
	    elseif (mid==right) then
	       exit
	    else
	    endif

            call solve_axb(mid**2)
            fm = xi(g%n)
            if (fl.lt.fr) then
               left = mid
               fl = fm
            else
               right = mid
               fr = fm
            endif
         enddo
      endif        
      print*,mid**2*(Rsun**3)/(cgrav*Msun),fm

   else
      continue
   endif
enddo

end subroutine


end module sweep
