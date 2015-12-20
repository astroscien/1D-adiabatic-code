# 1D-adiabatic-code
!
! mathematical and physical constants (in cgs)
!
! the 2006 CODATA recommended values of the physical constants
! by cohen & taylor

! based on mesa/const/public/const_defs.f

module defs
implicit none

integer, parameter :: sp=selected_real_kind(p=5)	! 5 digits
integer, parameter :: dp=selected_real_kind(p=15)	! 15 digits
integer, parameter :: i4=selected_int_kind(9)		! short integer, -10^9 (exclusive) to 10^9 (exclusive)
integer, parameter :: i8=selected_int_kind(14)		! long integer, -10^{14} (exclusive) to 10^{14} (exclusive)

real (dp), parameter :: pi=2.d0*asin(1.d0)
real (dp), parameter :: cgrav=6.67428d-8		! gravitational constant, g^-1 cm^3 s^-2
real (dp), parameter :: planck_h = 6.62606896D-27	! planck's constant, erg s
real (dp), parameter :: hbar = planck_h/(2*pi)		! hbar, erg s
real (dp), parameter :: qe = 4.80320440D-10		! proton charge, esu = (g cm^3 s^-2)^(1/2)
real (dp), parameter :: clight = 2.99792458d10		! speed of light, cm s^{-1}
real (dp), parameter :: boltzm = 1.3806504D-16		! k_B, Boltzmann's constant, erg K^{-1}
real (dp), parameter :: avo = 6.02214179d23		! avogadro's number
real (dp), parameter :: cgas = boltzm*avo		! the gas constant k_B/m_proton
real (dp), parameter :: kev = 8.617385d-5		! eV per Kelvin
real (dp), parameter :: amu = 1.660538782d-24		! atomic mass unit, g
real (dp), parameter :: mp = 1.6726231d-24 		! proton mass, g
real (dp), parameter :: me = 9.1093826D-28 		! electron mass, g
real (dp), parameter :: boltz_sigma = 5.670400D-5	! boltzmann's sigma = a*c/4, erg cm^-2 K^-4 s^-1
real (dp), parameter :: crad = boltz_sigma*4/clight	! a=radiation constant, erg cm^-3 K^-4

real (dp), parameter :: msun = 1.9892d33  		! solar mass (g)
real (dp), parameter :: rsun = 6.9598d10 		! solar radius (cm)
real (dp), parameter :: lsun = 3.8418d33		! solar luminosity (erg s^-1)
real (dp), parameter :: agesun = 4.57d9			! solar age (years)
real (dp), parameter :: mearth = 5.9764d27 		! earth mass (g)
real (dp), parameter :: rearth = 6.37d8 		! earth radius (cm)
real (dp), parameter :: au = 1.495978921d13 		! astronomical unit (cm)
real (dp), parameter :: mjup = 1.8986d30 		! jupiter mass (g)
real (dp), parameter :: rjup = 6.9911d9 		! jupiter mean radius (cm)

real (dp), parameter :: secyer = 3.1558149984d7 	! seconds per year
real (dp), parameter :: secday=86400.d0			! seconds per solar day

end module defs
