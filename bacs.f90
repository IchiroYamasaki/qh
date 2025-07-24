

!--------------------------------------
! Basic constant set
! 2009/10/16 v2.0
!--------------------------------------
module bacs
implicit none
public
integer, parameter :: single=selected_real_kind(6), double=selected_real_kind(15)
real(kind=double), parameter :: &
   & zero=0.0_double, one=1.0_double, two=2.0_double, three=3._double, four=4._double, &
   & ten=10._double, half=0.5_double, &
   & root2=1.414213562373095048802_double,      root3 =1.732050807568877293527_double,&
   & Pi   =3.14159265358979323846264338_double, &
   & twoPi=6.28318530717958647692528676_double, rootPi=1.7724538509055160272981675_double
complex(kind=double), parameter :: img=(0.0_double,1.0_double), &
     & cOne =(1.0_double,0.0_double), cZero=(0.0_double,0.0_double)

end module bacs
