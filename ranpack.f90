module ranpack
!---------------------------------------------
! random number generator
!---------------------------------------------
! HISTORY
!  19??/ ?/ ?  [RandomNumber2] written by someone
!  2007/ 4/30  [RandomNumber3] edited by Koji Kobayashi
!  2008/ 9/24  [ranpack v1.0]
!  2008/ 9/25  [v1.1]
!  2008/10/ 3  [v1.1.1] (add 'save'. not important)
!---------------------------------------------
! USAGE
!  call ran__dini(iseed)
!  call ran__dnext(100)   ! initialize DOUBLE random number generation routine
!  x = ran__drnd()        ! DOUBLE random number (0,1) generation
!  call ran__iini(iseed)  ! initialize INTEGER random number generation routine
!  k = ran__irnd()        ! INTEGER random number in [1,2^{31}-1] generation
!---------------------------------------------
implicit none
private
 integer,parameter:: double=selected_real_kind(15)
 integer,parameter:: ip=250,iq=103,irbit=31,ilcbit=20
 integer,save :: ir(ip), irq(ip), inext(ip), ipos
 integer,save :: iseed
public:: ran__drnd, ran__irnd, ran__dini, ran__iini, ran__dnext
contains

! *********************************************************

!---------------------------------
! initialize the random number generation
! tausworthe sequence
! IR_{n}=IR_{n-IP} .xor. IR_{n-IQ}
! the bit-width is given by irbit
!---------------------------------
 subroutine ran__dini(irseed)
 integer,intent(in) :: irseed
 integer            :: i,j,imask,ilc
   ! initialize I-IQ position
    do i = 1, ip
       irq(i) = i - iq
       if (irq(i)<1) irq(i) = irq(i) + ip
    end do
   ! initialize next position
    do i = 1, ip - 1
       inext(i) = i + 1
    end do
    inext(ip) = 1
   ! initialize random seed
    call ran__iini(irseed)
    ir(1:ip) = 0
    imask = 2**(ilcbit-1)
    do i = 1, irbit
       do j = 1, ip
          ilc = ran__irnd()   !INTEGER random number
          ilc = ilc/imask
          ilc = iand(ilc,1)
          ir(j) = ir(j)*2 + ilc
       end do
    end do
    ipos = 1
 end subroutine ran__dini

!--------------------------------
! initialize the INTEGER random number generation routine
!--------------------------------
 subroutine ran__iini(irseed)
 integer,intent(in) :: irseed
    iseed = irseed
 end subroutine ran__iini

!--------------------------------
! INTEGER random number in [1,2^{31}-1] generation
!--------------------------------
 function ran__irnd()
 integer :: ran__irnd
    iseed = iseed*48828125
    if (iseed<=0) iseed = (iseed+2147483647) + 1
    ran__irnd = iseed
 end function ran__irnd


!---------------------------------
! DOUBLE random number (0 ,1) generation
!---------------------------------
 function ran__drnd()
 real(kind=double) :: ran__drnd
 real(kind=double),parameter:: dconst=1.0_double/(2.0_double**irbit)
    ir(ipos) = ieor(ir(ipos),ir(irq(ipos)))
    ran__drnd = dble(ir(ipos))*dconst
    ipos = inext(ipos)
 end function ran__drnd

!---------------------------------
! idling 
!---------------------------------
 subroutine ran__dnext(n)
 integer :: i,n
   do i=1,n
    ir(ipos) = ieor(ir(ipos),ir(irq(ipos)))
    ipos = inext(ipos)
   end do
 end subroutine ran__dnext

end module ranpack
