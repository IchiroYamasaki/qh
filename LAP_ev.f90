module LAP_ev
!---------------------------------------------
! Eigenvalue problem solver
!  for hermite or symmetric matrics
! (f90 wrapper of LAPACK 'ZHEEV', 'ZHEEVR/X', 'DSYEV')
!---------------------------------------------
! Author
!  Koji Kobayashi
!---------------------------------------------
! History
!  2008/ 9/25 [LAP_zheev v1.0] <new>
!  2009/ 7/ 1 [v1.1] <clean> ini 
!  2011/ 1/ 6 [LAP_ev v1.2] <add> variations
!  2014/ 5/11 [v1.3] <add> zheevx
!  2014/ 5/12 [v1.4] <add> zheevr
!  2015/ 6/13 [v1.5] <add> LAP_dsyev_1 as 'LAP__ev'
!  2015/11/24 [v1.6] <fix> in LAP_zheev_1(_v) and LAP_dsyev_1(_v), workspace: stack -> heap
!                    <interface> LAP_zheevr as 'LAP__ev'
!  2016/ 5/16 [v1.7] <change> save original 'A' as default.
!  2017/ 1/18 [v1.8] <feature> OpenMP multithreading available
!  2017/ 8/21 [v1.8.1] <fix> routine name: LAP__ev_free
!---------------------------------------------
! Requirement
!  LAPACK
!---------------------------------------------
! about ZHEEV
!  diagonalizes dcomplex Hermite matrix 'A'
!   and returns its eigenvalues [and eigenvectors in 'A'].
!  Order of eigenvalues : from the smallest to largest
!---------------------------------------------
! Usage
!
!--- simplest case (for single use) ---
!  call LAP__ev(N,A,egVals[,egVecs]) !--- eigenvalues [& eigenvectors] of A
!
!  !-- or, to minimize memory, --
!  call LAP__ev(N,A,egVals,JOBZ='V/N') !--- 'V': A is eigenvectors, 'N': A is destroyed
!
!--- for iterative use ---
!  call LAP__zheev_ini(N)  !--- initialize
!  do
!     {A = ...}   !--- define A
!     call LAP__zheev(A,egVals[,JOBZ='V/N'])  !--- eigenvalues of A [A is destroyed]
!  end do
!
!--- find Imin-th to Imax-th eigenvalues ---
! call LAP__zheevr(A,Imin,Imax,egVals[,egVecs])
!
!
!--- for parallel use ---
!  call LAP__ev(N,A,egVals[,egVecs])
!
!  !-- or, (not sure...) --
!  !$OMP parallel
!  call LAP__zheev_ini(N)  !--- initialize
!  !$OMP do private(A,egVals)
!  do
!     {A = ...}   !--- define A
!     call LAP__zheev(A,egVals[,JOBZ='V/N'])  !--- eigenvalues of A [A is destroyed]
!  end do
!  !$OMP end do
!  !$OMP end parallel
!
!---------------------------------------------
! ZHEEVX, ZHEEVR
!  restrict the number of eigenvalues to be calculated (faster)
!  I don't know which (X/R) is better
!---------------------------------------------
! Usage
!  "call LAP__ev(A,Imin,Imax,egVals[,egVecs])"
! for Imin-th to Imax-th egVals  (stored in egVals(1:Imax-Imin+1) of egVals(1:N))
! or,
!  "call LAP__ev(A,Vmin,Vmax,nEig,egVals[,egVecs])"
! for all(=nEig) egVals in the range (Vmin,Vmax]
!
! (to use ZHEEVX, "call LAP__zheevx")
!----------------
!
!
!---------------------------------------------
! ZHEEV ARGUMENTS
!
!  A       [inout]  dcomplex (n,n) matrix
!          On entry, the matrix A.
!          If UPLO = 'U', the upper triangular part of A
!             contains the upper triangular part of the matrix A.
!          If UPLO = 'L', the lower triangular part of A 
!             contains the lower triangular part of the matrix A.
!          On exit:
!          If JOBZ = 'V', then the columns of A contain the orthonormal 
!             eigenvectors of the matrix A in the order of the eigenvalues.
!          If JOBZ = 'N', then the 'UPLO' triangle of A, 
!             including the diagonal, is destroyed.
!
!  egVals  [output]  dreal (n) vector
!          The eigenvalues in ascending order.
!
!  JOBZ    [(optional)input]  CHARACTER*1
!          = 'N':  Compute eigenvalues only; (default of this module)
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  UPLO    [(optional)input]  CHARACTER*1
!          = 'U':  Upper triangle of A is stored; (default of this module)
!          = 'L':  Lower triangle of A is stored.
!
!---------------------------------------------
 !$ use omp_lib 
 implicit none
 private
  integer, parameter :: double=selected_real_kind(15)
  !-- for iteration --
  integer :: LWORK, N
  complex(kind=double), allocatable :: zWork(:)  !--- for zheev, (LWORK)
  real   (kind=double), allocatable :: dWork(:)  !--- for zheev, (3*N-2)
  !$OMP threadprivate(LWORK,N,zWork,dWork)
 external :: zheev, zheevx, zheevr, dsyev

 INTERFACE LAP__ev
    module procedure LAP_zheev_1,LAP_zheev_1_v, LAP_dsyev_1,LAP_dsyev_1_v, LAP_zheevr_I,LAP_zheevr_V
 end INTERFACE
 INTERFACE LAP__zheev_ini
    module procedure ini_A, ini_N 
 end INTERFACE
 INTERFACE LAP__zheev
    module procedure LAP_zheev_run,LAP_zheev_1,LAP_zheev_1_v
 end INTERFACE
 INTERFACE LAP__zheevx
    module procedure LAP_zheevx_I,LAP_zheevx_V
 end INTERFACE
 INTERFACE LAP__zheevr
    module procedure LAP_zheevr_I,LAP_zheevr_V
 end INTERFACE
 INTERFACE LAP_heev  !--- legacy
    module procedure LAP_zheev_1,LAP_zheev_1_v
 end INTERFACE
 INTERFACE LAP__dsyev !--- legacy
    module procedure LAP_dsyev_1,LAP_dsyev_1_v
 end INTERFACE

 public :: LAP__ev, LAP__zheev_ini, LAP__zheev, LAP__zheevx, LAP__zheevr, LAP_heev, LAP__dsyev
 public :: LAP__ev_free

CONTAINS

 !-------------------------------------------------------
 ! zheev | single use ver.
 !  By setting 'JOBZ' explicitly, you can spare memory.
 !-------------------------------------------------------
 subroutine LAP_zheev_1(N,A,egVals,JOBZ,UPLO)
  integer, intent(in) :: N
  complex(kind=double), intent(inout)        :: A(N,N)
  real   (kind=double), intent(out)          :: egVals(N)
  character(LEN=1),     intent(in), optional :: JOBZ,  UPLO
  character(LEN=1)                           :: lJOBZ, lUPLO
  complex(kind=double), allocatable :: Aorg(:,:)
  logical :: saveA = .True.
  !--temp--
  complex(kind=double), allocatable :: zWork(:)  !--- for zheev(LWORK)
  real   (kind=double), allocatable :: dWork(:)  !--- for zheev(3*N-2)
  integer :: info
     allocate(zWork((64+1)*N),dWork(3*N-2))
     if( present(JOBZ) ) then
        lJOBZ = JOBZ
        saveA = .False.
     else ; lJOBZ = 'N'
     end if
     lUPLO = 'U'  ;  if( present(UPLO) ) lUPLO = UPLO
     if(saveA) then
        allocate(Aorg(N,N))
        Aorg = A
     end if
     call zheev(lJOBZ,lUPLO,N, A       , N ,egVals,zWork       ,(64+1)*N,dWork     ,info)
      !LAPACK  ("N/V","U/L",N,ZA(LDA,N),LDA,DW(N) ,ZWORK(LWORK),LWORK,DRWORK(3*N-2),INFO)
      !          in    in   in  inout   in  out    wout         in    w             out
     if(saveA) then
        A = Aorg
        deallocate(Aorg)
     end if
 end subroutine LAP_zheev_1
 !-------------------------------------------------------
 subroutine LAP_dsyev_1(N,A,egVals,JOBZ,UPLO)
  integer, intent(in) :: N
  real   (kind=double), intent(inout)        :: A(N,N)
  real   (kind=double), intent(out)          :: egVals(N)
  character(LEN=1),     intent(in), optional :: JOBZ,  UPLO
  character(LEN=1)                           :: lJOBZ, lUPLO
  real   (kind=double), allocatable :: Aorg(:,:)
  logical :: saveA = .True.
  !--temp--
  real   (kind=double), allocatable :: dWork(:)  !--- for zheev(3*N-2)
  integer :: info
     allocate(dWork((64+2)*N))
     if( present(JOBZ) ) then
        lJOBZ = JOBZ
        if(lJOBZ == 'N') saveA = .False.
     else ; lJOBZ = 'N'
     end if
     lUPLO = 'U'  ;  if( present(UPLO) ) lUPLO = UPLO
     if(saveA) then
        allocate(Aorg(N,N))
        Aorg = A
     end if
     call dsyev(lJOBZ,lUPLO,N, A       , N ,egVals,dWork      ,(64+2)*N,info)
      !LAPACK  ("N/V","U/L",N,DA(LDA,N),LDA,DW(N) ,WORK(LWORK),LWORK,INFO)
      !          in    in   in  inout   in  out    wout        in    out
     if(saveA) then
        A = Aorg
        deallocate(Aorg)
     end if
 end subroutine LAP_dsyev_1

 !-------------------------------------------------------
 ! zheev,dsyev | single use, egVecs output ver.
 !-------------------------------------------------------
 subroutine LAP_zheev_1_v(N,A,egVals,egVecs,UPLO)
  integer, intent(in) :: N
  complex(kind=double), intent(in)           :: A(N,N)
  complex(kind=double), intent(out)          :: egVecs(N,N)
  real   (kind=double), intent(out)          :: egVals(N)
  character(LEN=1),     intent(in), optional :: UPLO
  character(LEN=1)                           :: lUPLO
  !--temp--
  complex(kind=double), allocatable :: zWork(:)  !--- for zheev(LWORK)
  real   (kind=double), allocatable :: dWork(:)  !--- for zheev(3*N-2)
  integer :: info
     allocate(zWork((64+1)*N),dWork(3*N-2))
     lUPLO = 'U'  ;  if( present(UPLO) ) lUPLO = UPLO
     egVecs = A
     call zheev("V"  ,lUPLO,N, egVecs  , N ,egVals,zWork       ,(64+1)*N,dWork     ,info)
      !LAPACK  ("N/V","U/L",N,ZA(LDA,N),LDA,DW(N) ,ZWORK(LWORK),LWORK,DRWORK(3*N-2),INFO)
      !          in    in   in  inout   in  out    wout         in    w             out
 end subroutine LAP_zheev_1_v
 !-------------------------------------------------------
 subroutine LAP_dsyev_1_v(N,A,egVals,egVecs,UPLO)
  integer, intent(in) :: N
  real   (kind=double), intent(in)           :: A(N,N)
  real   (kind=double), intent(out)          :: egVecs(N,N)
  real   (kind=double), intent(out)          :: egVals(N)
  character(LEN=1),     intent(in), optional :: UPLO
  character(LEN=1)                           :: lUPLO
  !--temp--
  real   (kind=double), allocatable :: dWork(:)  !--- for zheev(3*N-2)
  integer :: info
     allocate(dWork((64+2)*N))
     lUPLO = 'U'  ;  if( present(UPLO) ) lUPLO = UPLO
     egVecs = A
     call dsyev("V"  ,lUPLO,N, egVecs  , N ,egVals,dWork      ,(64+2)*N,info)
      !LAPACK  ("N/V","U/L",N,DA(LDA,N),LDA,DW(N) ,WORK(LWORK),LWORK,INFO)
      !          in    in   in  inout   in  out    wout        in    out
 end subroutine LAP_dsyev_1_v

!-------------------------------------------------------------------------------------------------------



 !-------------------------------------------------------
 ! initialize
 !-------------------------------------------------------
 subroutine ini_A(A)
  complex(kind=double), intent(in) :: A(:,:)
    N = size(A,1)
    call this__ini
 end subroutine ini_A
 !-------------------------------------------------------
 subroutine ini_N(Nin)
  integer             , intent(in)  :: Nin
    N = Nin
    call this__ini
 end subroutine ini_N
 !-------------------------------------------------------
 subroutine this__ini
  complex(kind=double) :: Wvecs(1,1)
  real   (kind=double) :: Wvals(1)
  integer :: info
     call LAP__ev_free
     allocate (zWork(1),dWork(3*N-2))
     !-- calculate optimal LWORK
     !-- (When LWORK[IN]=-1, zWork(1)[OUT] returns optimal LWORK for )
     call zheev("N","U",N,Wvecs,N,Wvals,zWork, -1 ,dWork,info)
     LWORK = zWork(1)
     deallocate(zWork)
     allocate (zWork(LWORK))
 end subroutine this__ini
 !-------------------------------------------------------
 subroutine LAP__ev_free
     if(allocated(zWork)) deallocate(zWork)
     if(allocated(dWork)) deallocate(dWork)
 end subroutine

 !-------------------------------------------------------
 ! use zheev
 !-------------------------------------------------------
 subroutine LAP_zheev_run(A,egVals,JOBZ,UPLO)
  complex(kind=double), intent(inout)        :: A(N,N)
  real   (kind=double), intent(out)          :: egVals(N)
  character(LEN=1),     intent(in), optional :: JOBZ,  UPLO
  character(LEN=1)                           :: lJOBZ, lUPLO
  complex(kind=double), allocatable          :: Aorg(:,:)
  logical :: saveA = .True.
  integer :: info
     if( present(JOBZ) ) then
        lJOBZ = JOBZ
        saveA = .False.
     else ; lJOBZ = 'N'
     end if
     lUPLO = 'U'  ;  if( present(UPLO) ) lUPLO = UPLO
     if(saveA) then
        allocate(Aorg(N,N))
        Aorg = A
     end if
     call zheev(lJOBZ,lUPLO,N, A       , N ,egVals,zWork       ,LWORK,dWork        ,info)
      !LAPACK  ("N/V","U/L",N,ZA(LDA,N),LDA,DW(N) ,ZWORK(LWORK),LWORK,DRWORK(3*N-2),INFO)
       !          in    in   in  inout   in  out    wout         in    w             out
     if(saveA) then
        A = Aorg
        deallocate(Aorg)
     end if
 end subroutine LAP_zheev_run
 !-------------------------------------------------------


!-------------------------------------------------------------------------------------------------------

 !-------------------------------------------------------
 ! zheevx
 !---
 ! RANGE = 'A': all eigenvalues will be found.
 !       = 'V': all eigenvalues in the half-open interval (VL,VU]
 !              will be found.
 !       = 'I': the IL-th through IU-th eigenvalues will be found.
 !-------------------------------------------------------
 ! RANGE='I' version
 !---
 subroutine LAP_zheevx_I(A,Imin,Imax,egVals,egVecs,UPLO,highAccuracy)
  implicit none
  complex(kind=double), intent(inout)         :: A(:,:)
  integer             , intent(in)            :: Imin,Imax
  real   (kind=double)                        :: Vmin,Vmax
  real   (kind=double), intent(out)           :: egVals(:)
  complex(kind=double), intent(out), optional :: egVecs(:,:)
  character(LEN=1),     intent(in) , optional :: UPLO
  logical,              intent(in) , optional :: highAccuracy
  complex(kind=double), allocatable :: wVecs(:,:)
  complex(kind=double), allocatable :: AdiagOrg(:)
  integer             , allocatable :: iWork(:), iFail(:)
  real   (kind=double)              :: ABSTOL
  integer                           :: N, LzWork
  integer                           :: nEig
  integer                           :: i
  character(LEN=1)                  :: chUPLO
  real   (kind=double), external    :: DLAMCH
  !--temp--
  complex(kind=double), allocatable :: zWork(:)  !--- for zheev(LWORK)
  real   (kind=double), allocatable :: dWork(:)  !--- for zheev(3*N-2)
  integer :: info
   Vmin = -1.0d0
   Vmax = 1.0d0
   N = size(A,1)
   nEig = Imax-Imin+1
   if( present(egVecs) ) then
      if(size(egVecs,2) < nEig) then
         print *, "### zheevx error ### size of egVecs is too small."
         STOP
      end if
      allocate(wVecs(1,1))
   else
      allocate(wVecs(1,nEig))
   end if
   chUPLO = 'U'  ;  if( present(UPLO) ) chUPLO = UPLO
   ABSTOL = 0.0d0
   if( present(highAccuracy) ) then 
      if(highAccuracy) ABSTOL = 2*DLAMCH('S')  !--- twice the underflow threashold
   end if

   allocate (iFail(N),iWork(5*N),dWork(7*N))
   allocate (zWork(1))

   !--- calculate optimal LzWork
   !---  When LzWork[IN]=-1, zWork(1)[OUT] returns optimal LzWork (w.r.t. UPLO and N)
   LzWork = -1
   call zheevx('N'  ,'I'    ,chUPLO,N, A       , N ,Vmin,Vmax,Imin,Imax,ABSTOL, &
             & nEig,egVals,wVecs   ,N  ,zWork       ,LzWork,dWork     ,iWork     ,iFail   ,info)
   LzWork = zWork(1)
   deallocate(zWork)
   allocate (zWork(LzWork))
   !--- main part ---
   allocate (AdiagOrg(N))
   do i=1,N
      AdiagOrg(i) = A(i,i)
   end do
   if(present(egVecs)) then
      call zheevx('V'  ,'I'    ,chUPLO,N, A       , N ,Vmin,Vmax,Imin,Imax,ABSTOL, &
      !LAPACK    ("N/V","A/V/I","U/L" ,N,ZA(LDA,N),LDA,VL  ,VU  ,IL  ,IU  ,ABSTOL,
      !           in     in      in    in inout    in  in   in   in   in    in
                 & nEig,egVals,egVecs  ,N  ,zWork       ,LzWork,dWork     ,iWork     ,iFail   ,info)
      !            M   ,DW(N) ,Z(LDZ,M),LDZ,ZWORK(LWORK),LWORK,RWORK(7*N),IWORK(5*N),IFAIL(N),INFO)
      !            out   out    out     in   w(out)       in    w           w         out     out
   else
      call zheevx('N'  ,'I'    ,chUPLO,N, A       , N ,Vmin,Vmax,Imin,Imax,ABSTOL, &
                & nEig,egVals,wVecs   ,N  ,zWork       ,LzWork,dWork     ,iWork     ,iFail   ,info)
   end if
   call this__reconstructA(chUPLO,A,AdiagOrg)
   deallocate (AdiagOrg)
 end subroutine LAP_zheevx_I
!--------
! RANGE='V' version
!  !!! not confirmed yet
 subroutine LAP_zheevx_V(A,Vmin,Vmax,nEig,egVals,egVecs,UPLO,highAccuracy)
  implicit none
  complex(kind=double), intent(inout)         :: A(:,:)
  integer                                     :: Imin,Imax
  real   (kind=double), intent(in)            :: Vmin,Vmax
  integer             , intent(out)           :: nEig
  real   (kind=double), intent(out)           :: egVals(:)
  complex(kind=double), intent(out), optional :: egVecs(:,:)
  character(LEN=1),     intent(in) , optional :: UPLO
  logical,              intent(in) , optional :: highAccuracy
  complex(kind=double), allocatable :: wVecs(:,:)
  integer             , allocatable :: iWork(:), iFail(:)
  complex(kind=double), allocatable :: AdiagOrg(:)
  real   (kind=double)              :: ABSTOL
  integer                           :: N, LzWork
  integer                           :: i
  character(LEN=1)                  :: chUPLO
  real   (kind=double), external    :: DLAMCH
  !--temp--
  complex(kind=double), allocatable :: zWork(:)  !--- for zheev(LWORK)
  real   (kind=double), allocatable :: dWork(:)  !--- for zheev(3*N-2)
  integer :: info
   Imin = 1
   Imax = 2
   N = size(A,1)
   nEig = N !--- upper bound
   if( present(egVecs) ) then
      if(size(egVecs,2) < nEig) then
         print *, "### zheevx error ### size of egVecs is too small."
         STOP
      end if
      allocate(wVecs(1,1))
   else
      allocate(wVecs(1,nEig))
   end if
   chUPLO = 'U'
   if( present(UPLO) ) chUPLO = UPLO
   ABSTOL = 0.0d0
   if( present(highAccuracy) ) then 
      if(highAccuracy) ABSTOL = 2*DLAMCH('S')  !--- twice the underflow threashold
   end if

   allocate (iFail(N),iWork(5*N),dWork(7*N))
   allocate (zWork(1))

   !--- calculate optimal LzWork
   !---  When LzWork[IN]=-1, zWork(1)[OUT] returns optimal LzWork (w.r.t. UPLO and N)
   LzWork = -1
   call zheevx('N'  ,'V'    ,chUPLO,N, A       , N ,Vmin,Vmax,Imin,Imax,ABSTOL, &
             & nEig,egVals,wVecs   ,N  ,zWork       ,LzWork,dWork     ,iWork     ,iFail   ,info)
   LzWork = zWork(1)
   deallocate(zWork)
   allocate (zWork(LzWork))

   !--- main part ---
   allocate (AdiagOrg(N))
   do i=1,N
      AdiagOrg(i) = A(i,i)
   end do
   if(present(egVecs)) then
      call zheevx('V'  ,'V'    ,chUPLO,N, A       , N ,Vmin,Vmax,Imin,Imax,ABSTOL, &
      !LAPACK    ("N/V","A/V/I","U/L" ,N,ZA(LDA,N),LDA,VL  ,VU  ,IL  ,IU  ,ABSTOL,
      !            in    in      in    in inout    in  in   in   in   in    in
                 & nEig,egVals,egVecs  ,N  ,zWork       ,LzWork,dWork     ,iWork     ,iFail   ,info)
      !            M   ,DW(N) ,Z(LDZ,M),LDZ,ZWORK(LWORK),LWORK,RWORK(7*N),IWORK(5*N),IFAIL(N),INFO)
      !            out   out    out     in   w(out)       in    w           w         out     out
   else
      call zheevx('N'  ,'V'    ,chUPLO,N, A       , N ,Vmin,Vmax,Imin,Imax,ABSTOL, &
                & nEig,egVals,wVecs   ,N  ,zWork       ,LzWork,dWork     ,iWork     ,iFail   ,info)
   end if
   call this__reconstructA(chUPLO,A,AdiagOrg)
   deallocate (AdiagOrg)
 end subroutine LAP_zheevx_V


 !-------------------------------------------------------
 ! use zheevr
 !-------------------------------------------------------
 ! RANGE = 'A': all eigenvalues will be found.
 !       = 'V': all eigenvalues in the half-open interval (VL,VU]
 !              will be found.
 !       = 'I': the IL-th through IU-th eigenvalues will be found.
 !-------------------------------------------------------
 ! ISUPPZ : support of Z.
 !          I don't know how to use.
 !-------------------------------------------------------
 ! RANGE='I' version
 subroutine LAP_zheevr_I(A,Imin,Imax,egVals,egVecs,UPLO,highAccuracy)
  implicit none
  complex(kind=double), intent(inout)         :: A(:,:)
  integer             , intent(in)            :: Imin,Imax
  real   (kind=double)                        :: Vmin,Vmax
  real   (kind=double), intent(out)           :: egVals(:)
  complex(kind=double), intent(out), optional :: egVecs(:,:)
  character(LEN=1),     intent(in) , optional :: UPLO
  logical,              intent(in) , optional :: highAccuracy
  complex(kind=double), allocatable :: wVecs(:,:), zWork(:)
  real   (kind=double), allocatable :: dWork(:)
  integer             , allocatable :: iWork(:), iFail(:), iSuppZ(:)
  complex(kind=double), allocatable :: AdiagOrg(:)
  real   (kind=double)              :: ABSTOL
  integer                           :: N, LzWork, LdWork, LiWork
  integer                           :: nEig
  integer                           :: i
  character(LEN=1)                  :: chUPLO
  real   (kind=double), external    :: DLAMCH
  !--temp--
  integer :: info
   Vmin = -1.0d0
   Vmax = +1.0d0
   N = size(A,1)
   nEig = Imax-Imin+1
   if( present(egVecs) ) then
      if(size(egVecs,2) < nEig) then
         print *, "### zheevx error ### size of egVecs is too small."
         STOP
      end if
      allocate(wVecs(1,1))
   else
      allocate(wVecs(1,nEig))
   end if
   chUPLO = 'U'
   if( present(UPLO) ) chUPLO = UPLO
   ABSTOL = 0.0d0
   if( present(highAccuracy) ) then 
      if(highAccuracy) ABSTOL = 2*DLAMCH('S')  !--- twice the underflow threashold
   end if

   allocate (iFail(N),iSuppZ(2*nEig))
   allocate (zWork(1),dWork(1),iWork(1))

   !--- calculate optimal LWORK
   !---  When LWORK[IN]=-1, work(1)[OUT] returns optimal LWORK (w.r.t. UPLO and N)
   LzWork = -1
   LdWork = -1
   LiWork = -1
   call zheevr('N'  ,'I'    ,chUPLO,N, A       , N ,Vmin,Vmax,Imin,Imax,ABSTOL, &
              & nEig,egVals,wVecs   ,N  ,iSuppZ, &
              & zWork       ,LzWork,dWork        ,LdWork,iWork        ,LiWork,iFail   ,info)
   LzWork = zWork(1)
   LdWork = dWork(1)
   LiWork = iWork(1)
   deallocate(zWork,dWork,iWork)
   allocate (zWork(LzWork),dWork(LdWork),iWork(LiWork))
   !--- main part ---
   allocate (AdiagOrg(N))
   do i=1,N
      AdiagOrg(i) = A(i,i)
   end do
   if(present(egVecs)) then
      call zheevr('V'  ,'I'    ,chUPLO,N, A       , N ,Vmin,Vmax,Imin,Imax,ABSTOL, &
      !LAPACK    ("N/V","A/V/I","U/L" ,N,ZA(LDA,N),LDA,VL  ,VU  ,IL  ,IU  ,ABSTOL,
      !           in     in      in    in inout    in  in   in   in   in    in
                 & nEig,egVals,egVecs  ,N  ,iSuppZ, &
      !            M   ,DW(N) ,Z(LDZ,M),LDZ,ISUPPZ(2*M),
      !            out   out    out     in   out        
                 & zWork       ,LzWork,dWork        ,LdWork,iWork        ,LiWork,iFail   ,info)
      !            ZWORK(LWORK),LWORK ,RWORK(LRWORK),LRWORK,IWORK(LIWORK),LIWORK,IFAIL(N),INFO)
      !             w(out)       in     w(out)        in     w(out)        in     out     out
   else
      call zheevr('N'  ,'I'    ,chUPLO,N, A       , N ,Vmin,Vmax,Imin,Imax,ABSTOL, &
                 & nEig,egVals,wVecs   ,N  ,iSuppZ, &
                 & zWork       ,LzWork,dWork        ,LdWork,iWork        ,LiWork,iFail   ,info)
   end if
   call this__reconstructA(chUPLO,A,AdiagOrg)
   deallocate (AdiagOrg)
 end subroutine LAP_zheevr_I
 !-------------------------------------------------------
 ! RANGE='V' version
 subroutine LAP_zheevr_V(A,Vmin,Vmax,nEig,egVals,egVecs,UPLO,highAccuracy)
  implicit none
  complex(kind=double), intent(inout)         :: A(:,:)
  integer                                     :: Imin,Imax
  real   (kind=double), intent(in)            :: Vmin,Vmax
  integer             , intent(out)           :: nEig
  real   (kind=double), intent(out)           :: egVals(:)
  complex(kind=double), intent(out), optional :: egVecs(:,:)
  character(LEN=1),     intent(in) , optional :: UPLO
  logical,              intent(in) , optional :: highAccuracy
  complex(kind=double), allocatable :: zWork(:), wVecs(:,:)
  real   (kind=double), allocatable :: dWork(:)
  integer             , allocatable :: iWork(:), iFail(:), iSuppZ(:)
  complex(kind=double), allocatable :: AdiagOrg(:)
  real   (kind=double)              :: ABSTOL
  integer                           :: N, LzWork, LdWork, LiWork
  integer                           :: i
  character(LEN=1)                  :: chUPLO
  real   (kind=double), external    :: DLAMCH
  !--temp--
  integer :: info
   Imin = 1
   Imax = 2
   N = size(A,1)
   nEig = N !--- upper bound
   if( present(egVecs) ) then
      if(size(egVecs,2) < nEig) then
         print *, "### zheevx error ### size of egVecs is too small."
         STOP
      end if
      allocate(wVecs(1,1))
   else
      allocate(wVecs(1,nEig))
   end if
   chUPLO = 'U'
   if( present(UPLO) ) chUPLO = UPLO
   ABSTOL = 0.0d0
   if( present(highAccuracy) ) then 
      if(highAccuracy) ABSTOL = 2*DLAMCH('S')  !--- twice the underflow threashold
   end if

   allocate (iFail(N),iSuppZ(2*nEig))
   allocate (zWork(1),dWork(1),iWork(1))

   !--- calculate optimal LWORK
   !---  When LWORK[IN]=-1, work(1)[OUT] returns optimal LWORK (w.r.t. UPLO and N)
   LzWork = -1
   LdWork = -1
   LiWork = -1
   call zheevr('N'  ,'V'    ,chUPLO,N, A       , N ,Vmin,Vmax,Imin,Imax,ABSTOL, &
              & nEig,egVals,wVecs   ,N  ,iSuppZ, &
              & zWork       ,LzWork,dWork        ,LdWork,iWork        ,LiWork,iFail   ,info)
   LzWork = zWork(1)
   LdWork = dWork(1)
   LiWork = iWork(1)
   deallocate(zWork,dWork,iWork)
   allocate (zWork(LzWork),dWork(LdWork),iWork(LiWork))
   !--- main part ---
   allocate (AdiagOrg(N))
   do i=1,N
      AdiagOrg(i) = A(i,i)
   end do
   if(present(egVecs)) then
      call zheevr('V'  ,'V'    ,chUPLO,N, A       , N ,Vmin,Vmax,Imin,Imax,ABSTOL, &
      !LAPACK    ("N/V","A/V/I","U/L" ,N,ZA(LDA,N),LDA,VL  ,VU  ,IL  ,IU  ,ABSTOL,
      !           in     in      in    in inout    in  in   in   in   in    in
                 & nEig,egVals,egVecs  ,N  ,iSuppZ, &
      !            M   ,DW(N) ,Z(LDZ,M),LDZ,ISUPPZ(2*M),
      !            out   out    out     in   out        
                 & zWork       ,LzWork,dWork        ,LdWork,iWork        ,LiWork,iFail   ,info)
      !            ZWORK(LWORK),LWORK ,RWORK(LRWORK),LRWORK,IWORK(LIWORK),LIWORK,IFAIL(N),INFO)
      !             w(out)       in     w(out)        in     w(out)        in     out     out
   else
      call zheevr('N'  ,'V'    ,chUPLO,N, A       , N ,Vmin,Vmax,Imin,Imax,ABSTOL, &
                 & nEig,egVals,wVecs   ,N  ,iSuppZ, &
                 & zWork       ,LzWork,dWork        ,LdWork,iWork        ,LiWork,iFail   ,info)
   end if
   call this__reconstructA(chUPLO,A,AdiagOrg)
   deallocate (AdiagOrg)
 end subroutine LAP_zheevr_V

 !--- in LAPACK, 'UPLO' triangler part including diagonal is destroyed ---
 subroutine this__reconstructA(chUPLO,A,AdiagOrg)
  character(LEN=1)    , intent(in)    :: chUPLO
  complex(kind=double), intent(inout) :: A(:,:)
  complex(kind=double), intent(in)    :: AdiagOrg(:)
  integer :: i,j
   if(chUPLO=='U')then
    	do j=1,N
         A(j,j) = AdiagOrg(j)
         do i=j+1,N
            A(i,j) = A(j,i)
         end do
      end do
   else
    	do j=1,N
         do i=1,j-1
            A(i,j) = A(j,i)
         end do
         A(j,j) = AdiagOrg(j)
      end do
   end if
 end subroutine

!-------------------------------------------------------------------------------------------------------


end module LAP_ev
