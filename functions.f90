module functions

  use bacs
  implicit none

contains

  function arg(z)
    implicit none
    complex(kind=double), intent(in) :: z
    real(kind=double) :: arg

    arg = atan2(aimag(z),real(z))

  end function arg

  function index(m,n,spin,Nx,Ny)
    implicit none
    integer, intent(in) :: m, n, Ny, Nx, spin
    integer :: index
    
    index = Nx*Ny*(spin-1)+m*Ny+n+1
    
  end function index

  function index2(m,n,spin,dagger,Nx,Ny)
    implicit none
    integer, intent(in) :: m, n, spin, dagger, Nx, Ny
    integer :: index2

    index2 = 2*Nx*Ny*(dagger-1)+Nx*Ny*(spin-1)+m*Ny+n+1

  end function index2

  function index3(n,alpha_1,alpha_n,q,lx,ly,Nx,Ny)
    implicit none
    integer,intent(in)::n,alpha_1,alpha_n,q,lx,ly,Nx,Ny
    integer::index3
    integer::ns
    integer::Nn

    ns=n-alpha_1
    Nn=alpha_n-alpha_1+1
    
    index3 = ns*(Nx/q)*Ny+(Nx/q)*ly+lx+1
  end function index3

  function indexB2(q,m,spin)
    implicit none
    integer,intent(in)::q,m,spin
    integer indexB2

    indexB2 = q*(spin-1)+m+1
  end function indexB2
  
  
  function indexB(q,m,spin,dagger)
    implicit none
    integer,intent(in)::q,m,spin,dagger
    integer::indexB

    indexB = 2*q*(dagger-1)+q*(spin-1)+m+1
  end function indexB

  function eh(dagger)
    implicit none
    integer,intent(in)::dagger
    real::eh

    eh = (-1)**(dagger+1)
  end function eh
  
  

  !---Charn number (Fukui Hatsugai Suzuki---
  !---end FHS formula---
  !---Chern number (NTW formula)---
  !---instructions---
  !---ひねった境界条件の位相theta_x,theta_y は，{2*pi/N}*n (n = 0,1, 2, ... N-1)とすること．
  !---argument---
  !---dim...Hamiltonianの次元
  !---m,n...状態u_m から状態u_m までのバンドのChern数を調べる
  !---msはメッシュの数
  !---固有ベクトルE(dim,dim,lx,ly)...二つ目の引数はバンド
  !---chern数を返す，実数
  
  subroutine NTW(dim,m,n,ms,multiplet,Chern) !(dim, band1, matrix, )
    implicit none
    
    integer, intent(in) :: dim, m, n, ms !m,n は，調べたいバンドの初めと終わり
    real(kind=double), intent(inout) :: Chern
    complex(kind=double), intent(in) :: multiplet(dim,n-m+1,0:ms-1,0:ms-1)
    integer :: i, j, k
    integer :: qN
    integer :: x1, x2
    integer :: y1, y2
    integer :: ierr
    real(kind=double) :: berry
    character(50) :: filename, chq, chb, chu

    qN = n - m + 1
    Chern = 0.0d0
    write(chq,'(i0)') dim
    write(chb,'(i0)') m
    write(chu,'(i0)') n
    
    do x1 = 0, ms-1
       x2 = mod(x1+1,ms)
       do y1 = 0, ms-1
          y2 = mod(y1+1,ms)

          berry = arg&
               (det(qN, matmul(transpose(conjg(multiplet(:,:,x1,y1))),multiplet(:,:,x2,y1)))&
               *det(qN, matmul(transpose(conjg(multiplet(:,:,x2,y1))),multiplet(:,:,x2,y2)))&
               *det(qN, matmul(transpose(conjg(multiplet(:,:,x2,y2))),multiplet(:,:,x1,y2)))&
               *det(qN, matmul(transpose(conjg(multiplet(:,:,x1,y2))),multiplet(:,:,x1,y1))))

          
          Chern = Chern + berry/twopi

       end do
    end do

  end subroutine NTW

  !---end NTW formula--- 


  !---determinant---
  Function det(dim, matrix)
  Implicit none

  ! output and input parameters                                             

  Integer, Intent(in) :: dim
  complex(kind(0d0)), Intent(in) :: matrix(1:dim, 1:dim)
  complex(kind(0d0)) :: det
  
  ! Constants

  Double precision, parameter :: pi = acos(-1.0d0)
  complex(kind(0d0)), Parameter :: runit = (1.0d0, 0.0d0), &
       iunit = (0.0d0, 1.0d0)

  ! working parameters                                                         

  Integer :: i

  ! lapack blas variables                                                     

  Integer :: m, n, lda, info
  Integer, Allocatable :: Ipiv(:)
  
  ! definition to call zgetrf
  
  m = dim
  n = dim
  lda = dim
  Allocate(Ipiv(1:n))

  Call zgetrf(m, n, matrix, lda, Ipiv, info)

  det = runit
  Do i = 1, dim
     det = det * matrix(i, i)
     If (Ipiv(i) /= i) det = -1.0d0 * det
  End Do

  DeAllocate(Ipiv)
End Function det
!----------------

!---Pauli Matrix---
subroutine Pauli_Matrix(P)
  implicit none
  complex(kind=double),allocatable,  intent(inout) :: P(:,:,:)
  allocate(P(3,2,2))
  
  P(1,1,1) = 0.0d0 + img*0.0d0
  P(1,1,2) = 1.0d0 + img*0.0d0
  P(1,2,1) = 1.0d0 + img*0.0d0
  P(1,2,2) = 0.0d0 + img*0.0d0
  
  P(2,1,1) = 0.0d0 + img*0.0d0
  P(2,1,2) = 0.0d0 - img*1.0d0
  P(2,2,1) = 0.0d0 + img*1.0d0
  P(2,2,2) = 0.0d0 + img*0.0d0
    
  P(3,1,1) = 1.0d0 + img*0.0d0
  P(3,1,2) = 0.0d0 + img*0.0d0
  P(3,2,1) = 0.0d0 + img*0.0d0
  P(3,2,2) = -1.0d0 + img*0.0d0
    
end subroutine Pauli_Matrix

function fd(epsilon,temperature)
  implicit none
  real(kind=double), intent(in) :: temperature,epsilon
  real(kind=double) :: fd                                 
  
  fd=1/(exp(epsilon/temperature)+1)
end function fd
!-------------------------------


!---diagonarize ZTEEVR---
subroutine diagonalize(N, H, E, Ev)
  implicit none
  integer, intent(in) :: N
  COMPLEX(kind=double), intent(inout) :: H(N, N)
  real(kind=double), intent(out) :: E(N)
  COMPLEX(kind=double), intent(out) :: Ev(N, N)

  CHARACTER(1) :: JOBZ, RANGE, UPLO
  real(kind=double) :: VL, VU, ABSTOL
  integer :: IL, IU, M, LDA, LDZ, INFO

  integer, allocatable :: ISUPPZ(:), IWORK(:)
  complex(kind=double), allocatable :: WORK(:)
  real(kind=double), allocatable :: RWORK(:)
  integer :: LWORK, LRWORK, LIWORK

  ! 行列サイズ設定
  LDA = N
  LDZ = N

  ! 計算モード設定
  JOBZ = 'V'      ! 固有値と固有ベクトルを計算
  RANGE = 'A'     ! 全ての固有値を求める
  UPLO = 'U'      ! 行列の上三角部分を使用
  VL = 0.0D0      ! (RANGE='A' の場合不要)
  
  VU = 0.0D0
  IL = 0
  IU = 0
  ABSTOL = 1.0D-8 ! 絶対許容誤差

  ! 作業領域サイズを問い合わせ
  LWORK = -1
  LRWORK = -1
  LIWORK = -1
  allocate(WORK(1), RWORK(1), IWORK(1), ISUPPZ(2*N))

  
  call ZHEEVR(JOBZ, RANGE, UPLO, N, H, LDA, VL, VU, IL, IU, ABSTOL, M, E, Ev, LDZ, ISUPPZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)

  ! 最適な作業領域サイズを取得
  LWORK = int(real(WORK(1)))
  LRWORK = int(RWORK(1))
  LIWORK = IWORK(1)
  
  ! 作業配列を再確保
  deallocate(WORK, RWORK, IWORK)
  allocate(WORK(LWORK), RWORK(LRWORK), IWORK(LIWORK))
  
  ! 再度 ZHEEVR を実行して対角化
  call ZHEEVR(JOBZ, RANGE, UPLO, N, H, LDA, VL, VU, IL, IU, ABSTOL, M, E, Ev, LDZ, ISUPPZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)
  
  ! エラー処理
  if (INFO /= 0) then
     print *, "Error in ZHEEVR: INFO =", INFO
  end if
  
  ! メモリ解放
  deallocate(WORK, RWORK, IWORK, ISUPPZ)
end subroutine diagonalize


!---matmul---
!---行列✖行列---
subroutine MM(A, B, C, m, n, k)
  implicit none
  integer, intent(in) :: m, n, k
  real(kind=double), intent(in) :: A(m,k), B(k,n)
  real(kind=double), intent(out) :: C(m,n)
  
  ! DGEMM の引数
  real(kind=double) :: alpha, beta
  character(1) :: transa, transb
  
  ! 係数
  alpha = 1.0d0
  beta  = 0.0d0
  transa = 'N'
  transb = 'N'
  
  ! LAPACK (BLAS) DGEMM の呼び出し
  call DGEMM(transa, transb, m, n, k, alpha, A, m, B, k, beta, C, m)
  
end subroutine MM

!---行列✖ベクトル---

subroutine MV(A, x, y, M, N)
  implicit none
  integer, intent(in) :: M, N
  real(kind=double), intent(in) :: A(M, N), x(N)
  real(kind=double), intent(out) :: y(M)

  ! DGEMV の引数
  real(kind=double) :: alpha, beta
  character(1) :: trans

  ! 係数
  alpha = 1.0d0  ! A * x のスケール係数
  beta  = 0.0d0  ! y をゼロクリア
  trans = 'N'    ! 転置なしで A を使用

  ! BLAS DGEMV の呼び出し: y = A * x
  call DGEMV(trans, M, N, alpha, A, M, x, 1, beta, y, 1)

end subroutine MV
!---------

!------------------------------------------------------------
! Module: random_mod
!  - subroutine random_seed(seed): 整数 seed を受け取って乱数シードを設定
!  - subroutine rand01(r):            0<=r<1 の一様乱数を返す
!------------------------------------------------------------

  subroutine rand(seed)
    ! 引数 seed をもとに random_seed を初期化
    implicit none
    integer, intent(in)        :: seed
    integer, allocatable       :: put_seed(:)
    integer                    :: n, i

    ! ランダムシード配列の必要長を取得
    call random_seed(size = n)

    ! シード配列を確保
    allocate(put_seed(n))

    ! シード配列を埋める（単純に seed+i-1 を使う例）
    do i = 1, n
      put_seed(i) = seed + (i - 1)
    end do

    ! intrinsic の random_seed を呼んで内部状態を設定
    call random_seed(put = put_seed)

    ! 解放
    deallocate(put_seed)
  end subroutine rand
  

  subroutine rand01(r)
    ! 0 <= r < 1 の一様乱数を返す
    real, intent(out) :: r
    call random_number(r)
  end subroutine rand01

end module functions


