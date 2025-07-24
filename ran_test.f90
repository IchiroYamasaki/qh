program test_random
  use bacs
  use functions
  implicit none
  
  integer :: seed, i
  real    :: x

  ! ── ① メインで一度だけシードを決定 ────────────
  seed = 101

  ! ── ② サブルーチン呼び出しで乱数種を初期化 ───
  call rand(seed)

  ! ── ③ 乱数生成例 ────────────────────────────
  do i = 1, 10
    call rand01(x)
    print '(F0.6)', x
  end do
end program test_random
  
