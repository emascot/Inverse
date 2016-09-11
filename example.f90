program example
  use inverse
  implicit none
  integer, parameter :: n=200, dp=kind(0.d0)
  character(len=1), parameter :: mattype='h'
  integer :: i, size, time, ierr
  integer, allocatable :: seed(:)
  real(dp) :: re(n,n), im(n,n)
  complex(dp) :: m(n,n), minv(n,n)

  ! Use system time as random seed
  call random_seed(size = size)
  allocate(seed(size))
  call system_clock(count=time)
  seed = time + 37 * [(i-1, i=1, size)]
  call random_seed(put=seed)
  deallocate(seed)

  ! Generate random matrix
  call random_number(re)
  call random_number(im)
  m = dcmplx(re, im)

  ! Make symmetric, hermitian, or positive definite
  if (mattype .eq. 's') then
    m = 0.5 * (m + transpose(m))
  elseif (mattype .eq. 'h') then
    m = 0.5 * (m + conjg(transpose(m)))
  elseif (mattype .eq. 'p') then
    m = 0.5 * (m + conjg(transpose(m)))
    do i=1,n
      m(i,i) = m(i,i) + n
    enddo
  endif

  ! Print matrix
  write(*,*) "M"
  call print_matrix(n, m)

  ! Invert
  minv = m
  call invert(n, minv, 's', ierr)

  ! Print inverse
  write(*,*) "M^-1"
  call print_matrix(n, minv)

  ! Print identity
  write(*,*) "M.M^-1"
  call print_matrix(n, matmul(m, minv))
end program example


subroutine print_matrix(N,M)
  implicit none
  integer, parameter :: dp=kind(0.d0)
  integer, intent(in) :: N
  complex(dp), intent(in) :: M(N,N)
  complex(dp) :: MT(N,N)
  integer :: i, w, d, nmax
  character(len=64) :: fmt

  ! Show only corner if n is too large
  nmax = min(8, n)

  ! Use row major for proper printing
  MT = transpose(M)

  ! Set number of digits after decimal
  d = 5

  ! Calculate total width
  w = int(log10(maxval(abs([dreal(MT), dimag(MT), 1._dp]))))+3+d
  w = max(4+d, w)

  ! Construct format using width
  write(fmt,'("(",I0,A,I0,".",I0,A,I0,".",I0,A)') nmax,'(4x,ss,"(",f',w,d,',sp,f',w,d,',"i)"))'

  ! Print
  do i=1,nmax
    write(*,fmt) MT(1:nmax,i)
  enddo
  write(*,*)
end subroutine print_matrix