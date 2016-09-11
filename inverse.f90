module inverse
  use iso_fortran_env
  implicit none
  integer, parameter :: dp=kind(0.0d0)
  private
  public :: invert, invert2x2, invert4x4, &
            invert_symmetric, invert_hermitian, &
            invert_positive, invert_general, &
            is_symmetric, is_hermitian, is_positive
contains

!==============================================================
!  invert
!> Invert a complex*16 square matrix
!> \param[in]      N        Dimension of matrix
!> \param[in,out]  M        Matrix to invert
!> \param[in]      mattype  Matrix type
!>   - 's': symmetric indefinite
!>   - 'h': hermitian indefinite
!>   - 'p': hermitian positive definite
!> \param[out]     ierr     Status of inversion
!> \date September 2016
!> \author Eric Mascot

  subroutine invert(N, M, mattype, ierr)
    implicit none
    integer, intent(in) :: N
    complex(dp), intent(inout) :: M(N,N)
    character(len=1), optional, intent(in) :: mattype
    integer, optional, intent(out) :: ierr
    integer :: ierr_


    if (N.eq.2) then
      call invert2x2(M, ierr_)
    elseif (N.eq.3) then
      call invert3x3(M, ierr_)
    elseif (N.eq.4) then
      call invert4x4(M, ierr_)
    else
      if (present(mattype)) then
        select case(mattype)
        case('s')
          call invert_symmetric(N, M, ierr_)
        case('h')
          call invert_hermitian(N, M, ierr_)
        case('p')
          call invert_positive(N, M, ierr_)
        case default
          call invert_general(N, M, ierr_)
        endselect
      else
        call invert_general(N, M, ierr_)
      endif
    endif

    if (present(ierr)) then
      ierr = ierr_
    endif
  end subroutine invert

  subroutine invert2x2(m, ierr)
    implicit none
    complex(dp), intent(inout) :: m(2,2)
    integer, intent(out) :: ierr
    integer :: i,j
    complex(dp) :: det
    complex(dp), dimension(2,2) :: inv

    inv(1,1) = m(2,2)
    inv(1,2) =-m(1,2)
    inv(2,1) =-m(2,1)
    inv(2,2) = m(1,1)

    det = m(1,1) * m(2,2) - m(1,2) * m(2,1)

    if (det.eq.0) then
      ierr = 1
      write(ERROR_UNIT,*) 'Invert failed'
      return
    end if

    det = 1._dp / det

    do i = 1,2
      do j = 1,2
        m(i,j) = inv(i,j) * det
      enddo
    enddo

    ierr = 0
  end subroutine invert2x2

  subroutine invert3x3(m, ierr)
    implicit none
    complex(dp), intent(inout) :: m(3,3)
    integer, intent(out) :: ierr
    complex(dp) :: det, inv(3,3)
    integer :: i, j

    inv(1,1) = m(2,2)*m(3,3) - m(3,2)*m(2,3)
    inv(1,2) = m(3,2)*m(1,3) - m(1,2)*m(3,3)
    inv(1,3) = m(1,2)*m(2,3) - m(1,3)*m(2,2)
    inv(2,1) = m(2,3)*m(3,1) - m(2,1)*m(3,3)
    inv(2,2) = m(1,1)*m(3,3) - m(3,1)*m(1,3)
    inv(2,3) = m(2,1)*m(1,3) - m(1,1)*m(2,3)
    inv(3,1) = m(2,1)*m(3,2) - m(2,2)*m(3,1)
    inv(3,2) = m(3,1)*m(1,2) - m(1,1)*m(3,2)
    inv(3,3) = m(1,1)*m(2,2) - m(1,2)*m(2,1)

    det = inv(1,1)*m(1,1) + inv(1,2)*m(2,1) + inv(1,3)*m(3,1)

    if (det.eq.0) then
      ierr = 1
      write(ERROR_UNIT,*) 'Invert failed'
      return
    endif

    det = 1._dp / det

    do i = 1,3
      do j = 1,3
        m(i,j) = inv(i,j) * det
      enddo
    enddo

    ierr = 0
  end subroutine invert3x3

  subroutine invert4x4(m, ierr)
    implicit none
    complex(dp), intent(inout) :: m(4,4)
    integer, intent(out) :: ierr
    integer :: i,j
    complex(dp) :: det, inv(4,4)

    inv(1,1) = m(2,2) * m(3,3) * m(4,4) - &
               m(2,2) * m(3,4) * m(4,3) - &
               m(3,2) * m(2,3) * m(4,4) + &
               m(3,2) * m(2,4) * m(4,3) + &
               m(4,2) * m(2,3) * m(3,4) - &
               m(4,2) * m(2,4) * m(3,3)

    inv(2,1) =-m(2,1) * m(3,3) * m(4,4) + &
               m(2,1) * m(3,4) * m(4,3) + &
               m(3,1) * m(2,3) * m(4,4) - &
               m(3,1) * m(2,4) * m(4,3) - &
               m(4,1) * m(2,3) * m(3,4) + &
               m(4,1) * m(2,4) * m(3,3)

    inv(3,1) = m(2,1) * m(3,2) * m(4,4) - &
               m(2,1) * m(3,4) * m(4,2) - &
               m(3,1) * m(2,2) * m(4,4) + &
               m(3,1) * m(2,4) * m(4,2) + &
               m(4,1) * m(2,2) * m(3,4) - &
               m(4,1) * m(2,4) * m(3,2)

    inv(4,1) =-m(2,1) * m(3,2) * m(4,3) + &
               m(2,1) * m(3,3) * m(4,2) + &
               m(3,1) * m(2,2) * m(4,3) - &
               m(3,1) * m(2,3) * m(4,2) - &
               m(4,1) * m(2,2) * m(3,3) + &
               m(4,1) * m(2,3) * m(3,2)

    inv(1,2) =-m(1,2) * m(3,3) * m(4,4) + &
               m(1,2) * m(3,4) * m(4,3) + &
               m(3,2) * m(1,3) * m(4,4) - &
               m(3,2) * m(1,4) * m(4,3) - &
               m(4,2) * m(1,3) * m(3,4) + &
               m(4,2) * m(1,4) * m(3,3)

    inv(2,2) = m(1,1) * m(3,3) * m(4,4) - &
               m(1,1) * m(3,4) * m(4,3) - &
               m(3,1) * m(1,3) * m(4,4) + &
               m(3,1) * m(1,4) * m(4,3) + &
               m(4,1) * m(1,3) * m(3,4) - &
               m(4,1) * m(1,4) * m(3,3)

    inv(3,2) =-m(1,1) * m(3,2) * m(4,4) + &
               m(1,1) * m(3,4) * m(4,2) + &
               m(3,1) * m(1,2) * m(4,4) - &
               m(3,1) * m(1,4) * m(4,2) - &
               m(4,1) * m(1,2) * m(3,4) + &
               m(4,1) * m(1,4) * m(3,2)

    inv(4,2) = m(1,1) * m(3,2) * m(4,3) - &
               m(1,1) * m(3,3) * m(4,2) - &
               m(3,1) * m(1,2) * m(4,3) + &
               m(3,1) * m(1,3) * m(4,2) + &
               m(4,1) * m(1,2) * m(3,3) - &
               m(4,1) * m(1,3) * m(3,2)

    inv(1,3) = m(1,2) * m(2,3) * m(4,4) - &
               m(1,2) * m(2,4) * m(4,3) - &
               m(2,2) * m(1,3) * m(4,4) + &
               m(2,2) * m(1,4) * m(4,3) + &
               m(4,2) * m(1,3) * m(2,4) - &
               m(4,2) * m(1,4) * m(2,3)

    inv(2,3) =-m(1,1) * m(2,3) * m(4,4) + &
               m(1,1) * m(2,4) * m(4,3) + &
               m(2,1) * m(1,3) * m(4,4) - &
               m(2,1) * m(1,4) * m(4,3) - &
               m(4,1) * m(1,3) * m(2,4) + &
               m(4,1) * m(1,4) * m(2,3)

    inv(3,3) = m(1,1) * m(2,2) * m(4,4) - &
               m(1,1) * m(2,4) * m(4,2) - &
               m(2,1) * m(1,2) * m(4,4) + &
               m(2,1) * m(1,4) * m(4,2) + &
               m(4,1) * m(1,2) * m(2,4) - &
               m(4,1) * m(1,4) * m(2,2)

    inv(4,3) =-m(1,1) * m(2,2) * m(4,3) + &
               m(1,1) * m(2,3) * m(4,2) + &
               m(2,1) * m(1,2) * m(4,3) - &
               m(2,1) * m(1,3) * m(4,2) - &
               m(4,1) * m(1,2) * m(2,3) + &
               m(4,1) * m(1,3) * m(2,2)

    inv(1,4) =-m(1,2) * m(2,3) * m(3,4) + &
               m(1,2) * m(2,4) * m(3,3) + &
               m(2,2) * m(1,3) * m(3,4) - &
               m(2,2) * m(1,4) * m(3,3) - &
               m(3,2) * m(1,3) * m(2,4) + &
               m(3,2) * m(1,4) * m(2,3)

    inv(2,4) = m(1,1) * m(2,3) * m(3,4) - &
               m(1,1) * m(2,4) * m(3,3) - &
               m(2,1) * m(1,3) * m(3,4) + &
               m(2,1) * m(1,4) * m(3,3) + &
               m(3,1) * m(1,3) * m(2,4) - &
               m(3,1) * m(1,4) * m(2,3)

    inv(3,4) =-m(1,1) * m(2,2) * m(3,4) + &
               m(1,1) * m(2,4) * m(3,2) + &
               m(2,1) * m(1,2) * m(3,4) - &
               m(2,1) * m(1,4) * m(3,2) - &
               m(3,1) * m(1,2) * m(2,4) + &
               m(3,1) * m(1,4) * m(2,2)

    inv(4,4) = m(1,1) * m(2,2) * m(3,3) - &
               m(1,1) * m(2,3) * m(3,2) - &
               m(2,1) * m(1,2) * m(3,3) + &
               m(2,1) * m(1,3) * m(3,2) + &
               m(3,1) * m(1,2) * m(2,3) - &
               m(3,1) * m(1,3) * m(2,2)

    det = m(1,1) * inv(1,1) + m(1,2) * inv(2,1) + &
          m(1,3) * inv(3,1) + m(1,4) * inv(4,1)

    if (det.eq.0) then
      ierr = 1
      write(ERROR_UNIT,*) 'Invert failed'
      return
    end if

    det = 1._dp / det

    do i = 1,4
      do j = 1,4
        m(i,j) = inv(i,j) * det
      enddo
    enddo

    ierr = 0
  end subroutine invert4x4

  ! symmetric indefinite
  subroutine invert_symmetric(n,m,ierr)
    implicit none
    integer, intent(in) :: n
    complex(dp), intent(inout) :: m(n,n)
    integer, intent(out) :: ierr
    complex(dp) :: work(2*n)
    integer :: ipiv(n), i, j
    
    call zsytrf('U',n,m,n,ipiv,work,n,ierr)
    if (ierr.ne.0) write(ERROR_UNIT,*) &
      'Upper triangular factorization failed'

    call zsytri('U',n,m,n,ipiv,work,ierr)
    if (ierr.ne.0) write(ERROR_UNIT,*) &
      'Invert failed'

    do i=1,n
      do j=1,i-1
        m(i,j) = m(j,i)
      enddo
    enddo
  end subroutine invert_symmetric

  ! hermitian indefinite
  subroutine invert_hermitian(n,m,ierr)
    implicit none
    integer, intent(in) :: n
    complex(dp), intent(inout) :: m(n,n)
    integer, intent(out) :: ierr
    complex(dp) :: work(n)
    integer :: ipiv(n), i, j
    
    call zhetrf('U',n,m,n,ipiv,work,n,ierr)
    if (ierr.ne.0) write(ERROR_UNIT,*) &
      'Upper triangular factorization failed'

    call zhetri('U',n,m,n,ipiv,work,ierr)
    if (ierr.ne.0) write(ERROR_UNIT,*) &
      'Invert failed'

    do i=1,n
      do j=1,i-1
        m(i,j) = dconjg(m(j,i))
      enddo
    enddo
  end subroutine invert_hermitian

  ! hermitian positive definite
  subroutine invert_positive(n,m,ierr)
    implicit none
    integer, intent(in) :: n
    complex(dp), intent(inout) :: m(n,n)
    integer, intent(out) :: ierr
    integer :: i, j

    call zpotrf('U',n,m,n,ierr)
    if (ierr.ne.0) write(ERROR_UNIT,*) &
      'Upper triangular factorization failed'

    call zpotri('U',n,m,n,ierr)
    if (ierr.ne.0) write(ERROR_UNIT,*) &
      'Invert failed'

    do i=1,n
      do j=1,i-1
        m(i,j) = dconjg(m(j,i))
      enddo
    enddo
  end subroutine invert_positive

  ! general - LU factorization
  subroutine invert_general(n,m,ierr)
    implicit none
    integer, intent(in) :: n
    complex(dp), dimension(n,n), intent(inout) :: m
    integer, intent(out) :: ierr
    complex(dp) :: work(n)
    integer :: ipiv(n)

    call zgetrf(n,n,m,n,ipiv,ierr)
    if (ierr.ne.0) write(ERROR_UNIT,*) &
      'LU Decomposition failed'

    call zgetri(n,m,n,ipiv,work,n,ierr)
    if (ierr.ne.0) write(ERROR_UNIT,*) &
      'Invert failed'
  end subroutine invert_general

  ! Test if a matrix is symmetric
  function is_symmetric(n,m)
    implicit none
    integer, intent(in) :: n
    complex(dp), intent(in) :: m(n,n)
    complex(dp) :: mt(n,n)
    logical :: is_symmetric

    mt = transpose(m)

    is_symmetric = all(abs(mt-m) .le. 1e-8)
  end function is_symmetric

  ! Test if a matrix is hermitian
  function is_hermitian(n,m)
    implicit none
    integer, intent(in) :: n
    complex(dp), intent(in) :: m(n,n)
    complex(dp) :: mt(n,n)
    logical :: is_hermitian

    mt = conjg(transpose(m))

    is_hermitian = all(abs(mt-m) .le. 1e-8)
  end function is_hermitian

  ! Test if a matrix is positive definite
  function is_positive(n,m)
    implicit none
    integer, intent(in) :: n
    complex(dp), intent(in) :: m(n,n)
    complex(dp) :: mp(n,n)
    integer :: ierr
    logical :: is_positive

    ! Duplicate
    mp = m

    ! Cholesky factorization
    call zpotrf('U',n,mp,n,ierr)
    if (ierr.ne.0) then
      is_positive = .false.
    else
      is_positive = .true.
    endif
  end function is_positive
end module inverse