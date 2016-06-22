C      FORTRAN square matrix inversion by Eric Mascot
C
C      inverse(N,M)
C
C      N - Dimension of matrix
C        integer
C        intent(in)
C
C      M - Matrix to invert
C        complex(kind=8)
C        dimension(N,N)
C        intent(inout)
C
C============================================================
C
C     program example
C       implicit none
C       integer, parameter :: dp=kind(0.d0)
C       complex(dp), dimension(2,2) :: a
C       complex(dp), dimension(4,4) :: b
C       complex(dp), parameter :: ci = dcmplx(0.d0,1.d0)
C 
C       interface inverse
C         subroutine inverse(n,m)
C           integer, intent(in) :: n
C           complex(kind=8), intent(inout) :: m(n,n)
C         end subroutine inverse
C       end interface
C 
C       a = transpose(reshape((/ complex(dp) :: &
C         ci, 0, &
C          0, 1 /), shape(a)))
C 
C       b = transpose(reshape((/ complex(dp) :: &
C         ci, 0, 0, 0, &
C          0,-1, 0, 0, &
C          0, 0, 1, 0, &
C          0, 0, 0,-1 /), shape(b)))
C 
C       print *, "a="
C       print *, a
C       call inverse(2,a)
C       print *, "a^-1="
C       print *, a
C 
C       print *, "b="
C       print *, b
C       call inverse(4,b)
C       print *, "b^-1="
C       print *, b
C     end program example

      subroutine inverse(N, M)
      use iso_fortran_env
      implicit none

      integer, intent(in) :: N
      complex(kind=8), intent(inout) :: M(N,N)
      complex(kind=8), allocatable, dimension(:) :: work
      integer, allocatable, dimension(:) :: ipiv
      integer :: err

      interface inverse2x2
        subroutine inverse2x2(M, INFO)
          complex(kind=8), intent(inout) :: M(2,2)
          integer, intent(out) :: INFO
        end subroutine inverse2x2
      end interface

      interface inverse4x4
        subroutine inverse4x4(M, INFO)
          complex(kind=8), intent(inout) :: M(4,4)
          integer, intent(out) :: INFO
        end subroutine inverse4x4
      end interface

      interface zgetrf
        subroutine zgetrf(M, N, A, LDA, IPIV, INFO)
          integer, intent(in) :: LDA, M, N
          integer, intent(out) :: INFO, IPIV( * )
          complex(kind=8), intent(inout) :: A( LDA, * )
        end subroutine zgetrf
      end interface

      interface zgetri
        subroutine zgetri(N, A, LDA, IPIV, WORK, LWORK, INFO)
          integer, intent(in) :: LDA, LWORK, N, IPIV( * )
          integer, intent(out) :: INFO
          complex(kind=8), intent(inout) :: A( LDA, * )
          complex(kind=8), intent(out) :: WORK( * )
        end subroutine zgetri
      end interface

      if (N.eq.2) then
        call inverse2x2(M, err)
        if (err.ne.0) write(ERROR_UNIT,*) "Inverse failed"
      else if (N.eq.4) then
        call inverse4x4(M, err)
        if (err.ne.0) write(ERROR_UNIT,*) "Inverse failed"
      else
        allocate(work(N), ipiv(N), stat=err)
        if (err.ne.0) write(ERROR_UNIT,*) 
     +    "Work, ipiv: Allocation request denied"

        call zgetrf(N,N,M,N,ipiv,err)
        if (err.ne.0) write(ERROR_UNIT,*) "LU Decomposition failed"

        call zgetri(N,M,N,ipiv,work,N,err)
        if (err.ne.0) write(ERROR_UNIT,*) "Inverse failed"

        if (allocated(work)) deallocate(work, stat=err)
        if (err.ne.0) write(ERROR_UNIT,*)
     +    "Work: Deallocation request denied"

        if (allocated(ipiv)) deallocate(ipiv, stat=err)
        if (err.ne.0) write(ERROR_UNIT,*)
     +    "IPIV: Deallocation request denied"  
      end if
      end subroutine inverse


C      FORTRAN 2 by 2 matrix inversion by Eric Mascot
C
C      inverse2x2(m,info)
C
C      m - matrix to invert
C        complex(kind=8)
C        dimension(2,2)
C        intent(inout)
C        array of size 2,2
C
C      info - status of result
C        integer
C        intent(out)
C        0 - Normal exit
C        1 - Singular matrix
C
C============================================================
C
C      program example
C        implicit none
C        complex, dimension(2,2) :: a,b
C        integer :: status
C
C        complex, parameter :: ci = dcmplx(0.d0,1.d0)
C
C        a = (/ complex :: &
C          ci, 0, &
C           0,-1  /)
C
C        call inverse2x2(a,status)
C
C        if (status.eq.0) then
C          print *, a
C        else
C          print *, "Inversion failed"
C        end if
C
C      end program example

      subroutine inverse2x2(m, info)
      implicit none

      complex(kind=8), dimension(2,2), intent(inout) :: m
      integer, intent(out) :: info

      complex(kind=8) :: det
      complex(kind=8), dimension(2,2) :: inv

      inv(1,1) = m(2,2)
      inv(1,2) =-m(1,2)
      inv(2,1) =-m(2,1)
      inv(2,2) = m(1,1)

      det = m(1,1) * m(2,2) - m(1,2) * m(2,1)

      if (det.eq.0) then
          info = 1
          return
      end if

      det = 1 / det

      m(1,1) = inv(1,1) * det
      m(1,2) = inv(1,2) * det
      m(2,1) = inv(2,1) * det
      m(2,2) = inv(2,2) * det

      info = 0
      end subroutine inverse2x2


C      FORTRAN implementation of gluInvertMatrix by Eric Mascot
C      http://www.mesa3d.org/
C
C      Inverts a 4 by 4 matrix
C
C
C      interface inverse4x4(m, info)
C        subroutine inverse4x4
C          complex(kind=8), intent(inout) :: m(4,4)
C          integer :: info
C        end subroutine inverse4x4
C      end interface
C
C      m - matrix to invert
C        complex(kind=8)
C        dimension(4,4)
C        intent(inout)
C        array of size 4,4
C
C      info - status of result
C        integer
C        intent(out)
C        0 - Normal exit
C        1 - Singular matrix
C
C============================================================
C
C      program example
C        implicit none
C        complex, dimension(4,4) :: a
C        integer :: status
C
C        complex, parameter :: ci = dcmplx(0.d0,1.d0)
C
C        a = (/ complex :: &
C          ci, 0, 0, 0, &
C           0,-1, 0, 0, &
C           0, 0, 1, 0, &
C           0, 0, 0,-1 /)
C
C        call inverse4x4(a,status)
C
C        if (status.eq.0) then
C          print *, a
C        else
C          print *, "Inversion failed"
C        end if
C
C      end program example

      subroutine inverse4x4(m, info)
      implicit none

      complex(kind=8), dimension(4,4), intent(inout) :: m
      integer, intent(out) :: info

      integer :: i,j
      complex(kind=8) :: det
      complex(kind=8), dimension(4,4) :: inv

      inv(1,1) = m(2,2) * m(3,3) * m(4,4) -
     +           m(2,2) * m(3,4) * m(4,3) -
     +           m(3,2) * m(2,3) * m(4,4) +
     +           m(3,2) * m(2,4) * m(4,3) +
     +           m(4,2) * m(2,3) * m(3,4) -
     +           m(4,2) * m(2,4) * m(3,3)

      inv(2,1) =-m(2,1) * m(3,3) * m(4,4) +
     +           m(2,1) * m(3,4) * m(4,3) +
     +           m(3,1) * m(2,3) * m(4,4) -
     +           m(3,1) * m(2,4) * m(4,3) -
     +           m(4,1) * m(2,3) * m(3,4) +
     +           m(4,1) * m(2,4) * m(3,3)

      inv(3,1) = m(2,1) * m(3,2) * m(4,4) -
     +           m(2,1) * m(3,4) * m(4,2) -
     +           m(3,1) * m(2,2) * m(4,4) +
     +           m(3,1) * m(2,4) * m(4,2) +
     +           m(4,1) * m(2,2) * m(3,4) -
     +           m(4,1) * m(2,4) * m(3,2)

      inv(4,1) =-m(2,1) * m(3,2) * m(4,3) +
     +           m(2,1) * m(3,3) * m(4,2) +
     +           m(3,1) * m(2,2) * m(4,3) -
     +           m(3,1) * m(2,3) * m(4,2) -
     +           m(4,1) * m(2,2) * m(3,3) +
     +           m(4,1) * m(2,3) * m(3,2)

      inv(1,2) =-m(1,2) * m(3,3) * m(4,4) +
     +           m(1,2) * m(3,4) * m(4,3) +
     +           m(3,2) * m(1,3) * m(4,4) -
     +           m(3,2) * m(1,4) * m(4,3) -
     +           m(4,2) * m(1,3) * m(3,4) +
     +           m(4,2) * m(1,4) * m(3,3)

      inv(2,2) = m(1,1) * m(3,3) * m(4,4) -
     +           m(1,1) * m(3,4) * m(4,3) -
     +           m(3,1) * m(1,3) * m(4,4) +
     +           m(3,1) * m(1,4) * m(4,3) +
     +           m(4,1) * m(1,3) * m(3,4) -
     +           m(4,1) * m(1,4) * m(3,3)

      inv(3,2) =-m(1,1) * m(3,2) * m(4,4) +
     +           m(1,1) * m(3,4) * m(4,2) +
     +           m(3,1) * m(1,2) * m(4,4) -
     +           m(3,1) * m(1,4) * m(4,2) -
     +           m(4,1) * m(1,2) * m(3,4) +
     +           m(4,1) * m(1,4) * m(3,2)

      inv(4,2) = m(1,1) * m(3,2) * m(4,3) -
     +           m(1,1) * m(3,3) * m(4,2) -
     +           m(3,1) * m(1,2) * m(4,3) +
     +           m(3,1) * m(1,3) * m(4,2) +
     +           m(4,1) * m(1,2) * m(3,3) -
     +           m(4,1) * m(1,3) * m(3,2)

      inv(1,3) = m(1,2) * m(2,3) * m(4,4) -
     +           m(1,2) * m(2,4) * m(4,3) -
     +           m(2,2) * m(1,3) * m(4,4) +
     +           m(2,2) * m(1,4) * m(4,3) +
     +           m(4,2) * m(1,3) * m(2,4) -
     +           m(4,2) * m(1,4) * m(2,3)

      inv(2,3) =-m(1,1) * m(2,3) * m(4,4) +
     +           m(1,1) * m(2,4) * m(4,3) +
     +           m(2,1) * m(1,3) * m(4,4) -
     +           m(2,1) * m(1,4) * m(4,3) -
     +           m(4,1) * m(1,3) * m(2,4) +
     +           m(4,1) * m(1,4) * m(2,3)

      inv(3,3) = m(1,1) * m(2,2) * m(4,4) -
     +           m(1,1) * m(2,4) * m(4,2) -
     +           m(2,1) * m(1,2) * m(4,4) +
     +           m(2,1) * m(1,4) * m(4,2) +
     +           m(4,1) * m(1,2) * m(2,4) -
     +           m(4,1) * m(1,4) * m(2,2)

      inv(4,3) =-m(1,1) * m(2,2) * m(4,3) +
     +           m(1,1) * m(2,3) * m(4,2) +
     +           m(2,1) * m(1,2) * m(4,3) -
     +           m(2,1) * m(1,3) * m(4,2) -
     +           m(4,1) * m(1,2) * m(2,3) +
     +           m(4,1) * m(1,3) * m(2,2)

      inv(1,4) =-m(1,2) * m(2,3) * m(3,4) +
     +           m(1,2) * m(2,4) * m(3,3) +
     +           m(2,2) * m(1,3) * m(3,4) -
     +           m(2,2) * m(1,4) * m(3,3) -
     +           m(3,2) * m(1,3) * m(2,4) +
     +           m(3,2) * m(1,4) * m(2,3)

      inv(2,4) = m(1,1) * m(2,3) * m(3,4) -
     +           m(1,1) * m(2,4) * m(3,3) -
     +           m(2,1) * m(1,3) * m(3,4) +
     +           m(2,1) * m(1,4) * m(3,3) +
     +           m(3,1) * m(1,3) * m(2,4) -
     +           m(3,1) * m(1,4) * m(2,3)

      inv(3,4) =-m(1,1) * m(2,2) * m(3,4) +
     +           m(1,1) * m(2,4) * m(3,2) +
     +           m(2,1) * m(1,2) * m(3,4) -
     +           m(2,1) * m(1,4) * m(3,2) -
     +           m(3,1) * m(1,2) * m(2,4) +
     +           m(3,1) * m(1,4) * m(2,2)

      inv(4,4) = m(1,1) * m(2,2) * m(3,3) -
     +           m(1,1) * m(2,3) * m(3,2) -
     +           m(2,1) * m(1,2) * m(3,3) +
     +           m(2,1) * m(1,3) * m(3,2) +
     +           m(3,1) * m(1,2) * m(2,3) -
     +           m(3,1) * m(1,3) * m(2,2)

      det = m(1,1) * inv(1,1) + m(1,2) * inv(2,1) +
     +      m(1,3) * inv(3,1) + m(1,4) * inv(4,1)

      if (det.eq.0) then
          info = 1
          return
      end if

      det = 1 / det

      do 10 i = 1,4
      do 20 j = 1,4
          m(i,j) = inv(i,j) * det
20    continue
10    continue

      info = 0
      end subroutine inverse4x4