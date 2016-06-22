program main
  implicit none
  integer, parameter :: dp=kind(0.d0)
  complex(dp), dimension(2,2) :: a
  complex(dp), dimension(4,4) :: b
  complex(dp), parameter :: ci = dcmplx(0.d0,1.d0)

  interface inverse
    subroutine inverse(n,m)
      integer, intent(in) :: n
      complex(kind=8), intent(inout) :: m(n,n)
    end subroutine inverse
  end interface

  a = transpose(reshape((/ complex(dp) :: &
     1, 2, &
     3, 4 /), shape(a)))

  b = transpose(reshape((/ complex(dp) :: &
    ci, 0, 0, 0, &
     0,-1, 0, 0, &
     0, 0, 1, 0, &
     0, 0, 0,-1 /), shape(b)))

  print *, "a="
  print *, a(1,1), a(1,2)
  print *, a(2,1), a(2,2)
  call inverse(2,a)
  print *, "a^-1="
  print *, a(1,1), a(1,2)
  print *, a(2,1), a(2,2)

  print *, "b="
  print *, b
  call inverse(4,b)
  print *, "b^-1="
  print *, b
end program main