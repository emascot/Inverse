# Inverse

Invert a double complex square matrix in Fortran

4x4 inverse adapted from [gluInvertMatrix](http://www.mesa3d.org/)

|Methods          |Description                      |
|-----------------|---------------------------------|
|invert           |Main driver                      |
|invert2x2        |Invert 2x2                       |
|invert3x3        |Invert 3x3                       |
|invert4x4        |Invert 4x4                       |
|invert_symmetric |Invert symmetric                 |
|invert_hermitian |Invert hermitian                 |
|invert_positive  |Invert positive definite         |
|invert_general   |Invert general matrix            |
|is_symmetric     |Return true if symmetric         |
|is_hermitian     |Return true if hermitian         |
|is_positive      |Return true if positive definite |

##### Usage

```fortran
use inverse
...
call invert(n, m, mattype, ierr) ! mattype and ierr are optional
if (is_symmetric(m)) print *, "m is symmetric"
```
|Parameter|Description                        |
|---------|-----------------------------------|
|N        |Dimension of matrix                |
|M        |Matrix to invert                   |
|mattype  |(Optional) Matrix type             |
|         |- 's': symmetric indefinite        |
|         |- 'h': hermitian indefinite        |
|         |- 'p': hermitian positive definite |
|ierr     |(Optional) Status                  |
|         |- 0: successful exit               |