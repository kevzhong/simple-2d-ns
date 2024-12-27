subroutine write2DField(array,Nx,Ny,fpref,it)
    ! Write array(Nx,Ny) field as binary data
    implicit none
    integer :: Nx, Ny
    character(len=50) :: filename 
    real, dimension(Nx,Ny) :: array
    integer :: it
    character(len=*) :: fpref

    write(filename, '("outputdir/",A, "_it", I0, ".dat")') trim(fpref), it
    open(unit=10, file=filename, status='replace', access='stream', action='write')

    write(10) array
    close(10)

end subroutine write2DField


subroutine writeVTK(array, Nx, Ny, dx, dy,fpref,it)
    ! Write array(Nx,Ny) field as binary .vtk file for convenient viewing in Paraview
    implicit none

    character(len=50) :: filename  
    character(len=*), intent(in) :: fpref  
    integer, intent(in) :: Nx, Ny
    real(8), intent(in) :: array(Nx,Ny)   
    real(8) :: VAL 
    real, intent(in) :: dx, dy  
    integer, intent(in) :: it

    ! Local variables
    integer :: i, j, num_points
    integer :: file_unit


    write(filename, '("outputdir/", A, "_it", I0, ".vtk")') trim(fpref), it


    ! Calculate the total number of points
    num_points = Nx * Ny

    open(newunit=file_unit, file=filename, status='replace', form='formatted', action='write')

    ! Write VTK ASCII header
    write(file_unit, '(A)') '# vtk DataFile Version 3.0'
    write(file_unit, '(A)') '2D structured grid data'
    !write(file_unit, '(A)') 'ASCII'
    write(file_unit, '(A)') 'BINARY'
    write(file_unit, '(A)') 'DATASET STRUCTURED_POINTS'
    write(file_unit, '(A, I0, A, I0, A)') 'DIMENSIONS ', Nx, ' ', Ny, ' 1'
    write(file_unit, '(A, F10.5, A, F10.5, A)') 'SPACING ', dx, ' ', dy, ' 1.0'
    write(file_unit, '(A)') 'ORIGIN 0.0 0.0 0.0'
    write(file_unit, '(A, I0)') 'POINT_DATA ', num_points
    write(file_unit, '(A)') 'SCALARS field double'
    write(file_unit, '(A)') 'LOOKUP_TABLE default'

    ! BINARY WRITING
    close(file_unit)
    ! Open the file in unformatted (binary) mode to append the data
    open(newunit=file_unit, file=filename, position='append', form='unformatted', access='stream')
    do j = 1, Ny
        do i = 1, Nx
            write(file_unit) array(i, j)
        end do
    end do


    ! ! ASCII WRITING
    ! DO j = 1, Ny
    !     DO i = 1, Nx
    !        VAL = array(i, j)
    !        WRITE(file_unit, '(F10.5)', ADVANCE='NO') VAL
    !        WRITE(file_unit, '(A)', ADVANCE='YES') ''  ! New line for clarity
    !     END DO
    !  END DO


    ! Close the file
    close(file_unit)

end subroutine writeVTK



! SUBROUTINE WriteVTK(filename, matrix, format, aratio)
!     IMPLICIT NONE
!     CHARACTER(LEN=*), INTENT(IN) :: filename
!     REAL, INTENT(IN), ALLOCATABLE :: matrix(:,:,:)
!     CHARACTER(LEN=*), INTENT(IN) :: format
!     REAL, INTENT(IN) :: aratio(3)
!     INTEGER :: Nx, Ny, Nz
!     INTEGER :: fid, i, j, k
!     REAL :: value
!     INTEGER :: total_points

!     ! Get dimensions
!     Nx = SIZE(matrix, 1)
!     Ny = SIZE(matrix, 2)
!     Nz = SIZE(matrix, 3)
!     total_points = Nx * Ny * Nz

!     ! Open file
!     OPEN(UNIT=10, FILE=filename, STATUS='REPLACE', FORM='FORMATTED', ACTION='WRITE')
!     fid = 10

!     ! Write header
!     WRITE(fid, '(A)') '# vtk DataFile Version 2.0'
!     WRITE(fid, '(A)') 'Volume example'
!     IF (TRIM(format) == 'ascii') THEN
!        WRITE(fid, '(A)') 'ASCII'
!     ELSEIF (TRIM(format) == 'binary') THEN
!        WRITE(fid, '(A)') 'BINARY'
!     ELSE
!        STOP 'Unsupported format'
!     END IF
!     WRITE(fid, '(A)') 'DATASET STRUCTURED_POINTS'
!     WRITE(fid, '(A,I0,A,I0,A,I0)') 'DIMENSIONS ', Nx, ' ', Ny, ' ', Nz
!     WRITE(fid, '(A,F0.1,A,F0.1,A,F0.1)') 'ASPECT_RATIO ', aratio(1), ' ', aratio(2), ' ', aratio(3)
!     WRITE(fid, '(A,I0,A,I0,A,I0)') 'ORIGIN ', 0, ' ', 0, ' ', 0
!     WRITE(fid, '(A,I0)') 'POINT_DATA ', total_points
!     WRITE(fid, '(A)') 'SCALARS Pressure float 1'
!     WRITE(fid, '(A)') 'LOOKUP_TABLE default'

!     ! Write data
!     IF (TRIM(format) == 'ascii') THEN
!        DO k = 1, Nz
!           DO j = 1, Ny
!              DO i = 1, Nx
!                 value = matrix(i, j, k)
!                 WRITE(fid, '(F10.5)', ADVANCE='NO') value
!                 WRITE(fid, '(A)', ADVANCE='YES') ''  ! New line for clarity
!              END DO
!           END DO
!        END DO
!     ELSEIF (TRIM(format) == 'binary') THEN
!        CLOSE(fid)
!        OPEN(UNIT=10, FILE=filename, STATUS='APPEND', FORM='UNFORMATTED', ACCESS='STREAM')
!        DO k = 1, Nz
!           DO j = 1, Ny
!              DO i = 1, Nx
!                 WRITE(fid) matrix(i, j, k)
!              END DO
!           END DO
!        END DO
!     END IF

!     ! Close file
!     CLOSE(fid)
!   END SUBROUTINE WriteVTK

