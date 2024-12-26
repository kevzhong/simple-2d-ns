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
    real, intent(in) :: array(Nx,Ny)    
    real, intent(in) :: dx, dy  
    integer, intent(in) :: it

    ! Local variables
    integer :: i, j, num_points
    integer :: file_unit
    character(len=50) :: header_line
    character(len=10) :: nx_str, ny_str


    write(filename, '("outputdir/", A, "_it", I0, ".vtk")') trim(fpref), it


    ! Calculate the total number of points
    num_points = Nx * Ny

    ! Open the file for binary write
    open(newunit=file_unit, file=filename, status='replace', form='unformatted', action='write')

    ! Write VTK header as unformatted binary strings
    header_line = '# vtk DataFile Version 3.0'
    write(file_unit) header_line

    header_line = '2D structured grid data'
    write(file_unit) header_line

    header_line = 'BINARY'
    write(file_unit) header_line

    header_line = 'DATASET STRUCTURED_POINTS'
    write(file_unit) header_line

    ! Convert integers to strings for DIMENSIONS
    write(nx_str, '(I0)') Nx
    write(ny_str, '(I0)') Ny
    header_line = 'DIMENSIONS ' // trim(nx_str) // ' ' // trim(ny_str) // ' 1'
    write(file_unit) header_line

    ! Write ORIGIN line
    header_line = 'ORIGIN ' // trim(adjustl("0.0")) // ' ' // trim(adjustl("0.0")) // ' 0.0'
    write(file_unit) header_line

    ! Write SPACING line
    write(header_line, '(A, E16.8, A, E16.8, A, E16.8)') 'SPACING ', dx, ' ', dy, ' 1.0'
    write(file_unit) header_line

    ! Write SCALARS and LOOKUP_TABLE
    header_line = 'SCALARS scalars double'
    write(file_unit) header_line

    header_line = 'LOOKUP_TABLE default'
    write(file_unit) header_line

    write(*,*) "Writing array now!"
    ! Write the array in binary format (row-major order)
    do j = 1, Ny
        do i = 1, Nx
            write(file_unit) array(i, j)
        end do
    end do

    ! Close the file
    close(file_unit)

end subroutine writeVTK