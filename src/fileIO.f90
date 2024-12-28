subroutine write2DField(array,Nx,Ny,dx,dy,fpref,it)
    ! Write array(Nx,Ny) field as binary data
    implicit none
    integer :: Nx, Ny
    real(8), intent(in) :: dx, dy   ! Grid spacing: (dx, dy)
    character(len=50) :: filename 
    real, dimension(Nx,Ny) :: array
    integer :: it
    character(len=*) :: fpref

    write(filename, '("outputdir/",A, "_it", I0, ".dat")') trim(fpref), it
    open(unit=10, file=filename, status='replace', access='stream', action='write')
    write(10) array
    close(10)

    call write_xdmf(fpref, it, Nx, Ny, dx, dy) ! Writing header for viewing in Paraview
end subroutine write2DField

subroutine write_xdmf(fpref, it,nx, ny, dx,dy)
    implicit none
    character(len=50) :: filename  ! xdmf file name
    character(len=50) :: bfilename ! binary filename
    character(len=*) :: fpref
    integer :: it
    integer, intent(in)   :: nx, ny          
    real(8), intent(in) :: dx, dy   

    character(len=32) :: nx_str, ny_str
    integer :: unit       

    write(nx_str, '(I0)') Nx
    write(ny_str, '(I0)') Ny



    write(filename, '("outputdir/",A, "_it", I0, ".xmf")') trim(fpref), it
    write(bfilename, '(A, "_it", I0, ".dat")') trim(fpref), it


    ! Open the file
    open(newunit=unit, file=filename, status='replace', action='write')
    ! Write the XDMF header
    write(unit, '(A)') '<?xml version="1.0" ?>'
    write(unit, '(A)') '<Xdmf Version="3.0">'
    write(unit, '(A)') '  <Domain>'
    write(unit, '(A)') '    <Grid Name="StructuredGrid" GridType="Uniform">'
    write(unit, '(A)') '      <Topology TopologyType="2DCoRectMesh" Dimensions="' // &
                       trim(adjustl(ny_str)) // ' ' // &
                       trim(adjustl(nx_str)) // '"/>'
    write(unit, '(A)') '      <Geometry GeometryType="ORIGIN_DXDY">'
    write(unit, '(A)') '      <DataItem Format="XML" Dimensions="2">0.0 0.0</DataItem>'

    write(unit, '(A,F12.6, A, F12.6, A)') '        <DataItem Format="XML" Dimensions="2">', dx, ' ', dy, '</DataItem>'

    write(unit, '(A)') '      </Geometry>'
    write(unit, '(A)') '      <Attribute Name="' // fpref // '" AttributeType="Scalar" Center="Node">'
    write(unit, '(A)') '        <DataItem Format="Binary" Endian="Little" DataType="Float" Precision="8" Dimensions="' // &
                                trim(adjustl(ny_str)) // ' ' // &
                                trim(adjustl(nx_str)) // '">' // trim(adjustl(bfilename)) // '</DataItem>'

    write(unit, '(A)') '      </Attribute>'
    write(unit, '(A)') '    </Grid>'
    write(unit, '(A)') '  </Domain>'
    write(unit, '(A)') '</Xdmf>'

    ! Close the file
    close(unit)

end subroutine write_xdmf



! subroutine writeVTK(array, Nx, Ny, dx, dy,fpref,it)
!     ! Write array(Nx,Ny) field as binary .vtk file for convenient viewing in Paraview
!     implicit none

!     character(len=50) :: filename  
!     character(len=*), intent(in) :: fpref  
!     integer, intent(in) :: Nx, Ny
!     real(8), intent(in) :: array(Nx,Ny)   
!     real, intent(in) :: dx, dy  
!     integer, intent(in) :: it

!     ! Local variables
!     integer :: i, j, num_points
!     integer :: file_unit


!     write(filename, '("outputdir/", A, "_it", I0, ".vtk")') trim(fpref), it


!     ! Calculate the total number of points
!     num_points = Nx * Ny

!     open(newunit=file_unit, file=filename, status='replace', form='formatted', action='write')

!     ! Write VTK ASCII header
!     write(file_unit, '(A)') '# vtk DataFile Version 3.0'
!     write(file_unit, '(A)') '2D structured grid data'
!     !write(file_unit, '(A)') 'ASCII'
!     write(file_unit, '(A)') 'BINARY'
!     write(file_unit, '(A)') 'DATASET STRUCTURED_POINTS'
!     write(file_unit, '(A, I0, A, I0, A)') 'DIMENSIONS ', Nx, ' ', Ny, ' 1'
!     write(file_unit, '(A, F10.5, A, F10.5, A)') 'SPACING ', dx, ' ', dy, ' 1.0'
!     write(file_unit, '(A)') 'ORIGIN 0.0 0.0 0.0'
!     write(file_unit, '(A, I0)') 'POINT_DATA ', num_points
!     write(file_unit, '(A)') 'SCALARS field double'
!     write(file_unit, '(A)') 'LOOKUP_TABLE default'

!     ! BINARY WRITING
!     close(file_unit)
!     ! Open the file in unformatted (binary) mode to append the data
!     open(newunit=file_unit, file=filename, position='append', form='unformatted', access='stream')
!     do j = 1, Ny
!         do i = 1, Nx
!             write(file_unit) array(i, j)
!         end do
!     end do


!     ! ! ASCII WRITING
!     ! DO j = 1, Ny
!     !     DO i = 1, Nx
!     !        VAL = array(i, j)
!     !        WRITE(file_unit, '(F10.5)', ADVANCE='NO') VAL
!     !        WRITE(file_unit, '(A)', ADVANCE='YES') ''  ! New line for clarity
!     !     END DO
!     !  END DO


!     ! Close the file
!     close(file_unit)

! end subroutine writeVTK



