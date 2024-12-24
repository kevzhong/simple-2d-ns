! Setup the compuational grid
! Cell centres: xm, cell edges: xc

subroutine generate_grid
    use parameters, only: Lx, Lz, Nx, Nz
    use grid
    use ghost, only: halosize
    implicit none
    integer :: i

    ! Grid spacing
    dx = Lx / Nx
    dz = Lz / Nz


    ! Allocate memory for spatial locations

    allocate( xc( -halosize+1:Nx+halosize  )  )
    allocate( xm( -halosize+1:Nx+halosize  )  )

    allocate( zc( -halosize+1:Nz+halosize  )  )
    allocate( zm( -halosize+1:Nz+halosize  )  )

    do i = -halosize+1,Nx+halosize
        xc(i) = dble(i - 1) * dx
        xm(i) = ( dble(i-1) + 0.5 ) * dx
    enddo

    do i = -halosize+1,Nz+halosize
        zc(i) = dble(i - 1) * dz
        zm(i) = ( dble(i-1) + 0.5 ) * dz
    enddo

end subroutine generate_grid


subroutine dealloc_grid
    use grid

    if(allocated(xc)) deallocate(xc)
    if(allocated(xm)) deallocate(xm)

    if(allocated(zc)) deallocate(zc)
    if(allocated(zm)) deallocate(zm)

end subroutine dealloc_grid
