subroutine initialCondition
    use velfields
    use grid
    use parameters
    use ghost
    implicit none
    integer :: i, k
    real(8) :: PI = 2.d0*dasin(1.d0) 



    ! Taylor--Green vortices

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k) &
    !$omp shared(u,w,xc,xm,zc,zm,Nx,Nz,Lx,Lz,PI)
    do i = 1,Nx
        do k = 1,Nz
            u(i,k) =  sin( 2.0*PI * xc(i) / Lx ) * cos( 2.0*PI * zm(k) / Lz )
            w(i,k) = -cos( 2.0*PI * xm(i) / Lx ) * sin( 2.0*PI * zc(k) / Lz )
        enddo
    enddo
    !$omp end parallel do

    call update_ghost(u)
    call update_ghost(w)

end subroutine initialCondition