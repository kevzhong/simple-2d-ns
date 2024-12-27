subroutine initialCondition
    use velfields
    use grid
    use parameters
    use ghost
    use scalarfields
    implicit none
    integer :: i, k
    real(8) :: PI = 2.d0*dasin(1.d0) 



    ! Taylor--Green vortices

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k) &
    !$omp shared(u,w,xc,xm,zc,zm,Nx,Nz,Lx,Lz,PI)
    do k = 1,Nz
        do i = 1,Nx
            !u(i,k) =  sin( 2.0*PI * xc(i) / Lx ) * cos( 2.0*PI * zm(k) / Lz )
            !w(i,k) = -cos( 2.0*PI * xm(i) / Lx ) * sin( 2.0*PI * zc(k) / Lz )

            !u(i,k) =  sin(  xc(i) ) * cos(  zm(k)  )
            !w(i,k) = -cos(  xm(i) ) * sin(  zc(k)  )

            u(i,k) = 0.0
            w(i,k) = 0.0
        enddo
    enddo
    !$omp end parallel do
    call update_ghost_walls(u,w,ubot,utop,wbot,wtop)

    if (scalarmode .eqv. .true.) then
        !$omp parallel do &
        !$omp default(none) &
        !$omp private(i,k) &
        !$omp shared(temp,Nx,Nz)
        do k = 1,Nz
            do i = 1,Nx
                temp(i,k) = 0.0
        enddo
    enddo
    !$omp end parallel do

    call update_ghost_wallTemp(temp,Tbot,Ttop)
    
    endif

end subroutine initialCondition