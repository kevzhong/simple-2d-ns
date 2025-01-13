subroutine initialCondition
    use velfields
    use grid
    use parameters
    use ghost
    use scalarfields
    implicit none
    integer :: i, j,k
    real(8) :: PI = 2.d0*dasin(1.d0) 
    real :: eps



    ! Taylor--Green vortices

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,j,k) &
    !$omp shared(u,v,w,xc,xm,zc,zm,Nx,Ny,Nz,Lx,Lz,PI)
    do k = 1,Nz
        do j = 1,Ny
            do i = 1,Nx
                !u(i,k) =  sin( 2.0*PI * xc(i) / Lx ) * cos( 2.0*PI * zm(k) / Lz )
                !w(i,k) = -cos( 2.0*PI * xm(i) / Lx ) * sin( 2.0*PI * zc(k) / Lz )

                !u(i,k) =  sin(  xc(i) ) * cos(  zm(k)  )
                !w(i,k) = -cos(  xm(i) ) * sin(  zc(k)  )

                u(i,j,k) = 0.0
                v(i,j,k) = 0.0
                w(i,j,k) = 0.0
            enddo
        enddo
    enddo
    !$omp end parallel do
    call update_ghost_wallsU(u,bctype_ubot,bctype_utop,bcval_ubot,bcval_utop)
    call update_ghost_wallsU(v,bctype_vbot,bctype_vtop,bcval_vbot,bcval_vtop)
    call update_ghost_wallsW(w,bctype_wbot,bctype_wtop,bcval_wbot,bcval_wtop)
    !call update_ghost_pressure(p)

    if (scalarmode .eqv. .true.) then
        call random_seed()
        !$omp parallel do &
        !$omp default(none) &
        !$omp private(eps,i,j,k) &
        !$omp shared(temp,zm,bcval_Tbot,Lz,bcval_Ttop,Nx,Ny,Nz)
        do k = 1,Nz
            do j = 1,Ny
                do i = 1,Nx
                    temp(i,j,k) = (bcval_Ttop - bcval_Tbot) / Lz * zm(k) + bcval_Tbot 

                    ! Super-impose perturbation
                    call random_number(eps)
                    temp(i,j,k) = temp(i,j,k) + 0.2 * (eps - 0.5)
                enddo
            enddo
        enddo
        !$omp end parallel do

    call update_ghost_wallTemp(temp,bctype_Tbot,bctype_Ttop,bcval_Tbot,bcval_Ttop)
    
    endif

end subroutine initialCondition