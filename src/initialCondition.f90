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
    real :: Ub, Av, Aw

    !Ub = -mean_dpdx * (0.5*Lz)**2 / (3.0 * nu) !bulk Laminar velocity
    !Ub = -mean_dpdx * (0.5*Lz)**2 / (3.0 * 1.0) !bulk Laminar velocity
    Ub = 15.0

    Aw = 0.1 * Ub
    Av = -0.25 * Aw * Ly



    ! Taylor--Green vortices

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,j,k) &
    !$omp shared(u,v,w,xc,xm,zc,zm,Nx,Ny,Nz,Lx,Ly,Lz,PI) &
    !$omp shared(Ub,Av,Aw,nu,mean_dpdx,yc,ym)
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

                ! Pirozzoli (2019): Laminar Poiseuille w/ longitudinal rollers
                !u(i,j,k) = 1.0/nu*mean_dpdx*zm(k) * (0.5 * zm(k) - 0.5*Lz )
                !u(i,j,k) = 1.0*mean_dpdx*zm(k) * (0.5 * zm(k) - 0.5*Lz )
                !u(i,j,k) = -Ub*3.0*zm(k) * (0.5 * zm(k) / Lz - 1.0 )
                !v(i,j,k) = Av * cos( 0.5 * PI * zm(k) ) * sin(2.0 * PI * yc(j) / Ly )
                !w(i,j,k) = Aw * sin( 0.5 * PI * zc(k) ) * cos(2.0 * PI * ym(j) / Ly )


                !u(i,k) =  sin( 2.0*PI * xc(i) / Lx ) * cos( 2.0*PI * zm(k) / Lz )
                !w(i,k) = -cos( 2.0*PI * xm(i) / Lx ) * sin( 2.0*PI * zc(k) / Lz )

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


subroutine readRestart
    use velfields
    use grid
    use parameters
    use ghost
    use scalarfields
    implicit none

    call read3DField(u(1:Nx,1:Ny,1:Nz),Nx,Ny,Nz,'u',nt0)
    call read3DField(v(1:Nx,1:Ny,1:Nz),Nx,Ny,Nz,'v',nt0)
    call read3DField(w(1:Nx,1:Ny,1:Nz),Nx,Ny,Nz,'w',nt0)
    call update_ghost_wallsU(u,bctype_ubot,bctype_utop,bcval_ubot,bcval_utop)
    call update_ghost_wallsU(v,bctype_vbot,bctype_vtop,bcval_vbot,bcval_vtop)
    call update_ghost_wallsW(w,bctype_wbot,bctype_wtop,bcval_wbot,bcval_wtop)

    if (scalarmode .eqv. .true. ) then
        call read3DField(temp(1:Nx,1:Ny,1:Nz),Nx,Ny,Nz,'C',nt0)
        call update_ghost_wallTemp(temp,bctype_Tbot,bctype_Ttop,bcval_Tbot,bcval_Ttop)
    endif

end subroutine readRestart