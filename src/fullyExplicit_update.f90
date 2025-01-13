subroutine fullyExplicit_update
    use parameters
    implicit none

    call update_velocity_fullyExplicit
    if (scalarmode .eqv. .true.) call update_scalar_fullyExplicit ! New timestep value for temperature

end subroutine fullyExplicit_update

! Calculate intermediate velocity field, ustar, to be used for fractional step
! Explicit Euler manner
subroutine update_velocity_fullyExplicit
    use velfields
    use parameters
    use velMemory
    use ghost
    implicit none
    integer :: i, j, k

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,j,k) &
    !$omp shared(Nx,Ny,Nz) &
    !$omp shared(u,v,w,rhs_u,rhs_v,rhs_w)
    do k = 1,Nz
        do j = 1,Ny
            do i = 1, Nx
                u(i,j,k) = u(i,j,k) +  rhs_u(i,j,k)
                v(i,j,k) = v(i,j,k) +  rhs_v(i,j,k)
                w(i,j,k) = w(i,j,k) +  rhs_w(i,j,k)
            enddo
        enddo
    enddo
    !$omp end parallel do

    call update_ghost_wallsU(u,bctype_ubot,bctype_utop,bcval_ubot,bcval_utop)
    call update_ghost_wallsU(v,bctype_vbot,bctype_vtop,bcval_vbot,bcval_vtop)
    call update_ghost_wallsW(w,bctype_wbot,bctype_wtop,bcval_wbot,bcval_wtop)


end subroutine update_velocity_fullyExplicit


! Calculate intermediate velocity field, ustar, to be used for fractional step
subroutine update_scalar_fullyExplicit
    use scalarfields
    use scalarMemory
    use parameters
    use ghost
    implicit none
    integer :: i, j, k

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,j,k) &
    !$omp shared(Nx,Ny,Nz) &
    !$omp shared(temp,rhs_temp)
    do k = 1,Nz
        do j = 1,Ny
            do i = 1, Nx
                temp(i,j,k) = temp(i,j,k) + rhs_temp(i,j,k)
            enddo
        enddo
    enddo
    !$omp end parallel do


    call update_ghost_wallTemp(temp,bctype_Tbot,bctype_Ttop,bcval_Tbot,bcval_Ttop)

end subroutine update_scalar_fullyExplicit

