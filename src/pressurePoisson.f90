subroutine pressurePoisson
    use parameters
    use velfields, only: pseudo_p
    use fftMemory
    use ghost
    use bctypes
    implicit none

    call build_rhsPoisson
    !call solve_pressurePoisson    
    call solve_helmholtz(pseudo_p(1:Nx,1:Nz), rhs_poisson, 0.0, -1.0, PRESSUREBC, PRESSUREBC )
    call update_ghost_pressure(pseudo_p)
    call projectionUpdate


end subroutine pressurePoisson

subroutine build_rhsPoisson
    use velfields, only: u,w
    use grid
    use rk3
    use parameters
    use fftMemory, only: rhs_poisson
    implicit none
    integer :: i, k
    real :: divu


    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k,divu) &
    !$omp shared(u,w,rhs_poisson,dx,dz,Nx,Nz,aldt)
    do k = 1,Nz
        do i = 1,Nx
            divu = ( u(i+1,k) - u(i,k) ) / dx + &
                   ( w(i,k+1) - w(i,k) ) / dz 

            rhs_poisson(i,k) = divu / aldt
        enddo
    enddo
    !$omp end parallel do
end subroutine build_rhsPoisson

subroutine projectionUpdate
    use velfields
    use rk3
    use grid
    use ghost
    use parameters
    implicit none
    integer :: i, k
    real :: half_nualdt

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k) &
    !$omp shared(u,w,p,pseudo_p,dx,dz,Nx,Nz,aldt)
    do k = 1,Nz
        do i = 1,Nx
            u(i,k) = u(i,k) - aldt * ( pseudo_p(i,k) - pseudo_p(i-1,k) ) / dx
            w(i,k) = w(i,k) - aldt * ( pseudo_p(i,k) - pseudo_p(i,k-1) ) / dz

            p(i,k) = p(i,k) + pseudo_p(i,k)

        enddo
    enddo
    !$omp end parallel do


    ! Additional Laplacian contribution from implicit treatment
    if (implicitmode .eqv. .true.) then
        half_nualdt = 0.5 * nu * aldt

        !$omp parallel do &
        !$omp default(none) &
        !$omp private(i,k) &
        !$omp shared(p,pseudo_p,dx,dz,Nx,Nz,dt,half_nualdt)
        do k = 1,Nz
            do i = 1,Nx
                p(i,k) = p(i,k) - half_nualdt * ( ( pseudo_p(i-1,k) - 2.0 * pseudo_p(i,k) + pseudo_p(i+1,k) ) / dx**2 &
                                                + ( pseudo_p(i,k-1) - 2.0 * pseudo_p(i,k) + pseudo_p(i,k+1) ) / dz**2 )
            enddo
        enddo
        !$omp end parallel do
    endif

    call update_ghost_pressure(p)

end subroutine projectionUpdate