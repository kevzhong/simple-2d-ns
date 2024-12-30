! Accumulates the RHS terms in the linear system for velocity components
! That is, the advective, diffusive, and pressure-gradient terms
! Since the current treatment is fully explicit for all terms using OpenMP
! everything has been accumulated in a single (Nx,Nz) loop so that only a single instance
! of OMP thread creation/destruction is needed for brevity

subroutine rhsScalar
    implicit none
    call build_rhsScalar
    
end subroutine rhsScalar

! Calculate intermediate velocity field, ustar, to be used for fractional step
subroutine update_scalar_fullyExplicit
    use scalarfields
    use scalarMemory
    use parameters
    use ghost
    implicit none
    integer :: i, k

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k) &
    !$omp shared(Nx,Nz,dt) &
    !$omp shared(temp,rhs_temp)
    do k = 1,Nz
        do i = 1, Nx
            temp(i,k) = temp(i,k) + dt * rhs_temp(i,k)
        enddo
    enddo
    !$omp end parallel do


    call update_ghost_wallTemp(temp,Tbot,Ttop)

end subroutine update_scalar_fullyExplicit

! Build all RHS vectors for velocity in Explicit Euler manner
subroutine build_rhsScalar
    use scalarfields
    use velfields
    use parameters
    use grid
    use scalarMemory
    implicit none
    integer :: i, k
    real :: duTdx,  dwTdz,  kap_dd2_t


    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k) &
    !$omp private(duTdx,dwTdz,kap_dd2_t) &
    !$omp shared(Nx,Nz,dx,dz,nu,prandtl) &
    !$omp shared(u,w,temp,rhs_temp)
    do k = 1,Nz
        do i = 1,Nx


            ! d uT  |               1     [                        ]
            ! ----- |           =  ----   |   uT       -   uT      | 
            !  dx   |i, k           dx    [    i+1/2        i-1/2  ]


            duTdx=( u(i+1,k)*(temp(i+1,k)+temp(i,k)) & 
                -u(i,k)*(temp(i,k)+temp(i-1,k)) ) / (2.0 * dx)

            ! d wT  |               1     [                        ]
            ! ----- |           =  ----   |   wT    -     wT       | 
            !  dz   |i, k           dz    [    k+1/2        k-1/2  ] 

            dwTdz=( w(i,k+1)*(temp(i,k+1)+temp(i,k)) & 
                -w(i,k)*(temp(i,k)+temp(i,k-1)) ) / (2.0 * dz)


            ! d2T / dxj2
            kap_dd2_t = nu/prandtl * (  ( temp(i-1,k) -2.0*temp(i,k) + temp(i+1,k)  ) / dx**2  + &
                                        ( temp(i,k-1) -2.0*temp(i,k) + temp(i,k+1)  ) / dz**2 )

            !--------------------- BUILD RHS -----------------------!

            ! Explicit Euler
            rhs_temp(i,k) =   -(duTdx + dwTdz) + kap_dd2_t
        enddo
    enddo
    !$omp end parallel do


end subroutine build_rhsScalar