! Accumulates the RHS terms in the linear system for velocity components
! That is, the advective, diffusive, and pressure-gradient terms
! Since the current treatment is fully explicit for all terms using OpenMP
! everything has been accumulated in a single (Nx,Nz) loop so that only a single instance
! of OMP thread creation/destruction is needed for brevity

subroutine rhsScalar
    implicit none
    call build_rhsScalar
    
end subroutine rhsScalar

! Build all RHS vectors for velocity in Explicit Euler manner
subroutine build_rhsScalar
    use rk3
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
    !$omp shared(Nx,Nz,dx,dz,nu,prandtl,aldt,gamdt,zetdt) &
    !$omp shared(u,w,temp,rhs_temp,expl_c,expl_c_m1)
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

            expl_c(i,k) =  -(duTdx + dwTdz)

            rhs_temp(i,k) = gamdt*expl_c(i,k) + zetdt*expl_c_m1(i,k) & ! Explicit advection terms
                           +aldt*kap_dd2_t
        enddo
    enddo
    !$omp end parallel do


end subroutine build_rhsScalar