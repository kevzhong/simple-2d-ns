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
    integer :: i, j,k
    real :: duTdx, dvTdy,  dwTdz,  kap_dd2_t, kap_dd2xy_t
    real :: dzph, dzmh


    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,j,k) &
    !$omp private(dzph,dzmh,duTdx,dvTdy,dwTdz,kap_dd2_t,kap_dd2xy_t) &
    !$omp shared(Nx,Ny,Nz,dx,dy,dz,nu,prandtl,aldt,gamdt,zetdt) &
    !$omp shared(u,v,w,temp,rhs_temp,implicitXYmode,expl_c,expl_c_m1)
    do k = 1,Nz
        dzmh = 0.5*( dz(k-1) + dz(k  ) )
        dzph = 0.5*( dz(k  ) + dz(k+1) )
        do j = 1,Ny
            do i = 1,Nx


            ! d uT  |                1     [                        ]
            ! ----- |            =  ----   |   uT       -   uT      | 
            !  dx   |i,j,k           dx    [    i+1/2        i-1/2  ]


            duTdx=( u(i+1,j,k)*(temp(i+1,j,k)+temp(i  ,j,k)) & 
                   -u(i,  j,k)*(temp(i,  j,k)+temp(i-1,j,k)) ) / (2.0 * dx)

            ! d vT  |                1     [                        ]
            ! ----- |            =  ----   |   vT       -   vT      | 
            !  dy   |i,j,k           dy    [    j+1/2        j-1/2  ]


            dvTdy=( v(i,j+1,k)*(temp(i,j+1,k)+temp(i,j  ,k)) & 
                   -v(i,j  ,k)*(temp(i,j  ,k)+temp(i,j-1,k)) ) / (2.0 * dy)

            ! d wT  |                1     [                        ]
            ! ----- |            =  ----   |   wT    -     wT       | 
            !  dz   |i,j,k           dz    [    k+1/2        k-1/2  ] 

            dwTdz=( w(i,j,k+1)*(temp(i,j,k+1)+temp(i,j,k)) & 
                   -w(i,j,k)*(temp(i,j,k)+temp(i,j,k-1)) ) / (2.0 * dz(k) )


            ! d2T / dxj2
            kap_dd2xy_t = nu/prandtl * ( ( temp(i-1,j,k) -2.0*temp(i,j,k) + temp(i+1,j,k)  ) / dx**2  + &
                                         ( temp(i,j-1,k) -2.0*temp(i,j,k) + temp(i,j+1,k)  ) / dy**2 )

            kap_dd2_t = nu/prandtl * (   ( temp(i,j,k+1) - temp(i,j,k  ) ) / dzph  - &
                                         ( temp(i,j,k  ) - temp(i,j,k-1) ) / dzmh   ) / dz(k)

            !--------------------- BUILD RHS -----------------------!

            expl_c(i,j,k) =  -(duTdx + dvTdy + dwTdz)

            if (implicitXYmode .eqv. .true.) then
                kap_dd2_t = kap_dd2_t + kap_dd2xy_t
            else
                expl_c(i,j,k) = expl_c(i,j,k) + kap_dd2xy_t
            endif

            rhs_temp(i,j,k) = gamdt*expl_c(i,j,k) + zetdt*expl_c_m1(i,j,k) & ! Explicit advection terms
                             +aldt*kap_dd2_t

            enddo
        enddo
    enddo
    !$omp end parallel do


end subroutine build_rhsScalar