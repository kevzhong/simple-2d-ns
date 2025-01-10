! Accumulates the RHS terms in the linear system for velocity components
! That is, the advective, diffusive, and pressure-gradient terms
! Since the current treatment is fully explicit for all terms using OpenMP
! everything has been accumulated in a single (Nx,Nz) loop so that only a single instance
! of OMP thread creation/destruction is needed for brevity

subroutine rhsVelocity
    implicit none
    call build_rhsVelocity
end subroutine rhsVelocity

! Build all RHS vectors for velocity
subroutine build_rhsVelocity
    use velfields
    use rk3
    use parameters
    use grid
    use velMemory
    use scalarfields
    implicit none
    integer :: i, k
    real :: duudx, duwdz, dwudx, dwwdz, nudd2_u, nudd2_w, dpdx, dpdz
    real :: Tu_interp, Tw_interp

    ! Staggered grid arrangement, the stagger is by -1/2, index convention:

    !                     v(i,j+1,k)
    !                      ^
    !                      |
    !                      |
    !    __________________|___________________
    !   |                                      |
    !   |                                      |
    !   |                                      |
    !   |                                      |
    !   |  u(i,j,k)                            |
    ! ----->                               ------>  u(i+1,j,k)
    !   |                  O                   |
    !   |              temp(i,j,k)             |
    !   |              or p(i,j,k)             |
    !   |                                      |
    !   |                                      |
    !   |                  ^                   |
    !   |__________________|___________________
    !                      |
    !                      |
    !                     v(i,j,k)

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k) &
    !$omp private(duudx,duwdz,dwudx,dwwdz,nudd2_u,nudd2_w,dpdx,dpdz,Tu_interp,Tw_interp) &
    !$omp shared(Nx,Nz,dx,dz,nu,mean_dpdx,scalarmode,implicitXmode,beta_gx,beta_gz,gamdt,zetdt,aldt) &
    !$omp shared(u,w,p,temp,rhs_u,rhs_w,expl_u,expl_w,expl_u_m1,expl_w_m1)
    do k = 1,Nz
        do i = 1,Nx

            !--------------------- U VELOCITY -----------------------!

            ! d uu  |               1     [                  ]
            ! ----- |           =  ----   |   uu    -   uu   | 
            !  dx   |i+1/2, k       dx    [    i+1        i  ]

            duudx =( (u(i+1,k)+u(i,k)) &
                     *(u(i+1,k)+u(i,k)) &
                     -(u(i-1,k)+u(i,k)) &
                     *(u(i-1,k)+u(i,k)) &
                     ) / (4.0 * dx)

            ! d uw  |               1     [                  ]
            ! ----- |           =  ----   |   uw    -   uw   | 
            !  dz   |i+1/2, k       dz    [    i+1        i  ] 

            duwdz = ( (w(i,k+1)+w(i-1,k+1))*(u(i,k+1)+u(i,k)) &
                    -(w(i,k)+w(i-1,k))*(u(i,k)+u(i,k-1)) &
                     ) / (4.0 * dz)      


            ! d2u / dxj2
            nudd2_u = nu * ( u(i,k-1) -2.0*u(i,k) + u(i,k+1)  ) / dz**2 
                              

            ! dp/dx | i+1/2
             dpdx = ( -p(i-1,k) + p(i,k) ) / dx


            !--------------------- W VELOCITY -----------------------!
                     
            ! d ww  |               1     [                  ]
            ! ----- |           =  ----   |   ww    -   ww   | 
            !  dz   |i, k+1/2       dz    [    k+1        k  ] 

             dwwdz =( (w(i,k+1)+w(i,k)) &
                     *(w(i,k+1)+w(i,k)) &
                     -(w(i,k-1)+w(i,k)) &
                     *(w(i,k-1)+w(i,k)) &
                     ) / (4.0 * dz)

            ! d wu  |               1     [                  ]
            ! ----- |           =  ----   |   wu    -   wu   | 
            !  dx   |i, k+1/2       dx    [    k+1        k  ]

            dwudx = ( ( (u(i+1,k)+u(i+1,k-1)) &
                    *(w(i+1,k)+w(i,k))) &
                    -((u(i,k)+u(i,k-1)) &
                    *(w(i,k)+w(i-1,k)))) / (4.0 * dx)

            ! d2w / dxj2
            nudd2_w = nu  * ( w(i,k-1) -2.0*w(i,k) + w(i,k+1)  ) / dz**2 

            ! dp/dz | k+1/2
            dpdz = ( -p(i,k-1) + p(i,k) ) / dz


            !--------------------- BUILD RHS -----------------------!

            expl_u(i,k) = -( duudx + duwdz)
            expl_w(i,k) = -( dwudx + dwwdz)

            if (implicitXmode .eqv. .true.) then
                nudd2_u = nudd2_u + nu *   ( u(i-1,k) -2.0*u(i,k) + u(i+1,k)  ) / dx**2 
                nudd2_w = nudd2_w + nu *   ( w(i-1,k) -2.0*w(i,k) + w(i+1,k)  ) / dx**2 
            else
                expl_u(i,k) = expl_u(i,k) + nu * ( u(i-1,k) -2.0*u(i,k) + u(i+1,k)  ) / dx**2
                expl_w(i,k) = expl_w(i,k) + nu * ( w(i-1,k) -2.0*w(i,k) + w(i+1,k)  ) / dx**2 
            endif

            if (scalarmode .eqv. .true.) then

                Tu_interp = 0.5 * ( temp(i-1,k) + temp(i,k) )
                Tw_interp = 0.5 * ( temp(i,k-1) + temp(i,k) )

                expl_u(i,k) = expl_u(i,k) + ( beta_gx * Tu_interp )
                expl_w(i,k) = expl_w(i,k) + ( beta_gz * Tw_interp )

            endif

            rhs_u(i,k) = gamdt*expl_u(i,k) + zetdt*expl_u_m1(i,k)  & ! Fully-explicit terms (advection + buoyancy)
                        -aldt * dpdx  & ! Pressure gradient
                        +aldt * nudd2_u & ! Diffusive term: semi-implicit part
                        -aldt * mean_dpdx ! Mean pressure gradient body forcing

            rhs_w(i,k) = gamdt*expl_w(i,k) + zetdt*expl_w_m1(i,k)  & ! Fully-explicit terms (advection + buoyancy)
                        -aldt * dpdz  & ! Pressure gradient
                        +aldt * nudd2_w ! Diffusive term: semi-implicit part

        enddo
    enddo
    !$omp end parallel do


end subroutine build_rhsVelocity