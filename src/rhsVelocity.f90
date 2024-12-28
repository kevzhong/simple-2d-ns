! Accumulates the RHS terms in the linear system for velocity components
! That is, the advective, diffusive, and pressure-gradient terms
! Since the current treatment is fully explicit for all terms using OpenMP
! everything has been accumulated in a single (Nx,Nz) loop so that only a single instance
! of OMP thread creation/destruction is needed for brevity

subroutine rhsVelocity
    implicit none
    call build_rhsVelocity

    ! For later: e.g. RK3, terms are accumulated separately
    !call advect_u
    !call advect_w
    !call gradp

end subroutine rhsVelocity

! Calculate intermediate velocity field, ustar, to be used for fractional step
subroutine update_velocity
    use velfields
    use parameters
    use velMemory
    use ghost
    implicit none
    integer :: i, k

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k) &
    !$omp shared(Nx,Nz,dt) &
    !$omp shared(u,w,rhs_u,rhs_w)
    do k = 1,Nz
        do i = 1, Nx
            u(i,k) = u(i,k) + dt * rhs_u(i,k)
            w(i,k) = w(i,k) + dt * rhs_w(i,k)
        enddo
    enddo
    !$omp end parallel do

    !call update_ghost_periodic(u)
    !call update_ghost_periodic(w)

    call update_ghost_walls(u,w,ubot,utop,wbot,wtop)

end subroutine update_velocity

! Build all RHS vectors for velocity in Explicit Euler manner
subroutine build_rhsVelocity
    use velfields
    use parameters
    use grid
    use velMemory
    use scalarfields
    implicit none
    integer :: i, k
    real :: duudx, duwdz, dwudx, dwwdz, nudd2_u, nudd2_w, dpdx, dpdz
    real :: Tu_interp, Tw_interp

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k) &
    !$omp private(duudx,duwdz,dwudx,dwwdz,nudd2_u,nudd2_w,dpdx,dpdz,Tu_interp,Tw_interp) &
    !$omp shared(Nx,Nz,dx,dz,nu,mean_dpdx,scalarmode,beta_gx,beta_gz) &
    !$omp shared(u,w,p,temp,rhs_u,rhs_w)
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
            nudd2_u = nu * (  ( u(i-1,k) -2.0*u(i,k) + u(i+1,k)  ) / dx**2  + &
                              ( u(i,k-1) -2.0*u(i,k) + u(i,k+1)  ) / dz**2 )

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
            nudd2_w = nu * (  ( w(i-1,k) -2.0*w(i,k) + w(i+1,k)  ) / dx**2  + &
                              ( w(i,k-1) -2.0*w(i,k) + w(i,k+1)  ) / dz**2 )

            ! dp/dz | k+1/2
            dpdz = ( -p(i,k-1) + p(i,k) ) / dz


            !--------------------- BUILD RHS -----------------------!

            ! Explicit Euler
            rhs_u(i,k) =   -(dpdx + duudx + duwdz) + nudd2_u - mean_dpdx
            rhs_w(i,k) =   -(dpdz + dwudx + dwwdz) + nudd2_w

            if (scalarmode .eqv. .true.) then

                Tu_interp = 0.5 * ( temp(i-1,k) + temp(i,k) )
                Tw_interp = 0.5 * ( temp(i,k-1) + temp(i,k) )

                rhs_u(i,k) = rhs_u(i,k) - beta_gx * Tu_interp
                rhs_w(i,k) = rhs_w(i,k) - beta_gz * Tw_interp

            endif

        enddo
    enddo
    !$omp end parallel do


end subroutine build_rhsVelocity

! subroutine advect_u
!     use fields
!     use grid
!     use auxMemory, only: adv_u
!     use params
!     implicit none
!     integer :: i, k
!     real(8) :: duudx, duwdz

!     !$omp parallel do &
!     !$omp default(none) &
!     !$omp private(i,k,duudx, duwdz) &
!     !$omp shared(u,w,adv_u,dx,dz,Nx,Nz)
!     do i = 1,Nx
!         do k = 1,Nz
!             ! duu / dx
!             duudx =( (u(i+1,k)+u(i,k)) &
!                     *(u(i+1,k)+u(i,k)) &
!                     -(u(i-1,k)+u(i,k)) &
!                     *(u(i-1,k)+u(i,k)) &
!                     ) / (4.0 * dx)

!             ! duw / dz
!             duwdz = ( (w(i,k+1)+w(i-1,k+1))*(u(i,k+1)+u(i,k)) &
!                       -(w(i,k)+w(i-1,k))*(u(i,k)+u(i,k-1)) &
!                     ) / (4.0 * dz)


!             adv_u(i,k) = -(duudx + duwdz) 
!         enddo
!     enddo
!     !$omp end parallel do
! end subroutine advect_u

! subroutine advect_w
!     use fields
!     use grid
!     use auxMemory, only: adv_w
!     use params
!     implicit none
!     integer :: i, k
!     real(8) :: dwudx, dwwdz

!     !$omp parallel do &
!     !$omp default(none) &
!     !$omp private(i,k,dwudx, dwwdz) &
!     !$omp shared(u,w,adv_w,dx,dz,Nx,Nz,gadt)
!     do i = 1,Nx
!         do k = 1,Nz
!             ! dww / dz
!             dwwdz =( (w(i,k+1)+w(i,k)) &
!                     *(w(i,k+1)+w(i,k)) &
!                     -(w(i,k-1)+w(i,k)) &
!                     *(w(i,k-1)+w(i,k)) &
!                     ) / (4.0 * dz)

!             ! dwu / dx

!              dwudx = ( ( (u(i+1,k)+u(i+1,k-1)) &
!                         *(w(i+1,k)+w(i,k))) &
!                        -((u(i,k)+u(i,k-1)) &
!                         *(w(i,k)+w(i-1,k)))) / (4.0 * dx)

!             adv_w(i,k) = -(dwudx + dwwdz) 
!         enddo
!     enddo
!     !$omp end parallel do
! end subroutine advect_w


! subroutine gradp
!     use fields, only: p
!     use grid
!     use auxMemory, only: rhs_u, rhs_w
!     use params
!     implicit none
!     integer :: i, k
!     real(8) :: dp

!     !$omp parallel do &
!     !$omp default(none) &
!     !$omp private(i,k,dp) &
!     !$omp shared(p,rhs_u,rhs_w,dx,dz,Nx,Nz,aldt)
!     do i = 1,Nx
!         do k = 1,Nz

!             ! dp/dx
!             dp = ( p(i,k)-p(i-1,k) ) / dx
!             rhs_u(i,k) = - dp * aldt

!              ! dp/dz
!             dp = ( p(i,k)-p(i,k-1) ) / dz
!             rhs_w(i,k) = - dp * aldt

!         enddo
!     enddo
!     !$omp end parallel do
! end subroutine gradp