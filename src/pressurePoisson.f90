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
    use parameters
    use fftMemory, only: rhs_poisson
    implicit none
    integer :: i, k
    real :: divu


    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k,divu) &
    !$omp shared(u,w,rhs_poisson,dx,dz,Nx,Nz,dt)
    do k = 1,Nz
        do i = 1,Nx
            divu = ( u(i+1,k) - u(i,k) ) / dx + &
                   ( w(i,k+1) - w(i,k) ) / dz 

            rhs_poisson(i,k) = divu / dt
        enddo
    enddo
    !$omp end parallel do
end subroutine build_rhsPoisson

! subroutine solve_pressurePoisson
!     use fftw3
!     use velfields, only: pseudo_p
!     use grid
!     use ghost
!     use fftMemory
!     use parameters
!     implicit none
!     integer :: i, k
!     real(8), allocatable :: am(:), ac(:), ap(:) ! tridiagonal coefficients
!     real(8), allocatable :: rbuffer(:), cbuffer(:) ! solution vector to rhs
!     real(8) :: a, b, c

!     allocate( am(1:Nz) ) ; allocate( ac(1:Nz) ) ; allocate( ap(1:Nz) ) ; allocate( rbuffer(1:Nz) ) ; allocate( cbuffer(1:Nz) )

!     ! Tri-diagonal inversion for each kx wavenumber

!     !$omp parallel do &
!     !$omp default(none) &
!     !$omp private(k) &
!     !$omp shared(fftw_plan_fwd,rhs_poisson,rhs_hat,Nz)
!     do k = 1,Nz
!         ! FFT in x direction for each z-location k
!         call fftw_execute_dft_r2c(fftw_plan_fwd, rhs_poisson(:,k), rhs_hat(:,k))
!     enddo
!     !$omp end parallel do

!     ! Build and solve tri-diagonal system for each kx wavenumber

!     a = 1.0 / dz**2
!     c = 1.0 / dz**2

!     !$omp parallel do &
!     !$omp default(none) &
!     !$omp private(i,b,k,am,ac,ap,rbuffer,cbuffer) &
!     !$omp shared(a,c,lmb_x_on_dx2,rhs_hat,soln_hh_hat,dz,Nx,Nz)
!     do i = 1,Nx/2+1

!         b = lmb_x_on_dx2(i) - 2.0 / dz**2

!         do k = 1,Nz

!             rhs_hat(i,k) = rhs_hat(i,k) / Nx

!             am(k) = a
!             ac(k) = b
!             ap(k) = c

!             rbuffer(k) = real(rhs_hat(i,k))
!             cbuffer(k) = aimag(rhs_hat(i,k))
!         enddo

!         ! Arbitrary Dirichlet for 0 mode, over-ride
!         if (i .eq. 1) then
!             ap(1) = 0.0
!             ac(1) = 1.0
!             rbuffer(1) = 0.0
!             cbuffer(1) = 0.0
!         endif


!         call tridiag(am,ac,ap,rbuffer,Nz)
!         call tridiag(am,ac,ap,cbuffer,Nz)

!         do k = 1,Nz
!             soln_hh_hat(i,k) = CMPLX( rbuffer(k), cbuffer(k) )
!         enddo
!     enddo
!     !$omp end parallel do


!     !$omp parallel do &
!     !$omp default(none) &
!     !$omp private(i,k) &
!     !$omp shared(fftw_plan_bwd,soln_hh_hat,soln_hh,Nz)
!     do k = 1,Nz
!         call fftw_execute_dft_c2r(fftw_plan_bwd, soln_hh_hat(:,k), soln_hh(:,k))
!     enddo
!     !$omp end parallel do


!     !$omp parallel do &
!     !$omp default(none) &
!     !$omp private(i,k) &
!     !$omp shared(soln_hh,pseudo_p,Nx,Nz)
!     do k = 1,Nz
!         do i = 1,Nx
!         pseudo_p(i,k) = soln_hh(i,k)
!         enddo
!     enddo
!     !$omp end parallel do

!     call update_ghost_pressure(pseudo_p)

!     deallocate(am) ; deallocate(ac) ; deallocate(ap) ; deallocate(rbuffer) ; deallocate(cbuffer)
    
! end subroutine solve_pressurePoisson



subroutine projectionUpdate
    use velfields
    use grid
    use ghost
    use parameters
    implicit none
    integer :: i, k
    real :: half_nudt

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k) &
    !$omp shared(u,w,p,pseudo_p,dx,dz,Nx,Nz,dt)
    do k = 1,Nz
        do i = 1,Nx
            u(i,k) = u(i,k) - dt * ( pseudo_p(i,k) - pseudo_p(i-1,k) ) / dx
            w(i,k) = w(i,k) - dt * ( pseudo_p(i,k) - pseudo_p(i,k-1) ) / dz

            p(i,k) = p(i,k) + pseudo_p(i,k)

        enddo
    enddo
    !$omp end parallel do


    ! Additional Laplacian contribution from implicit treatment
    if (implicitmode .eqv. .true.) then
        half_nudt = 0.5 * nu * dt

        !$omp parallel do &
        !$omp default(none) &
        !$omp private(i,k) &
        !$omp shared(p,pseudo_p,dx,dz,Nx,Nz,dt,half_nudt)
        do k = 1,Nz
            do i = 1,Nx
                p(i,k) = p(i,k) - half_nudt * ( ( pseudo_p(i-1,k) - 2.0 * pseudo_p(i,k) + pseudo_p(i+1,k) ) / dx**2 &
                                              + ( pseudo_p(i,k-1) - 2.0 * pseudo_p(i,k) + pseudo_p(i,k+1) ) / dz**2 )
            enddo
        enddo
        !$omp end parallel do
    endif

    call update_ghost_pressure(p)

end subroutine projectionUpdate