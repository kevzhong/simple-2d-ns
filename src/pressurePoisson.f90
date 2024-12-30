subroutine pressurePoisson
    implicit none

    call build_rhsPoisson
    call solve_pressurePoisson    
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

subroutine solve_pressurePoisson
    use fftw3
    use velfields, only: pseudo_p
    use grid
    use ghost
    use fftMemory
    use parameters
    implicit none
    integer :: i, k
    real(8), allocatable :: am(:), ac(:), ap(:) ! tridiagonal coefficients
    real(8), allocatable :: rbuffer(:), cbuffer(:) ! solution vector to rhs
    real(8) :: a, b, c

    allocate( am(1:Nz) ) ; allocate( ac(1:Nz) ) ; allocate( ap(1:Nz) ) ; allocate( rbuffer(1:Nz) ) ; allocate( cbuffer(1:Nz) )

    ! Tri-diagonal inversion for each kx wavenumber

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(k) &
    !$omp shared(fftw_plan_fwd,rhs_poisson,rhs_hat,Nz)
    do k = 1,Nz
        ! FFT in x direction for each z-location k
        call fftw_execute_dft_r2c(fftw_plan_fwd, rhs_poisson(:,k), rhs_hat(:,k))
    enddo
    !$omp end parallel do

    ! Build and solve tri-diagonal system for each kx wavenumber

    a = 1.0 / dz**2
    c = 1.0 / dz**2

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,b,k,am,ac,ap,rbuffer,cbuffer) &
    !$omp shared(a,c,lmb_x_on_dx2,rhs_hat,pseudo_phat,dz,Nx,Nz)
    do i = 1,Nx/2+1

        b = lmb_x_on_dx2(i) - 2.0 / dz**2

        do k = 1,Nz

            rhs_hat(i,k) = rhs_hat(i,k) / Nx

            am(k) = a
            ac(k) = b
            ap(k) = c

            rbuffer(k) = real(rhs_hat(i,k))
            cbuffer(k) = aimag(rhs_hat(i,k))
        enddo

        ! Arbitrary Dirichlet for 0 mode, over-ride
        if (i .eq. 1) then
            ap(1) = 0.0
            ac(1) = 1.0
            rbuffer(1) = 0.0
            cbuffer(1) = 0.0
        endif


        call tridiag(am,ac,ap,rbuffer,Nz)
        call tridiag(am,ac,ap,cbuffer,Nz)

        do k = 1,Nz
            pseudo_phat(i,k) = CMPLX( rbuffer(k), cbuffer(k) )
        enddo
    enddo
    !$omp end parallel do


    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k) &
    !$omp shared(fftw_plan_bwd,pseudo_phat,soln_poisson,Nz)
    do k = 1,Nz
        call fftw_execute_dft_c2r(fftw_plan_bwd, pseudo_phat(:,k), soln_poisson(:,k))
    enddo
    !$omp end parallel do


    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k) &
    !$omp shared(soln_poisson,pseudo_p,Nx,Nz)
    do k = 1,Nz
        do i = 1,Nx
        pseudo_p(i,k) = soln_poisson(i,k)
        enddo
    enddo
    !$omp end parallel do

    call update_ghost_pressure(pseudo_p)

    deallocate(am) ; deallocate(ac) ; deallocate(ap) ; deallocate(rbuffer) ; deallocate(cbuffer)
    
end subroutine solve_pressurePoisson



subroutine projectionUpdate
    use velfields
    use grid
    use ghost
    use parameters
    implicit none
    integer :: i, k
    real :: half_nudt

    half_nudt = 0.5 * nu * dt

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k) &
    !$omp shared(u,w,p,pseudo_p,dx,dz,Nx,Nz,dt, half_nudt)
    do k = 1,Nz
        do i = 1,Nx

           !if (i .eq. 1) then
                !u(i,k) = u(i,k) - dt * ( pseudo_p(i,k) - pseudo_p(Nx,k) ) / dx
            !else
                u(i,k) = u(i,k) - dt * ( pseudo_p(i,k) - pseudo_p(i-1,k) ) / dx
           !endif

            ! Neumann implied dp/dz = dp/dn = 0 for k = 1
            ! w_n+1 = w* = w_wall at boundary
           !if (k .ne. 1) then
               w(i,k) = w(i,k) - dt * ( pseudo_p(i,k) - pseudo_p(i,k-1) ) / dz
           !endif

            p(i,k) = p(i,k) + pseudo_p(i,k)

        enddo
    enddo
    !$omp end parallel do

    !call update_ghost_periodic(u)
    !call update_ghost_periodic(w)
    !call update_ghost_periodic(p)

    !call update_ghost_walls(u,w,ubot,utop,wbot,wtop)

    call update_ghost_pressure(p)

end subroutine projectionUpdate


! subroutine projectionUpdate
!      use velfields, only: u,w,p
!      use grid
!      use fftMemory, only: pseudo_p
!      use ghost
!      use parameters
!      implicit none
!      integer :: i, k
!      real :: half_nudt

!      half_nudt = 0.5 * nu * dt

!      !$omp parallel do &
!      !$omp default(none) &
!      !$omp private(i,k) &
!      !$omp shared(u,w,p,pseudo_p,dx,dz,Nx,Nz,dt, half_nudt)
!      do k = 1,Nz
!          do i = 1,Nx
!             if (i .eq. 1) then
!                 u(i,k) = u(i,k) - dt * ( pseudo_p(i,k) - pseudo_p(Nx,k) ) / dx

!                 !p(i,k) = p(i,k) +  ( pseudo_p(Nx,k) - 2.0 * pseudo_p(i,k) + pseudo_p(i+1,k) ) / dx ! periodic

!                 else
!                 u(i,k) = u(i,k) - dt * ( pseudo_p(i,k) - pseudo_p(i-1,k) ) / dx

!                 !p(i,k) = p(i,k) +  ( pseudo_p(i-1,k) - 2.0 * pseudo_p(i,k) + pseudo_p(i+1,k) ) / dx

!             endif

!              ! Neumann implied dp/dz = dp/dn = 0 for k = 1
!              ! w_n+1 = w* = w_wall at boundary
!             if (k .eq. 1) then
!                 ! w(i,k) = w(i,k)
!                 p(i,k) = p(i,k) +  ( -1.0 * pseudo_p(i,k) + pseudo_p(i,k+1) ) / dz
!             elseif (k .eq. Nz) then
!                 w(i,k) = w(i,k) - dt * ( pseudo_p(i,k) - pseudo_p(i,k-1) ) / dz

!                 p(i,k) = p(i,k) +  ( pseudo_p(i,k-1) - pseudo_p(i,k)  ) / dz

!             else
!                 w(i,k) = w(i,k) - dt * ( pseudo_p(i,k) - pseudo_p(i,k-1) ) / dz

!                 p(i,k) = p(i,k) +  ( pseudo_p(i,k-1) - 2.0 * pseudo_p(i,k) + pseudo_p(i,k+1) ) / dz

!             endif

            
            
!              p(i,k) = p(i,k) + pseudo_p(i,k)

!          enddo
!      enddo
!      !$omp end parallel do

!      !call update_ghost_periodic(u)
!      !call update_ghost_periodic(w)
!      !call update_ghost_periodic(p)

!      !call update_ghost_walls(u,w,ubot,utop,wbot,wtop)

!      call update_ghost_pressure(p)

! end subroutine projectionUpdate