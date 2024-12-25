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
    use fftMemory
    use parameters
    implicit none
    integer :: i, k
    real :: wavenum_sq

    ! Perform single FFT 2D R2C transform
    ! Not OMP-accelerated but OK for now
    call dfftw_execute_dft_r2c(fftw_plan_fwd, rhs_poisson(:,:), rhs_hat(:,:))

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k,wavenum_sq) &
    !$omp shared(pseudo_phat,rhs_hat,Nx,Nz,lmb_x_on_dx2,lmb_z_on_dz2)
    do k = 1, Nz
        do i = 1, Nx/2+1
            wavenum_sq = lmb_x_on_dx2(i) + lmb_z_on_dz2(k)
            pseudo_phat(i,k) = rhs_hat(i,k) / wavenum_sq / dble(Nx * Nz)
        enddo
    enddo
    !$omp end parallel do

    pseudo_phat(1,1) = 0.0 ! Arbitrary

    call dfftw_execute_dft_c2r(fftw_plan_bwd, pseudo_phat(:,:), pseudo_p(:,:))
    
end subroutine solve_pressurePoisson

subroutine projectionUpdate
     use velfields, only: u,w,p
     use grid
     use fftMemory, only: pseudo_p
     use ghost
     use parameters
     implicit none
     integer :: i, k

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

     call update_ghost(u)
     call update_ghost(w)
     call update_ghost(p)

end subroutine projectionUpdate