subroutine init_fft
    use parameters, only: Nx, Nz
    use grid, only : dx, dz
    use fftw3
    use fftMemory
    use velMemory
    implicit none

    integer :: i,k
    real :: PI = 2.d0*dasin(1.d0) 


    !---------------------------- Modified wavenumbers --------------------------------
    allocate( lmb_x_on_dx2(Nx/2+1) )
    allocate( lmb_z_on_dz2(  Nz  ) )

    do i = 1,Nx/2+1
        lmb_x_on_dx2(i) = ( 2.0 * cos(2.0*PI*(real(i)-1.0) / Nx ) - 2.0 ) / dx**2
    enddo

    do k = 1,Nz
        lmb_z_on_dz2(k) = ( 2.0 * cos(2.0*PI*(real(k)-1.0) / Nz ) - 2.0 ) / dz**2
    enddo

    !--------------- FFT plans for Helmholtz solver using C-type allocations --------

    ! Solution array real(Nx,Nz)
    ptr1 = fftw_alloc_real(int(Nx * Nz, C_SIZE_T))
    call c_f_pointer(ptr1, pseudo_p, [Nx,Nz])

    ptr2 = fftw_alloc_complex(int((Nx/2+1) * Nz, C_SIZE_T))
    call c_f_pointer(ptr2, pseudo_phat, [Nx/2+1,Nz])

    ptr3 = fftw_alloc_complex(int((Nx/2+1) * Nz, C_SIZE_T))
    call c_f_pointer(ptr3, rhs_hat, [Nx/2+1,Nz])

    ptr4 = fftw_alloc_real(int(Nx * Nz, C_SIZE_T))
    call c_f_pointer(ptr4, rhs_poisson, [Nx,Nz])

    ! 2D real to complex plan:  double(Nx,Nz) ---> complex(Nx/2+1, Nz)
    fftw_plan_fwd = fftw_plan_dft_r2c_2d(Nx, Nz, rhs_poisson(:,:), rhs_hat(:,:), FFTW_ESTIMATE)

    ! 2D complex to real transform: complex(Nx/2+1, Nz) ----> double(Nx,Nz)
    fftw_plan_bwd = fftw_plan_dft_c2r_2d(Nx, Nz, rhs_hat(:,:), pseudo_p(:,:), FFTW_ESTIMATE)


end subroutine init_fft


subroutine dealloc_fft
    use fftMemory
    implicit none

    if(allocated(lmb_x_on_dx2)) deallocate(lmb_x_on_dx2)
    if(allocated(lmb_z_on_dz2)) deallocate(lmb_z_on_dz2)


    ! Destroy the FFTW plan
    call fftw_destroy_plan(fftw_plan_fwd)
    call fftw_destroy_plan(fftw_plan_bwd)
    call fftw_free(ptr1)
    call fftw_free(ptr2)
    call fftw_free(ptr3)
    call fftw_free(ptr4)

end subroutine dealloc_fft
