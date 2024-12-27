subroutine init_fft
    use parameters, only: Nx, Nz
    use grid, only : dx
    use fftw3
    use fftMemory
    use velMemory
    implicit none

    integer :: i
    real :: PI = 2.d0*dasin(1.d0) 


    !---------------------------- Modified wavenumbers --------------------------------
    allocate( lmb_x_on_dx2(Nx/2+1) )

    do i = 1,Nx/2+1
        lmb_x_on_dx2(i) = ( 2.0 * cos(2.0*PI*(real(i)-1.0) / Nx ) - 2.0 ) / dx**2
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

    ! 1D real to complex plan:  double(Nx,1) ---> complex(Nx/2+1, 1)
    fftw_plan_fwd = fftw_plan_dft_r2c_1d(Nx, rhs_poisson(:,1), rhs_hat(:,1), FFTW_ESTIMATE)

    ! 1D complex to real transform: complex(Nx/2+1, 1) ----> double(Nx,1)
    fftw_plan_bwd = fftw_plan_dft_c2r_1d(Nx, rhs_hat(:,1), pseudo_p(:,1), FFTW_ESTIMATE)


end subroutine init_fft


subroutine dealloc_fft
    use fftMemory
    implicit none

    if(allocated(lmb_x_on_dx2)) deallocate(lmb_x_on_dx2)

    ! Destroy the FFTW plan
    call fftw_destroy_plan(fftw_plan_fwd)
    call fftw_destroy_plan(fftw_plan_bwd)
    call fftw_free(ptr1)
    call fftw_free(ptr2)
    call fftw_free(ptr3)
    call fftw_free(ptr4)

end subroutine dealloc_fft
