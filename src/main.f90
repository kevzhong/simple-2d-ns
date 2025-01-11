!----OMP-accelerated 2D Navier--Stokes finite difference solver on a staggered grid -----
! USAGE: Go to params.f90 and modify parameters module as desired

program main
    use parameters
    use rk3
    use fftw3
    use omp_lib
    use grid
    use velfields
    use scalarfields
    use velMemory
    implicit none
    integer :: i,nrk3
    real :: time
  
    ! For timing
    real :: start_time, end_time

    !$omp parallel
    !$omp master
    write(*,*) "Executing solver with ", OMP_GET_NUM_THREADS() , " threads!"
    !$omp end master
    !$omp end parallel


    ! Variable initialisation
    call allocFields
    call generate_grid
    call init_fft

    call initialCondition
    time = 0.0d0

    !--------- Begin time-marching ------------------
    do i = 1,nt
        call cpu_time(start_time)

        do nrk3 = 1,3
            call next_rk3(nrk3)
            call updateFields
            call pressurePoisson ! Build div(ustar), solve Poisson, projection update to n+1
            ! call check
        enddo

        call cpu_time(end_time)

        ! call postpro

        time = time + dt

        if ( mod(i,traw) .eq. 0 ) then
            write(*,*) "Timestep ", i," CPU-time per step = ", end_time - start_time

            call write2DField(u(1:Nx,1:Nz),Nx,Nz,dx,dz(1),'u',i)
            call write2DField(w(1:Nx,1:Nz),Nx,Nz,dx,dz(1),'w',i)
            !call write2DField(p(1:Nx,1:Nz),Nx,Nz,dx,dz,'p',i)
            if (scalarmode .eqv. .true.) call write2DField(temp(1:Nx,1:Nz),Nx,Nz,dx,dz(1),'c',i)
        endif
    enddo
    !--------- End time-marching -------------------

    ! Memory deallocation
    call deallocFields
    call dealloc_grid
    call dealloc_fft

  end program main