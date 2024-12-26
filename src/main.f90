!----OMP-accelerated 2D Navier--Stokes finite difference solver on a staggered grid -----
! USAGE: Go to params.f90 and modify parameters module as desired

program main
    use parameters
    use fftw3
    use omp_lib
    use grid
    use velfields
    implicit none
    integer :: i
    real :: time
  

    write(*,*) "Executing solver with ", OMP_GET_NUM_THREADS() , " threads!"


    ! Variable initialisation
    call allocFields
    call generate_grid
    call init_fft

    call initialCondition
    time = 0.0d0

    !--------- Begin time-marching ------------------
    do i = 1,nt
        call rhsVelocity ! Build RHS vector for velocity solution
        call update_velocity ! Calculate intermediate velocity, ustar
        call pressurePoisson ! Build div(ustar), solve Poisson, projection update to n+1

        ! call postpro
        ! call check
        time = time + dt
        !write(*,*) u

        if ( mod(i,traw) .eq. 0 ) then
            write(*,*) i

            !call writeVTK(u(1:Nx,1:Nz), Nx, Nz, dx, dz,'u',i)

            call write2DField(u(1:Nx,1:Nz),Nx,Nz,'u',i)
            call write2DField(w(1:Nx,1:Nz),Nx,Nz,'w',i)
            call write2DField(p(1:Nx,1:Nz),Nx,Nz,'p',i)

        endif
    enddo
    !--------- End time-marching -------------------

    ! Memory deallocation
    call deallocFields
    call dealloc_grid
    call dealloc_fft

  end program main