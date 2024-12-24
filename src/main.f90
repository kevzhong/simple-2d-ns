!----OMP-accelerated 2D Navier--Stokes finite difference solver on a staggered grid -----
! USAGE: Go to params.f90 and modify parameters module as desired

program main
    use parameters, only: nt, dt
    use fftw3
    use omp_lib
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
        write(*,*) time
        write(*,*) u
    enddo
    !--------- End time-marching -------------------

    ! Memory deallocation
    call deallocFields
    call dealloc_grid
    call dealloc_fft

  end program main