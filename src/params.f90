!----- PARAMETERS: MODIFY TO CHANGE SIMULATION SETUP ----------------
module parameters
    implicit none

    ! Grid points
    integer :: Nx = 256
    integer :: Nz = 256

    ! Domain size
    !real :: Lx = 1.0
    !real :: Lz = 1.0
    
    real :: Lx = 2.0 * 3.141592653589793
    real :: Lz = 2.0 * 3.141592653589793


    ! Time-stepping
    integer :: Nt = 100000 ! No. of timesteps
    real :: dt = 1.0e-5 ! Timestep

    real :: nu = 1.0d0

    ! Data writing
    integer :: traw = 1000

end module parameters

!------------------ END PARAMETER SETUP MODIFICATION ----------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------



! Computational grid
module grid
    implicit none

    real :: dx, dz
    real, allocatable, dimension(:) :: xc,xm
    real, allocatable, dimension(:) :: zc,zm

end module grid


! Field variables
module velfields
    implicit none
    real, allocatable, dimension(:,:) :: u,w,p
end module velfields

! Working memory for velocity fields
module velMemory
    implicit none

    ! Work arrays for RHSs
    real, allocatable, dimension(:,:) :: rhs_u, rhs_w
    !real, allocatable, dimension(:,:) :: adv_u, adv_w
end module velMemory

! Working memory for FFT / Poisson solver
module fftMemory
    use fftw3
    implicit none
    ! Work arrays for Helmholtz/FFT solver
    real, allocatable, dimension(:) :: lmb_x_on_dx2
    real, allocatable, dimension(:) :: lmb_z_on_dz2

    type(C_PTR) :: fftw_plan_fwd
    type(C_PTR) :: fftw_plan_bwd
    type(C_PTR) :: ptr1, ptr2, ptr3, ptr4

    complex(C_DOUBLE_COMPLEX), pointer :: rhs_hat(:,:), pseudo_phat(:,:)
    real(C_DOUBLE), pointer :: pseudo_p(:,:), rhs_poisson(:,:) ! Poisson equation

end module fftMemory