!----- PARAMETERS: MODIFY TO CHANGE SIMULATION SETUP ----------------
module parameters
    implicit none

    ! Grid points
    integer :: Nx = 256
    integer :: Nz = 128

    ! Domain size
    real :: Lx = 2.0
    real :: Lz = 1.0
    
    !real :: Lx = 2.0 * 3.141592653589793
    !real :: Lz = 2.0 * 3.141592653589793

    ! Time-stepping
    integer :: Nt = 400000 ! No. of timesteps
    real :: dt = 1.0e-4 ! Timestep

    real :: nu = 1.0 / 10000.0
    real :: mean_dpdx = 0.0


    logical :: implicitmode = .true.

    ! Add-ons
    logical :: scalarmode = .true.

    real :: prandtl = 1.0
    real :: beta_gx = 0.0 ; real :: beta_gz = 1.0 ! buoyancy forcing


    ! Wall boundary-conditions
    real :: utop = 0.0
    real :: ubot = 0.0
    real :: wtop = 0.0
    real :: wbot = 0.0
    
    real :: Ttop = 0.0
    real :: Tbot = 1.0

    ! Data writing
    integer :: traw = 2000

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
    real, allocatable, dimension(:,:) :: pseudo_p

end module velfields

module scalarfields
    implicit none
    real, allocatable, dimension(:,:) :: temp
end module scalarfields

! Working memory for velocity fields
module velMemory
    implicit none

    ! Work arrays for RHSs
    real, allocatable, dimension(:,:) :: rhs_u, rhs_w
    !real, allocatable, dimension(:,:) :: adv_u, adv_w
end module velMemory

module scalarMemory
    implicit none
    ! Work arrays for RHSs
    real, allocatable, dimension(:,:) :: rhs_temp
end module scalarMemory

! Working memory for FFT / Poisson solver
module fftMemory
    use fftw3
    implicit none
    ! Work arrays for Helmholtz/FFT solver
    real, allocatable, dimension(:) :: lmb_x_on_dx2

    type(C_PTR) :: fftw_plan_fwd
    type(C_PTR) :: fftw_plan_bwd
    type(C_PTR) :: ptr1, ptr2, ptr3, ptr4

    complex(C_DOUBLE_COMPLEX), pointer :: rhs_hat(:,:), pseudo_phat(:,:)
    real(C_DOUBLE), pointer :: soln_poisson(:,:), rhs_poisson(:,:) ! Poisson equation

end module fftMemory

! Working memory for implicit solve
module implicit
    implicit none
    
    ! TDM coefficients
    real, allocatable, dimension(:) :: ami(:), aci(:), api(:)
    real, allocatable, dimension(:) :: amk(:), ack(:), apk(:)

    ! RHS/solution buffers for TDM
    real, allocatable, dimension(:) :: tdm_rhsX1(:), tdm_rhsX2(:) ! Two buffers for two solves in Sherman--Morrison
    real, allocatable, dimension(:) :: tdm_rhsZ(:)

    ! Provisional solution array (delta u etc to pass in x->y->z)
    real, allocatable, dimension(:,:) :: impl_delta

end module implicit