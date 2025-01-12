! Solve for the provisional velocity ustar (or scalar at new time-step) implicitly by solving a 1D Helmholtz equation
!
! [              (   2         2   )  ]
! [              (  d         d    )  ]
! [  p    - q    ( ----   +  ----  )  ] DT      =   rhs
! [              (    2         2  )  ]   ik           ik
! [              (  dx        dz   )  ]
!
!
!                n+1        n
! where DT    = T     -   T            and T can be either {u,v,w,temp}, etc
!         ik     ik        ik
!
! Application of an FFT in the x-periodic direction yields the 1D Helmholtz equation to solve
! 
! 
! [              (  \Lambda          2   )  ]  
! [              (          l       d    )  ]  ^            ^
! [  p    - q    ( ----------   +  ----  )  ]  DT      =   rhs
! [              (      2             2  )  ]    lk           lk
! [              (    DX            dz   )  ] 
!
!
! For each x-wavenumber l
subroutine solve_helmholtz(field,rhs,p,q,bc_type_bot,bc_type_top)
    use fftw3
    use grid
    use ghost
    use fftMemory
    use parameters
    use bctypes
    use implicit
    implicit none
    integer :: i, k
    real, intent(inout) :: field(Nx,Nz)
    real :: rhs(Nx,Nz)
    real :: p, q
    real :: dzmh, dzph
    integer, intent(in) :: bc_type_bot, bc_type_top



    ! Tri-diagonal inversion for each kx wavenumber

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(k) &
    !$omp shared(fftw_plan_fwd,rhs,rhs_hat,Nz)
    do k = 1,Nz
        ! FFT in x direction for each z-location k
        call fftw_execute_dft_r2c(fftw_plan_fwd, rhs(:,k), rhs_hat(:,k))
    enddo
    !$omp end parallel do

    !a = -q  / dz**2
    !c = -q  / dz**2

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k,amk,ack,apk,tdm_rhsZ_r,tdm_rhsZ_c,dzmh,dzph) &
    !$omp shared(p,q,lmb_x_on_dx2,rhs_hat,soln_hh_hat,dz,Nx,Nz,bc_type_bot,bc_type_top)
    do i = 1,Nx/2+1

        !b = p - q * ( lmb_x_on_dx2(i) - 2.0 / dz**2 )

        do k = 1,Nz

            rhs_hat(i,k) = rhs_hat(i,k) / Nx
            !amk(k) = a
            !ack(k) = b
            !apk(k) = c

            dzmh = 0.5*( dz(k-1) + dz(k  ) )
            dzph = 0.5*( dz(k  ) + dz(k+1) )

            amk(k) = -q / (dz(k) * dzmh )
            ack(k) = p - q * (lmb_x_on_dx2(i) - (dzmh + dzph) / ( dz(k) * dzmh * dzph )  )
            apk(k) = -q / (dz(k) * dzph )

            tdm_rhsZ_r(k) = real(rhs_hat(i,k))
            tdm_rhsZ_c(k) = aimag(rhs_hat(i,k))
        enddo

        !!!!!!!!   BOUNDARY CONDITIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Arbitrary Dirichlet for 0 mode when considering pressure-Poisson
        if (bc_type_bot .eq. PRESSUREBC) then
            if (i .eq. 1) then
                apk(1) = 0.0
                ack(1) = 1.0
                tdm_rhsZ_r(1) = 0.0
                tdm_rhsZ_c(1) = 0.0
            endif
        endif

        ! Dirichlet wall conditions
        if (bc_type_bot .eq. DIRICHLET) then
            ack(1) = ack(1) - amk(1)
            !apk(1) = c
        endif

        if (bc_type_top .eq. DIRICHLET) then
            ack(Nz) = ack(Nz) - apk(Nz) ! KZ: check this later
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call tridiag(amk,ack,apk,tdm_rhsZ_r,Nz)
        call tridiag(amk,ack,apk,tdm_rhsZ_c,Nz)

        do k = 1,Nz
            soln_hh_hat(i,k) = CMPLX( tdm_rhsZ_r(k), tdm_rhsZ_c(k) )
        enddo
    enddo
    !$omp end parallel do


    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k) &
    !$omp shared(fftw_plan_bwd,soln_hh_hat,soln_hh,Nz)
    do k = 1,Nz
        call fftw_execute_dft_c2r(fftw_plan_bwd, soln_hh_hat(:,k), soln_hh(:,k))
    enddo
    !$omp end parallel do


    if (bc_type_bot .eq. PRESSUREBC) then ! Pseudo-pressure Poisson solve, soln field is not a delta vector

        !$omp parallel do &
        !$omp default(none) &
        !$omp private(i,k) &
        !$omp shared(soln_hh,field,Nx,Nz)
        do k = 1,Nz
            do i = 1,Nx
                field(i,k) = soln_hh(i,k)
            enddo
        enddo
        !$omp end parallel do

    else ! Solution was Delta form for {DU,DT} etc

        !$omp parallel do &
        !$omp default(none) &
        !$omp private(i,k) &
        !$omp shared(soln_hh,field,Nx,Nz)
        do k = 1,Nz
            do i = 1,Nx
                field(i,k) = field(i,k) + soln_hh(i,k)
            enddo
        enddo
        !$omp end parallel do

    endif


end subroutine solve_helmholtz

! Special stencil for the staggered w-cells
subroutine solve_helmholtzW(field,rhs,p,q,bc_type_bot,bc_type_top)
    use fftw3
    use grid
    use ghost
    use fftMemory
    use parameters
    use bctypes
    use implicit
    implicit none
    integer :: i, k
    real, intent(inout) :: field(Nx,Nz)
    real :: rhs(Nx,Nz)
    real :: p, q
    real :: dzmh, dzph
    integer, intent(in) :: bc_type_bot, bc_type_top



    ! Tri-diagonal inversion for each kx wavenumber

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(k) &
    !$omp shared(fftw_plan_fwd,rhs,rhs_hat,Nz)
    do k = 1,Nz
        ! FFT in x direction for each z-location k
        call fftw_execute_dft_r2c(fftw_plan_fwd, rhs(:,k), rhs_hat(:,k))
    enddo
    !$omp end parallel do

    !a = -q  / dz**2
    !c = -q  / dz**2

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k,amk,ack,apk,tdm_rhsZ_r,tdm_rhsZ_c,dzmh,dzph) &
    !$omp shared(p,q,lmb_x_on_dx2,rhs_hat,soln_hh_hat,dz,Nx,Nz,bc_type_bot,bc_type_top)
    do i = 1,Nx/2+1

        !b = p - q * ( lmb_x_on_dx2(i) - 2.0 / dz**2 )

        do k = 1,Nz

            rhs_hat(i,k) = rhs_hat(i,k) / Nx
            !amk(k) = a
            !ack(k) = b
            !apk(k) = c

            dzmh = 0.5*( dz(k-1) + dz(k  ) )
            dzph = 0.5*( dz(k  ) + dz(k+1) )

            amk(k) = -q / (dz(k-1) * dzmh )
            ack(k) = p - q * (lmb_x_on_dx2(i) - (dz(k) + dzmh) / ( dz(k-1) * dz(k) * dzmh )  )
            apk(k) = -q / (dz(k) * dzmh )

            tdm_rhsZ_r(k) = real(rhs_hat(i,k))
            tdm_rhsZ_c(k) = aimag(rhs_hat(i,k))
        enddo

        !!!!!!!!   BOUNDARY CONDITIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Dirichlet wall conditions
        if (bc_type_bot .eq. DIRICHLET) then
            ack(1) = ack(1) - amk(1)
            !apk(1) = c
        endif

        if (bc_type_top .eq. DIRICHLET) then
            ack(Nz) = ack(Nz) - apk(Nz) ! KZ: check this later
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call tridiag(amk,ack,apk,tdm_rhsZ_r,Nz)
        call tridiag(amk,ack,apk,tdm_rhsZ_c,Nz)

        do k = 1,Nz
            soln_hh_hat(i,k) = CMPLX( tdm_rhsZ_r(k), tdm_rhsZ_c(k) )
        enddo
    enddo
    !$omp end parallel do


    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k) &
    !$omp shared(fftw_plan_bwd,soln_hh_hat,soln_hh,Nz)
    do k = 1,Nz
        call fftw_execute_dft_c2r(fftw_plan_bwd, soln_hh_hat(:,k), soln_hh(:,k))
    enddo
    !$omp end parallel do


    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k) &
    !$omp shared(soln_hh,field,Nx,Nz)
    do k = 1,Nz
        do i = 1,Nx
            field(i,k) = field(i,k) + soln_hh(i,k)
        enddo
    enddo
    !$omp end parallel do



end subroutine solve_helmholtzW


 !subroutine apply_bc_helmholtz(i,Nz,bc_type_bot, bc_type_top, rbuffer, cbuffer, amk, ack, apk,a,b,c)
 !    implicit none
 !end subroutine apply_bc_helmholtz