! Solve for the provisional velocity ustar (or scalar at new time-step) implicitly using the alternating direction implicit (ADI) / approximate factorization method
subroutine ADI_implicitUpdate
    use rk3
    use parameters
    use velfields
    use velMemory
    use scalarfields
    use scalarMemory
    use ghost
    use grid
    use implicit
    implicit none
    real :: lapl_prefac

    !--------- u velocity --------------------------------------
    lapl_prefac = 0.5 * nu * aldt / dx**2
    call ADI_periodicSolveX(lapl_prefac,rhs_u)
    lapl_prefac = 0.5 * nu * aldt / dz**2
    call ADI_wallSolveZ(lapl_prefac,impl_delta(:,:),u(1:Nx,1:Nz),bctype_ubot,bctype_utop)
    call update_ghost_wallsU(u,bctype_ubot,bctype_utop,bcval_ubot,bcval_utop)

    !--------- w velocity --------------------------------------
    lapl_prefac = 0.5 * nu * aldt / dx**2
    call ADI_periodicSolveX(lapl_prefac,rhs_w)
    lapl_prefac = 0.5 * nu * aldt / dz**2
    call ADI_wallSolveZ(lapl_prefac,impl_delta(:,:),w(1:Nx,1:Nz),bctype_wbot,bctype_wtop)
    call update_ghost_wallsW(w,bctype_wbot,bctype_wtop,bcval_wbot,bcval_wtop)


    !--------- scalar -------------------------------------------
    if (scalarmode .eqv. .true.) then
        lapl_prefac = 0.5 * nu/prandtl * aldt / dx**2
        call ADI_periodicSolveX(lapl_prefac,rhs_temp)
        lapl_prefac = 0.5 * nu/prandtl * aldt / dz**2
        call ADI_wallSolveZ(lapl_prefac,impl_delta(:,:),temp(1:Nx,1:Nz),bctype_Tbot,bctype_Ttop)

        call update_ghost_wallTemp(temp,bctype_Tbot,bctype_Ttop,bcval_Tbot,bcval_Ttop)
    endif

end subroutine ADI_implicitUpdate

subroutine ADI_periodicSolveX(half_nualdt_on_dx2,rhs)
    use velfields
    use parameters
    use velMemory
    use ghost
    use implicit
    implicit none
    integer :: i, k
    real, intent(in) :: half_nualdt_on_dx2
    real, dimension(Nx,Nz), intent(in) :: rhs
    real :: a,b,c,d

    ! For Sherman--Morrison
    real :: vT_rhsX1, vT_rhsX2, const


    ! TRIDIAGONAL SYSTEM TO INVERT
    ! For { du_i-1, du_i, du_i+1 }
    !
    !
    !  [             ]           [           ]         [             ]           
    !  |      -G     | du      + |  1 + 2*G  | du   +  |      -G     | du      =  RHS
    !  [             ]  i-1      [           ]  i      [             ]  i+1          i

    ! This is normalised to obtain unity diagonals, helping with diagonal dominance and conditioning:
    !
    !
    !
    !  [       G     ]           [     ]         [       G     ]           [       1     ]
    !  | -  -------- | du      + |  1  | du   +  | -  -------- | du      = |   --------- | RHS
    !  [    1 + 2*G  ]  i-1      [     ]  i      [    1 + 2*G  ]  i+1      [    1 + 2*G  ]    i


    ! Here, G := half_nualdt_on_dx2 = nu * dt / (0.5 * dx^2)
    !      du := ustar - u_n

    a = -half_nualdt_on_dx2 / (1.0 + 2.0 * half_nualdt_on_dx2)
    b = 1.0
    c = a
    d = 1.0 / (1.0 + 2.0 * half_nualdt_on_dx2) ! Normalisation factor for RHS


    ! Coefficients with no normalisation
    !a = -half_nualdt_on_dx2
    !b = 1.0 + 2.0 * half_nualdt_on_dx2
    !c = -half_nualdt_on_dx2
    !d = 1.0

    ! Sherman--Morrison pseudo-code
    !
    ! 0) Construct B matrix from A: perturb to get diagonal system
    ! 1) TDM solve B * rhsX1 = rhsX1 (By = rhs in original notation)
    ! 2) TDM solve B * rhsX2 = rhsX2 (Bq = u in original notation)
    !       Here, rhsX2 (u in original notation) is [-aci_1, 0, 0, ...., api_N] = [-b, 0, 0, ...., c]
    !             rhsX1 is solved with the original rhs vector of the problem
    ! 3) Compute vT_rhsX1 = rhsX1(1) - ami(1)/aci(1) * rhsX1(N)     (vTy in original notation)
    ! 4) Compute vT_rhsX2 = rhsX2(1) - ami(1)/aci(1) * rhsX2(N)     (vTq in original notation)
    ! 5) Final solution to periodic tridiag problem is:
    !       du = rhsX1 - const * rhsX2
    !           where const = vT_rhsX1 / (1 + vT_rhsX2)

    do i = 1,Nx
        ami(i) = a
        aci(i) = b
        api(i) = c
    enddo

    ! Perturbed matrix system: Sherman Morrison: A-->B
    aci(1) = 2.0 * b
    aci(Nx) = b + a*c/b 

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k,tdm_rhsX1,tdm_rhsX2,vT_rhsX1,vT_rhsX2,const) &
    !$omp shared(Nx,Nz,a,b,c,d,rhs,impl_delta,ami,aci,api)
    do k = 1,Nz

        ! Reset RHS
        do i = 1,Nx
            tdm_rhsX1(i) = rhs(i,k) * d
            tdm_rhsX2(i) = 0.0
        enddo
        tdm_rhsX2(1) = -b
        tdm_rhsX2(Nx) = c
        
        call tridiag(ami,aci,api,tdm_rhsX1,Nx) ! Step 1)
        call tridiag(ami,aci,api,tdm_rhsX2,Nx) ! Step 2)

        vT_rhsX1 =  tdm_rhsX1(1) - a/b * tdm_rhsX1(Nx) ! Step 3)
        vT_rhsX2 =  tdm_rhsX2(1) - a/b * tdm_rhsX2(Nx) ! Step 4)

        if ( abs(1.0 + vT_rhsX2) .lt. EPSILON(1.0d0) ) then
            const =  ( vT_rhsX1 / ( 1.0 + vT_rhsX2 + EPSILON(1.0d0) ) ) 
            write(*,*) " Singularity detected in SM"
        else
            const =  ( vT_rhsX1 / ( 1.0 + vT_rhsX2 ) ) 
        endif

        do i = 1,Nx
            impl_delta(i,k) = tdm_rhsX1(i) - const * tdm_rhsX2(i)
        enddo

    enddo
    !$omp end parallel do



end subroutine ADI_periodicSolveX

subroutine ADI_wallSolveZ(half_nualdt_on_dz2,rhs,field,bc_type_bot,bc_type_top)
    use velfields
    use bctypes
    use parameters
    use velMemory
    use ghost
    use implicit
    implicit none
    integer :: i, k
    real, intent(in) :: half_nualdt_on_dz2
    real, dimension(Nx,Nz), intent(in) :: rhs
    real, dimension(Nx,Nz), intent(inout) :: field
    integer :: bc_type_bot, bc_type_top
    real :: a,b,c,d, d_bcb, d_bct

    a = -half_nualdt_on_dz2 / (1.0 + 2.0 * half_nualdt_on_dz2)
    b = 1.0
    c = a
    d = 1.0 / (1.0 + 2.0 * half_nualdt_on_dz2) ! Normalisation factor for RHS

    do k = 1,Nz
        amk(k) = a
        ack(k) = b
        apk(k) = c
    enddo

    if (bc_type_bot .eq. DIRICHLET) then
        d_bcb = 1.0 / (1.0 + 3.0 * half_nualdt_on_dz2)
        !ack(1) = 1.0
        apk(1) = -half_nualdt_on_dz2 * d_bcb ! Bottom wall Dirichlet
    else
        d_bcb = 1.0 / (1.0 + 3.0 * half_nualdt_on_dz2) ! DUMMY TO AVOID WARNING
    endif

    if (bc_type_top .eq. DIRICHLET) then
        d_bct = 1.0 / (1.0 + 3.0 * half_nualdt_on_dz2)
        amk(Nz-1) = -half_nualdt_on_dz2 * d_bct ! Top wall Dirichlet
    else
        d_bct = 1.0 / (1.0 + 3.0 * half_nualdt_on_dz2) ! DUMMY TO AVOID WARNING
    endif

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,k,tdm_rhsZ_r) &
    !$omp shared(Nx,Nz,d,rhs,field,amk,ack,apk,d_bcb,d_bct)
    do i = 1,Nx
        
        ! Reset RHS
        do k = 2,Nz-1
            tdm_rhsZ_r(k) = rhs(i,k) * d
        enddo
        ! Boundary conditions
        tdm_rhsZ_r(1) = ( rhs(i,1)  ) * d_bcb
        tdm_rhsZ_r(Nz) = ( rhs(i,Nz) ) * d_bct


        call tridiag(amk,ack,apk,tdm_rhsZ_r,Nz) ! Step 1)

        do k = 1,Nz
            field(i,k) = field(i,k) + tdm_rhsZ_r(k)
        enddo

    enddo
    !$omp end parallel do


end subroutine ADI_wallSolveZ
