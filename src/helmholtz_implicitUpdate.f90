
subroutine helmholtz_implicitUpdate
    use parameters
    use rk3
    use velfields
    use velMemory
    use scalarfields
    use scalarMemory
    use ghost
    use grid
    use implicit
    use bctypes
    implicit none
    real :: lapl_prefac

    if (implicitXmode .eqv. .true. ) then
        !--------- u velocity --------------------------------------
        lapl_prefac =  0.5 * nu * aldt 
        call solve_helmholtz(u(1:Nx,1:Nz), rhs_u, 1.0 , lapl_prefac, bctype_ubot, bctype_utop )
        call update_ghost_wallsU(u,bctype_ubot,bctype_utop,bcval_ubot,bcval_utop)

        !--------- w velocity --------------------------------------
        call solve_helmholtz(w(1:Nx,1:Nz), rhs_w, 1.0 , lapl_prefac, bctype_wbot, bctype_wtop )
        call update_ghost_wallsW(w,bctype_wbot,bctype_wtop,bcval_wbot,bcval_wtop)

        !--------- scalar -------------------------------------------
        if (scalarmode .eqv. .true.) then
            lapl_prefac = 0.5 * nu/prandtl * aldt
            call solve_helmholtz(temp(1:Nx,1:Nz), rhs_temp, 1.0 , lapl_prefac, bctype_Tbot, bctype_Ttop )
            call update_ghost_wallTemp(temp,bctype_Tbot,bctype_Ttop,bcval_Tbot,bcval_Ttop)
        endif

    else ! No implicit in wall-parallel is just a single TDMA solve in z: no FFT needed
        !--------- u velocity --------------------------------------
        lapl_prefac = 0.5 * nu * aldt / dz**2
        call ADI_wallSolveZ(lapl_prefac,rhs_u(:,:),u(1:Nx,1:Nz),bctype_ubot,bctype_utop)
        call update_ghost_wallsU(u,bctype_ubot,bctype_utop,bcval_ubot,bcval_utop)

        !--------- w velocity --------------------------------------
        lapl_prefac = 0.5 * nu * aldt / dz**2
        call ADI_wallSolveZ(lapl_prefac,rhs_w(:,:),w(1:Nx,1:Nz),bctype_wbot,bctype_wtop)
        call update_ghost_wallsW(w,bctype_wbot,bctype_wtop,bcval_wbot,bcval_wtop)

        if (scalarmode .eqv. .true.) then
            !--------- scalar -------------------------------------------
            lapl_prefac = 0.5 * nu/prandtl * aldt / dz**2
            call ADI_wallSolveZ(lapl_prefac,rhs_temp(:,:),temp(1:Nx,1:Nz),bctype_Tbot,bctype_Ttop)
            call update_ghost_wallTemp(temp,bctype_Tbot,bctype_Ttop,bcval_Tbot,bcval_Ttop)
        endif

    endif




end subroutine helmholtz_implicitUpdate
