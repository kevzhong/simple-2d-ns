module ghost
    implicit none
    integer :: halosize = 1

    contains

    subroutine update_ghost_wallsU(u,bctype_ubot,bctype_utop,bcval_ubot,bcval_utop)
        use parameters, only: Nx, Nz
        use bctypes
        implicit none
        real(8), intent(inout), dimension(halosize-1:Nx+halosize , halosize-1:Nz+halosize ) :: u
        real(8) :: bcval_ubot, bcval_utop
        integer, intent(in) :: bctype_ubot, bctype_utop

        ! u Bottom wall
        if (bctype_ubot .eq. DIRICHLET) then
            u(:,1-halosize) = 2.0 * bcval_ubot - u(:,1)
        endif

        ! u Top wall
        if (bctype_utop .eq. DIRICHLET) then
            u(:,Nz+halosize) = 2.0 * bcval_utop - u(:,Nz)
        endif

        ! Periodicity in x
        u(1-halosize,:) = u(Nx,:)
        u(Nx+halosize,:) = u(1,:)

    end subroutine update_ghost_wallsU


    subroutine update_ghost_wallsW(w,bctype_wbot,bctype_wtop,bcval_wbot,bcval_wtop)
        use parameters, only: Nx, Nz
        use bctypes
        implicit none
        real(8), intent(inout), dimension(halosize-1:Nx+halosize , halosize-1:Nz+halosize ) :: w
        real(8) ::  bcval_wbot, bcval_wtop
        integer, intent(in) :: bctype_wbot, bctype_wtop

        ! w Bottom wall
        if (bctype_wbot .eq. DIRICHLET) then
            w(:,1) = bcval_wbot ! ghost cell redundant w(:,1-halosize)
        endif

        ! w top wall
        if (bctype_wtop .eq. DIRICHLET) then
            w(:,Nz+halosize) = bcval_wtop
        endif

        ! Periodicity in x

        w(1-halosize,:) = w(Nx,:)
        w(Nx+halosize,:) = w(1,:)

    end subroutine update_ghost_wallsW

    subroutine update_ghost_wallTemp(temp,bctype_Tbot,bctype_Ttop,bcval_Tbot,bcval_Ttop)
        use parameters, only: Nx, Nz
        use bctypes
        implicit none
        real(8), intent(inout), dimension(halosize-1:Nx+halosize , halosize-1:Nz+halosize ) :: temp
        real(8) :: bcval_Tbot, bcval_Ttop
        integer, intent(in) :: bctype_Tbot, bctype_Ttop

        ! Bottom wall
        if (bctype_Tbot .eq. DIRICHLET) then
            temp(:,1-halosize) = 2.0 * bcval_Tbot - temp(:,1)
        endif

        ! Top wall
        if (bctype_Ttop .eq. DIRICHLET) then
            temp(:,Nz+halosize) = 2.0 * bcval_Ttop - temp(:,Nz)
        endif

        ! Periodicity in x
        temp(1-halosize,:) = temp(Nx,:)
        temp(Nx+halosize,:) = temp(1,:)

    end subroutine update_ghost_wallTemp


    subroutine update_ghost_pressure(p)
        use parameters, only: Nx, Nz
        implicit none
        real(8), intent(inout), dimension(halosize-1:Nx+halosize , halosize-1:Nz+halosize ) :: p

        ! Wall BCs: dpdn = 0
        p(:,1-halosize) = p(:,1)
        p(:,Nz+halosize) = p(:,Nz)


        ! Periodicity in x
        p(1-halosize,:) = p(Nx,:)
        p(Nx+halosize,:) = p(1,:)

    end subroutine update_ghost_pressure

    ! subroutine update_ghost_periodic(array)
    !     use parameters, only: Nx, Nz
    !     implicit none
    !     real(8), intent(inout), dimension(halosize-1:Nx+halosize , halosize-1:Nz+halosize ) :: array

    !     ! Periodicity application

    !     ! x-halo swap
    !     array(1-halosize,:) = array(Nx,:)
    !     array(Nx+halosize,:) = array(1,:)

    !     ! z-halo swap
    !     array(:,halosize-1) = array(:,Nz)
    !     array(:,Nz+halosize) = array(:,1)

    ! end subroutine update_ghost_periodic

end module ghost