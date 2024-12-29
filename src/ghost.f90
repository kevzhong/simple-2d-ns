module ghost
    implicit none
    integer :: halosize = 1

    contains

    subroutine update_ghost_walls(u,w,ubot,utop,wbot,wtop)
        use parameters, only: Nx, Nz
        implicit none
        real(8), intent(inout), dimension(halosize-1:Nx+halosize , halosize-1:Nz+halosize ) :: u,w
        real(8) :: ubot, utop, wbot, wtop


        ! Bottom wall
        u(:,1-halosize) = 2.0 * ubot - u(:,1)
        w(:,1) = wbot ! ghost cell redundant w(:,1-halosize)

        ! Top wall
        u(:,Nz+halosize) = 2.0 * utop - u(:,Nz)
        w(:,Nz+halosize) = wtop

        ! Periodicity in x
        u(1-halosize,:) = u(Nx,:)
        u(Nx+halosize,:) = u(1,:)

        w(1-halosize,:) = w(Nx,:)
        w(Nx+halosize,:) = w(1,:)

    end subroutine update_ghost_walls


    subroutine update_ghost_wallTemp(temp,Tbot,Ttop)
        use parameters, only: Nx, Nz
        implicit none
        real(8), intent(inout), dimension(halosize-1:Nx+halosize , halosize-1:Nz+halosize ) :: temp
        real(8) :: Tbot, Ttop


        ! Bottom wall
        temp(:,1-halosize) = 2.0 * Tbot - temp(:,1)

        ! Top wall
        temp(:,Nz+halosize) = 2.0 * Ttop - temp(:,Nz)

        ! Periodicity in x
        temp(1-halosize,:) = temp(Nx,:)
        temp(Nx+halosize,:) = temp(1,:)

    end subroutine update_ghost_wallTemp

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