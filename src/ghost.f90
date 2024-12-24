module ghost
    implicit none
    integer :: halosize = 1

    contains

    subroutine update_ghost(array)
        use parameters, only: Nx, Nz
        implicit none
        real(8), intent(inout), dimension(halosize-1:Nx+halosize , halosize-1:Nz+halosize ) :: array

        ! Periodicity application

        ! x-halo swap
        array(1-halosize,:) = array(Nx,:)
        array(Nx+halosize,:) = array(1,:)

        ! z-halo swap
        array(:,halosize-1) = array(:,Nz)
        array(:,Nz+halosize) = array(:,1)

    end subroutine update_ghost

end module ghost