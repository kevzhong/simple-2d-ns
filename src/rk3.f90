module rk3
    implicit none

    real, dimension(3) :: al_coeff =  [ 8.0/15.0,        2.0/15.0,  1.0/3.0   ]
    real, dimension(3) :: gam_coeff = [ 8.0/15.0,        5.0/12.0,  3.0/4.0   ]
    real, dimension(3) :: zet_coeff = [      0.0,      -17.0/60.0, -5.0/12.0  ]

    real :: aldt, gamdt, zetdt

    contains

    subroutine next_rk3(substep)
        use parameters
        use velMemory
        use scalarMemory
        implicit none
        integer, intent(in) :: substep

        aldt = al_coeff(substep) ; zetdt = zet_coeff(substep) ; gamdt = gam_coeff(substep)
        
        call rk3_memSwap(expl_u, expl_u_m1)
        call rk3_memSwap(expl_w, expl_w_m1)
        
        if (scalarmode .eqv. .true.) call rk3_memSwap(expl_c, expl_c_m1)

    end subroutine next_rk3

    subroutine rk3_memSwap(expl, expl_m1)
        use parameters, only: Nx, Nz
        implicit none
        real, dimension(Nx,Nz), intent(inout) :: expl_m1
        real, dimension(Nx,Nz), intent(in) :: expl
        integer :: i,k

        !$omp parallel do &
        !$omp default(none) &
        !$omp private(i,k) &
        !$omp shared(expl,expl_m1,Nx,Nz)
        do k = 1,Nz
            do i = 1,Nx
                expl_m1(i,k) = expl(i,k)
            enddo
        enddo
        !$omp end parallel do

    end subroutine rk3_memSwap

end module rk3