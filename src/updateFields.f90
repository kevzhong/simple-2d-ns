subroutine updateFields
    use parameters
    implicit none

    call rhsVelocity ! Build RHS vector for velocity solution
    if (scalarmode .eqv. .true.) call rhsScalar ! Build RHS vector for temp solution


    if (implicitmode .eqv. .true.) then    
        if (implicit_type .eq. ADI ) then    
            call ADI_implicitUpdate
        else ! implicit_type .eq. HELMHOLTZ
            call helmholtz_implicitUpdate
        endif
    else
        !call fullyExplicit_update ! Deprecated for rk3 treatment
    endif
end subroutine updateFields