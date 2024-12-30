! Allocate memory arrays: fields and working memory arrays
subroutine allocFields
    use parameters
    use ghost
    use velfields
    use velMemory
    use scalarfields
    use scalarMemory
    use implicit
    implicit none


    ! Fields
    allocate( u( -halosize+1:Nx+halosize , -halosize+1:Nz+halosize )  )
    allocate( w( -halosize+1:Nx+halosize , -halosize+1:Nz+halosize )  )
    allocate( p( -halosize+1:Nx+halosize , -halosize+1:Nz+halosize )  )
    allocate( pseudo_p( -halosize+1:Nx+halosize , -halosize+1:Nz+halosize )  )

    ! Working memory: RHS vectors
    allocate( rhs_u( 1:Nx , 1:Nz )  )
    allocate( rhs_w( 1:Nx , 1:Nz )  )

    ! Implicit solver working memory
    allocate(ami( 1:Nx )) ; allocate(aci( 1:Nx )) ; allocate(api( 1:Nx ))
    allocate(amk( 1:Nz )) ; allocate(ack( 1:Nz )) ; allocate(apk( 1:Nz ))
    allocate(tdm_rhsX1( 1:Nx )) ; allocate(tdm_rhsX2( 1:Nx )) ; 
    allocate(tdm_rhsZ( 1:Nz ))
    allocate(impl_delta( 1:Nx, 1:Nz ))


 
    !allocate( adv_u( 1:Nx , 1:Nz )  )
    !allocate( adv_w( 1:Nx , 1:Nz )  )

    if (scalarmode .eqv. .true.) then
        allocate( temp( -halosize+1:Nx+halosize , -halosize+1:Nz+halosize )  )
        allocate( rhs_temp( 1:Nx , 1:Nz )  )
    endif

end subroutine allocFields


subroutine deallocFields
    use velfields
    use velMemory
    use scalarfields
    use scalarMemory
    use implicit
    implicit none

    if(allocated(u)) deallocate(u)
    if(allocated(w)) deallocate(w)
    if(allocated(p)) deallocate(p)
    if(allocated(pseudo_p)) deallocate(pseudo_p)

    if(allocated(rhs_u)) deallocate(rhs_u)
    if(allocated(rhs_w)) deallocate(rhs_w)
    !if(allocated(adv_u)) deallocate(adv_u)
    !if(allocated(adv_w)) deallocate(adv_w)

    if(allocated(ami)) deallocate(ami) ; if(allocated(aci)) deallocate(aci) ; if(allocated(api)) deallocate(api)
    if(allocated(amk)) deallocate(amk) ; if(allocated(ack)) deallocate(ack) ; if(allocated(apk)) deallocate(apk)
    if(allocated(tdm_rhsX1)) deallocate(tdm_rhsX1) ; if(allocated(tdm_rhsX2)) deallocate(tdm_rhsX2) 
    if(allocated(tdm_rhsZ)) deallocate(tdm_rhsZ)
    if(allocated(impl_delta)) deallocate(impl_delta)

    if(allocated(temp)) deallocate(temp)
    if(allocated(rhs_temp)) deallocate(rhs_temp)


end subroutine deallocFields
