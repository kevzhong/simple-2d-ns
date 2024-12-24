! Allocate memory arrays: fields and working memory arrays
subroutine allocFields
    use parameters, only: Nx, Nz
    use ghost
    use velfields
    use velMemory
    implicit none


    ! Fields
    allocate( u( -halosize+1:Nx+halosize , -halosize+1:Nz+halosize )  )
    allocate( w( -halosize+1:Nx+halosize , -halosize+1:Nz+halosize )  )
    allocate( p( -halosize+1:Nx+halosize , -halosize+1:Nz+halosize )  )

    ! Working memory: RHS vectors
    allocate( rhs_u( 1:Nx , 1:Nz )  )
    allocate( rhs_w( 1:Nx , 1:Nz )  )
    !allocate( adv_u( 1:Nx , 1:Nz )  )
    !allocate( adv_w( 1:Nx , 1:Nz )  )

end subroutine allocFields


subroutine deallocFields
    use velfields
    use velMemory

    if(allocated(u)) deallocate(u)
    if(allocated(w)) deallocate(w)
    if(allocated(p)) deallocate(p)

    if(allocated(rhs_u)) deallocate(rhs_u)
    if(allocated(rhs_w)) deallocate(rhs_w)
    !if(allocated(adv_u)) deallocate(adv_u)
    !if(allocated(adv_w)) deallocate(adv_w)

end subroutine deallocFields
