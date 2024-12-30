subroutine tridiag(ami,aci,api,rhs,N)
        implicit none
        integer :: i, N
        real, dimension(N) :: ami,aci,api,rhs
        real, dimension(N) :: scratch
        real    :: val
    
        ! Invert a standard tridiagonal system
        ! KZ This is an implementation following Naoki's documentation
        ! https://naokihori.github.io/SimpleNSSolver/numerical_method/tdm.html
        ! See also the wikipedia page on tridiagonal matrix algorithms
    
        ! Solution vector is stored in rhs at end
    
        scratch(1) = api(1) / aci(1)
        rhs(1) = rhs(1) / aci(1)
    
        ! Forward eliminiation step
        do i = 2,(N-1)
            val = 1.0 / ( aci(i) - ami(i) * scratch(i-1) )
            scratch(i) = api(i) * val
            rhs(i) = val * ( rhs(i) - ami(i) * rhs(i-1) )
        enddo
    
    
        val = aci(N) - ami(N) * scratch(N-1)
    
        ! ! Final entry of forward elimination, test to see if matrix is singular
        ! if ( abs(val) .gt. epsilon(1.0d0) ) then
        !     rhs(N) = 1.0 / val * ( rhs(N) - ami(N) * rhs(N-1) )
        ! else
        !     rhs(N) = 0.0
        ! endif

        rhs(N) = 1.0 / val * ( rhs(N) - ami(N) * rhs(N-1) )

    
        ! Backward substitution
        do i = N-1, 1, -1
            rhs(i) = rhs(i) - scratch(i) * rhs(i+1)
        enddo
    
    
        return
   
end subroutine tridiag
    