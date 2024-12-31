module bctypes
    implicit none
    integer, parameter :: DIRICHLET = 0
    !integer, parameter :: NEUMANN = 1
    integer, parameter :: PRESSUREBC = 2
end module bctypes

module implicit_types
    implicit none
    integer, parameter :: ADI = 0
    integer, parameter :: HELMHOLTZ = 1
end module implicit_types