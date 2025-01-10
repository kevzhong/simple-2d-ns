# simple-2d-ns

Hello World to incompressible Navier--Stokes. The basic features of this solver are:
- Two-dimensional flow
- Second-order finite difference, staggered grid
- Fractional-step method
- Pressure-Poisson equation solved via FFTs
- OpenMP acceleration
- Fully explicit Euler time-integration
- Uniform grid spacing

Various branches have been created in what is (generally) increasing order of implementation complexity:
- `main`: Doubly-periodic domain, Pressure-Poisson equation solved via 2D FFTs
- `walls`: Implementation of wall-boundary conditions in one direction. The Poisson equation is solved with 1D FFTs in the periodic x-direction and a tridiagonal inversion in the remaining wall-normal direction. A scalar field is also added so that 2D Rayleigh-Benard flow can be simulated.
- `implicit`: Implicit time-stepping treatment of the diffusive terms using an Alternating Direct Implicit (ADI) approach or Helmholtz solver. 
- `RK3`: Implementation of RK3 sub-stepping. Wall-normal viscous terms are always treated implicitly, while an option is available to treat the wall-parallel terms implicitly.
- `non-uniform-grid`: Grid-stretching in the wall-normal direction. WIP.

Fields are dumped as binary files. In the later branches, I have routines which dump additional XDMF instructions such that they can be conveniently imported into Paraview.

## Dependencies

- FFTW3
- Fortran compiler
- OpenMP

## Acknowledgements

* [Dr. Naoki Hori](https://naokihori.github.io/NaokiHori/) for his overall guidance in kick-starting this project
* My PhD advisor [A/Prof. Daniel Chung](https://people.eng.unimelb.edu.au/chungd1/), of which I have used his code as a reference for various intricacies related to incompressible flow solvers
* The contributors to [AFiD](https://github.com/chowland/AFiD-MuRPhFi), which I have also used as a reference