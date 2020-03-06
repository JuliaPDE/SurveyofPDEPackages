# Survey of PDE Packages

State of the ecosystem as of: 03/05/2020

## General PDE approximation methods

### [https://github.com/JuliaApproximation/ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl)

ApproxFun is a package for approximating functions. It is in a similar vein to the Matlab package Chebfun and the Mathematica package RHPackage. Active and high quality project.

### [https://juliadiffeq.org/](https://juliadiffeq.org/)

### https://github.com/JuliaDiffEq/FEniCS.jl



## Finite element, finite volume,  spectral element methods

### [https://github.com/KristofferC/JuAFEM.jl](https://github.com/KristofferC/JuAFEM.jl)

A simple finite element toolbox written in Julia. It is actually quite powerful, and it is being actively updated. Some parallels with deal.II might help with the learning curve.

### [https://github.com/JuliaFEM/](https://github.com/JuliaFEM/)

The JuliaFEM project develops open-source software for reliable, scalable, distributed Finite Element Method.

Maintains an overview page with documentation at
[http://www.juliafem.org/](http://www.juliafem.org/).

- NESSie.jl – Efficient and intuitive finite element and boundary element methods for nonlocal protein electrostatics in the Julia language [Link](https://www.sciencedirect.com/science/article/pii/S187775031730738X)
- [https://www.researchgate.net/project/GaLerKia-a-unified-Julia-implementation-of-mesh-and-meshfree-based-Galerkin-methods](https://www.researchgate.net/project/GaLerKia-a-unified-Julia-implementation-of-mesh-and-meshfree-based-Galerkin-methods)
- 
- https://github.com/pjabardo/Makhno.jl
### [https://github.com/pjabardo/HPFEM.jl](https://github.com/pjabardo/HPFEM.jl)

HP Finite elements in Julia. One-dimensional. Might have been abandoned.

### [https://github.com/pjabardo/SpectralElements.jl](https://github.com/pjabardo/SpectralElements.jl)

Not much activity.

### https://github.com/JuliaDiffEq/FEniCS.jl

### https://github.com/simulkade/JFVM.jl

Finite volume tool for the transport phenomena in chemical and petroleum engineering and similar fields (linear transient advection-diffusion PDE).
Updated for Julia 1.0.

### [https://github.com/gerhardtulzer/EllipticFEM.jl/](https://github.com/gerhardtulzer/EllipticFEM.jl/)

FEM Solver for Elliptic, Parabolic and Hyperbolic PDEs Written in Julia. No update in the past three years.

- [https://github.com/krcools/BEAST.jl](https://github.com/krcools/BEAST.jl)
- https://github.com/krcools/ClusterTrees.jl
- [https://github.com/ga96tik/SauterSchwabQuadrature.jl](https://github.com/ga96tik/SauterSchwabQuadrature.jl)


### [https://github.com/ranocha/HyperbolicDiffEq.jl](https://github.com/ranocha/HyperbolicDiffEq.jl) 

Active with updates for Julia 1.3.

### [https://github.com/ranocha/PolynomialBases.jl](https://github.com/ranocha/PolynomialBases.jl)

A library of functions for polynomial bases used in spectral element methods.

- [https://github.com/madsjulia/FiniteVolume.jl](https://github.com/madsjulia/FiniteVolume.jl)
- [https://github.com/Paulms/jFEMTools.jl](https://github.com/Paulms/jFEMTools.jl)
### [https://github.com/gridap/Gridap.jl](https://github.com/gridap/Gridap.jl)

Gridap provides a set of tools for the grid-based approximation of partial differential equations (PDEs) written in the Julia programming language. 

### [https://github.com/j-fu/VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl)

Solver for coupled nonlinear partial differential equations based on the Voronoi finite volume method.

### [https://github.com/PetrKryslUCSD/FinEtools.jl.git](https://github.com/PetrKryslUCSD/FinEtools.jl.git)

### https://github.com/rigetti/DiscretePDEs.jl

DiscretePDEs.jl is a package for discretizing partial differential equations using DiscreteExteriorCalculus.jl. 
In addition to functionality for discretizing arbitrary PDEs, DiscretePDEs.jl also has functionality specifically for modeling electromagnetism.

### [https://github.com/OptimalDesignLab/PDESolver.jl](https://github.com/OptimalDesignLab/PDESolver.jl)

PDESolver is a multi-physics solver primarily focused on Computational Fluid Dynamics. 

## Boundary element, Boundary integral methods

### [https://github.com/krcools/BEAST.jl](https://github.com/krcools/BEAST.jl)

This package contains common basis functions and assembly routines for the implementation of boundary element methods. Examples are included for the 2D and 3D Helmholtz equations and for the 3D Maxwell equations.

## Solvers, sparse and hierarchical matrix libraries

### https://github.com/JuliaParallel/PETSc.jl

This package provides a high level interface for PETSc.

### https://github.com/OptimalDesignLab/PETSc2.jl
This package provides thin wrappers for PETSc, as well as a few convenience functions that take advantage of multiple dispatch.

### https://github.com/ranocha/PositiveFactorizations.jl

### [https://github.com/JuliaDiff/FiniteDiff.jl](https://github.com/JuliaDiff/FiniteDiff.jl)

This package is for calculating derivatives, gradients, Jacobians, Hessians, etc. numerically. 

## Mesh and Grid Generation

### [https://github.com/JuliaFEM/Gmsh.jl](https://github.com/JuliaFEM/Gmsh.jl)

## Postprocessing, visualization

### [https://github.com/jipolanco/WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl).

This module allows to write VTK XML files, that can be visualised for example with ParaView. Seems pretty complete, writes compressed files.



## In the need of sorting

https://github.com/dpeschka/jPDE

https://github.com/ZenanH/juSFEM
https://www.sciencedirect.com/science/article/pii/S0898122120300523

https://www.researchgate.net/publication/335398490_CFD_Julia_A_Learning_Module_Structuring_an_Introductory_Course_on_Computational_Fluid_Dynamics

https://github.com/JuliaParallel/MPI.jl

OptimalDesignLab/PUMI.jl

https://www.researchgate.net/publication/312610697_Programming_the_material_point_method_in_Julia

https://github.com/johnfgibson/julia-pde-benchmark

