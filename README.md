# Survey of PDE Packages

State of the ecosystem as of: 03/05/2020

## Table of contents

- [General PDE approximation methods](#general)
- [Finite element, finite volume,  spectral element methods](#femfvm)
- [Boundary element, Boundary integral methods](#bie)
- [Solvers, sparse and hierarchical matrix libraries](#solvers)
- [Mesh and Grid Generation](#grids)
- [Postprocessing, visualization](#post)

## <a name="general"></a>General PDE approximation methods

### [https://github.com/JuliaApproximation/ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl)

ApproxFun is a package for approximating functions. It is in a similar vein to the Matlab package Chebfun and the Mathematica package RHPackage. Active and high quality project.

### [https://juliadiffeq.org/](https://juliadiffeq.org/)

### https://github.com/JuliaDiffEq/FEniCS.jl

## <a name="femfvm"></a>Finite difference, finite element, finite volume, spectral element methods

### [https://github.com/JuliaDiffEq/DiffEqOperators.jl](https://github.com/JuliaDiffEq/DiffEqOperators.jl)

Automatic construction of arbitrary order finite difference stencils on regular and irregular grids. Utilizes stencil compilers and matrix-free implementations for low memory high efficiency implementation.

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

### [https://github.com/krcools/BEAST.jl](https://github.com/krcools/BEAST.jl)

The package can support finite element discretization.
Also listed as a boundary-element code. Also see [TeaTalk.jl](https://github.com/krcools/TeaTalk.jl).

### [https://github.com/ga96tik/SauterSchwabQuadrature.jl](https://github.com/ga96tik/SauterSchwabQuadrature.jl)


### [https://github.com/ranocha/HyperbolicDiffEq.jl](https://github.com/ranocha/HyperbolicDiffEq.jl) 

Active with updates for Julia 1.3.

### [https://github.com/ranocha/PolynomialBases.jl](https://github.com/ranocha/PolynomialBases.jl)

A library of functions for polynomial bases used in spectral element methods.

### [https://github.com/madsjulia/FiniteVolume.jl](https://github.com/madsjulia/FiniteVolume.jl)

Finite Volume code. Not much information on the site.

### [https://github.com/Paulms/jFEMTools.jl](https://github.com/Paulms/jFEMTools.jl)

Tools for FEM and VEM (Virtual Element) methods.

### [https://github.com/gridap/Gridap.jl](https://github.com/gridap/Gridap.jl)

Gridap provides a rich set of tools for the grid-based approximation of PDEs, mainly finite element methods, written in the Julia programming language. Some features of the library are:

 - **Discretization mentods:** Continuous and discontinuous Galerkin methods with Lagrangian, Raviart-Thomas, and Nédélec interpolations of arbitrary order and dimension.
 - **Problem types:** Linear, and non-linear, single-field, and multi-physics problems, both volume-coupled and surface-coupled.
 - **Mesh generation:** Built-in Cartesian mesh generator in arbitrary dimensions; interface with GMSH for unstructured grids via the plugin [GridapGmsh](https://github.com/gridap/GridapGmsh.jl).
 - **Linear and non-linear solvers**: Interfaces with Pardiso and PETSc via the plugins [GridapPardiso](https://github.com/gridap/GridapPardiso.jl) and [GridapPETSc](https://github.com/gridap/GridapPETSc.jl).

### [https://github.com/j-fu/VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl)

Solver for coupled nonlinear partial differential equations based on the Voronoi finite volume method.

### [https://github.com/PetrKryslUCSD/FinEtools.jl.git](https://github.com/PetrKryslUCSD/FinEtools.jl.git)

FinEtools is a package for basic operations on finite element meshes. It also supports a number of packages for applications to heat conduction, acoustics, static and dynamic linear and nonlinear stress analysis, vibrations and fluids,  model reduction, and flexible beams.

### https://github.com/rigetti/DiscretePDEs.jl

DiscretePDEs.jl is a package for discretizing partial differential equations using DiscreteExteriorCalculus.jl. 
In addition to functionality for discretizing arbitrary PDEs, DiscretePDEs.jl also has functionality specifically for modeling electromagnetism.

### [https://github.com/OptimalDesignLab/PDESolver.jl](https://github.com/OptimalDesignLab/PDESolver.jl)

PDESolver is a multi-physics solver primarily focused on Computational Fluid Dynamics. 

## <a name="bie"></a>Boundary element, Boundary integral methods

### [https://github.com/krcools/BEAST.jl](https://github.com/krcools/BEAST.jl)

This package contains common basis functions and assembly routines for the implementation of boundary element methods. Examples are included for the 2D and 3D Helmholtz equations and for the 3D Maxwell equations.

## <a name="nonclassical">Non-classical methods

### [https://github.com/JuliaDiffEq/NeuralNetDiffEq.jl](https://github.com/JuliaDiffEq/NeuralNetDiffEq.jl)

A library for solving (partial) differential equations with neural networks. Currently supports parabolic differential equations, though a generic NN-based PDE solver is in progress. Can solve very high dimensional (hundred or thousand) partial differential equations through [universal differential equation](https://arxiv.org/abs/2001.04385) approaches.

## <a name="solvers">Solvers, sparse and hierarchical matrix libraries
  
### https://github.com/JuliaDiffEq/DifferentialEquations.jl

A package for solving time-stepping of differential equations which result from PDE discretizations. Heavy emphasis on large-scale stiff differential equations with sparse Jacobians, i.e. DEs from PDEs. See [the stiff ODE tutorial](https://docs.juliadiffeq.org/dev/tutorials/advanced_ode_example/) for more details.

### https://github.com/JuliaParallel/PETSc.jl

This package provides a high level interface for PETSc.

### https://github.com/OptimalDesignLab/PETSc2.jl
This package provides thin wrappers for PETSc, as well as a few convenience functions that take advantage of multiple dispatch.

### https://github.com/ranocha/PositiveFactorizations.jl

PositiveFactorizations is a package for computing a positive definite matrix decomposition (factorization) from an arbitrary symmetric input. The motivating application is optimization.

### [https://github.com/JuliaDiff/FiniteDiff.jl](https://github.com/JuliaDiff/FiniteDiff.jl)

This package is for calculating derivatives, gradients, Jacobians, Hessians, etc. numerically. 

### [https://github.com/JuliaDiff/SparseDiffTools.jl](https://github.com/JuliaDiff/SparseDiffTools.jl)

This package is for exploiting sparsity in Jacobians and Hessians to accelerate computations. Matrix-free Jacobian-vector product and Hessian-vector product operators are provided that are compatible with AbstractMatrix-based libraries like IterativeSolvers.jl for easy and efficient Newton-Krylov implementation. It is possible to perform matrix coloring, and utilize coloring in Jacobian and Hessian construction.

### [https://github.com/JuliaMatrices/HierarchicalMatrices.jl](https://github.com/JuliaMatrices/HierarchicalMatrices.jl)

This package provides a flexible framework for hierarchical matrices in Julia.

### [https://github.com/JuliaMatrices/LowRankApprox.jl](https://github.com/JuliaMatrices/LowRankApprox.jl)

This Julia package provides fast low-rank approximation algorithms for BLAS/LAPACK-compatible matrices based on some of the latest technology in adaptive randomized matrix sketching. 

### [https://github.com/JuliaMatrices/LazyBandedMatrices.jl](https://github.com/JuliaMatrices/LazyBandedMatrices.jl)

This package supports lazy banded and block-banded matrices.

### [https://github.com/JuliaMatrices/BlockBandedMatrices.jl](https://github.com/JuliaMatrices/BlockBandedMatrices.jl)

A Julia package for representing block-block-banded matrices and banded-block-banded matrices.

### [https://bitbucket.org/cgeoga/kernelmatrices.jl](https://bitbucket.org/cgeoga/kernelmatrices.jl)

This software suite is a companion to the manuscript Scalable Gaussian Process Computations using Hierarchical Matrices.

### https://github.com/krcools/ClusterTrees.jl

Tree data structures for fast multipole methods and H-matrices.

## <a name="grids">Mesh and Grid Generation

### [https://github.com/JuliaFEM/Gmsh.jl](https://github.com/JuliaFEM/Gmsh.jl)

Gmsh.jl contains API for Gmsh: a three-dimensional finite element mesh generator. With the help of Gmsh.jl, it is possible add parametric model construction and/or automatic mesh generation to a FEM/FVM pipeline.

### [https://github.com/JuliaGeometry/Triangulate.jl](https://github.com/JuliaGeometry/Triangulate.jl)

Julia wrapper for Johnathan Richard Shewchuk's Triangle mesh generator. 

### [https://github.com/JuliaGeometry/TetGen.jl](https://github.com/JuliaGeometry/TetGen.jl)

The TetGen.jl package is a Julia wrapper for the C++ project TetGen. This wrapper enables TetGen based tetrahedral meshing, and (constrained) 3D Delaunay and Voronoi tesselation.

### [https://github.com/krcools/CompScienceMeshes.jl](https://github.com/krcools/CompScienceMeshes.jl)

Geometry types and algorithms for computational science. Meshes, charts, and neighborhoods.

## <a name="post"></a>Postprocessing, visualization

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

https://github.com/JuliaGeometry/VoronoiDelaunay.jl



https://github.com/nep-pack/NonlinearEigenproblems.jl

https://github.com/JuliaNLSolvers/Manifolds.jl
