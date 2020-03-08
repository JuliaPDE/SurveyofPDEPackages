# Survey of PDE Packages

State of the ecosystem as of: 03/05/2020

If you think something was missed, if you’d like to amend or complement the information, or if, for any reason, you wish your software not to be included, file an issue, or even better, make it a PR.

For some projects the actual software is not available (or it is not clear how to get it). In that case I would at least include a link to the paper or other source of information.

## Table of contents

- [General PDE approximation methods](#general)
- [Transform methods](#tm)
- [Finite difference methods](#fdm)
- [Finite element methods](#femfvm)
- [Finite volume methods](#femfvm)
- [Spectral element methods](#sem)
- [Boundary element, Boundary integral methods](#bie)
- [Mesh free methods](#mfe)
- [Virtual element methods](#vem)
- [Multi-method packages](#mm)
- [Non-classical methods](#nonclassical)
- [Solvers, sparse and hierarchical matrix libraries](#solvers)
- [Mesh and Grid Generation](#grids)
- [Postprocessing, visualization](#post)
- [HPC, Parallel processing](#hpc)
- [Miscellanea](#misc)

## <a name="general"></a>General PDE approximation methods

### [https://github.com/JuliaApproximation/ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl)

ApproxFun is a package for approximating functions. It is in a similar vein to the Matlab package Chebfun and the Mathematica package RHPackage. Active and high quality project.

### [https://juliadiffeq.org/](https://juliadiffeq.org/)

JuliaDiffEq is a Github organization created to unify the packages for solving differential equations in Julia. By providing a diverse set of tools with a common interface, we provide a modular, easily-extendable, and highly performant ecosystem for solving various forms of differential equations.

## <a name="tm"></a>Transform methods

### [https://github.com/johnfgibson/julia-pde-benchmark](https://github.com/johnfgibson/julia-pde-benchmark)

Benchmarking a simple PDE integration algorithm in Julia and other languages.
Fourier approach.

### [https://github.com//JuliaMolSim/DFTK.jl/](https://github.com//JuliaMolSim/DFTK.jl/)

The density-functional toolkit is a library of Julia routines for experimentation with plane-wave-based density-functional theory (DFT): it is an engine to solve nonlinear eigenvector equations discretized in a Fourier basis, applied to the Kohn-Sham equations of electronic structure theory (as well as a couple of others).

## <a name="fdm"></a>Finite difference methods

### [https://github.com/JuliaDiffEq/DiffEqOperators.jl](https://github.com/JuliaDiffEq/DiffEqOperators.jl)

Automatic construction of arbitrary order finite difference stencils on regular and irregular grids. Utilizes stencil compilers and matrix-free implementations for low memory high efficiency implementation.

### [https://github.com/ooreilly/sbp.jl](https://github.com/ooreilly/sbp.jl)

Finite difference method, unregistered.

### [https://github.com/ranocha/SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)

A library of classical summation-by-parts (SBP) operators used in finite difference methods to get provably stable semidiscretisations, paying special attention to boundary conditions.

### [https://github.com/matthieugomez/EconPDEs.jl](https://github.com/matthieugomez/EconPDEs.jl)
This package  solves (systems of) nonlinear ODEs/PDEs arising in economic models (i.e. parabolic/elliptic PDEs arising from HJB equations)
The underlying algorithm is based on a combination of upwinding and fully implicit time stepping, using sparse Jacobians.

## <a name="fem"></a>Finite  element methods

### [https://github.com/JuliaDiffEq/FEniCS.jl](https://github.com/JuliaDiffEq/FEniCS.jl)

FEniCS.jl is a wrapper for the FEniCS library for finite element discretizations of PDEs. 

### [https://github.com/KristofferC/JuAFEM.jl](https://github.com/KristofferC/JuAFEM.jl)

A simple finite element toolbox written in Julia. It is actually quite powerful, and it is being actively updated. Some parallels with deal.II might help with the learning curve.

### [https://github.com/JuliaFEM/](https://github.com/JuliaFEM/)

The JuliaFEM project develops open-source software for reliable, scalable, distributed Finite Element Method.

Maintains an overview page with documentation at
[http://www.juliafem.org/](http://www.juliafem.org/).

### NESSie.jl 

Paper: Efficient and intuitive finite element and boundary element methods for nonlocal protein electrostatics in the Julia language [Link](https://www.sciencedirect.com/science/article/pii/S187775031730738X)

### [researchgate.net/.../GaLerKia-a-unified-Julia-implementation-of-mesh-and-meshfree-based-Galerkin-methods](https://www.researchgate.net/project/GaLerKia-a-unified-Julia-implementation-of-mesh-and-meshfree-based-Galerkin-methods)

Provide a modern and unified simulation platform for finite element methods, meshfree Galerkin methods and virtual element methods. The platform is being built over three Julia(*) independent libraries: FEMLia (for finite elements), MEMLia (for meshfree methods) and VEMLia (for virtual elements).

### [https://github.com/pjabardo/Makhno.jl](https://github.com/pjabardo/Makhno.jl)

### [https://github.com/pjabardo/HPFEM.jl](https://github.com/pjabardo/HPFEM.jl)

HP Finite elements in Julia. One-dimensional. Might have been abandoned.

### [https://github.com/gerhardtulzer/EllipticFEM.jl/](https://github.com/gerhardtulzer/EllipticFEM.jl/)

FEM Solver for Elliptic, Parabolic and Hyperbolic PDEs Written in Julia. No update in the past three years.

### [https://github.com/krcools/BEAST.jl](https://github.com/krcools/BEAST.jl)

The package can support finite element discretization.
Also listed as a boundary-element code. Also see [TeaTalk.jl](https://github.com/krcools/TeaTalk.jl).

### [https://github.com/ga96tik/SauterSchwabQuadrature.jl](https://github.com/ga96tik/SauterSchwabQuadrature.jl)


### [https://github.com/ranocha/HyperbolicDiffEq.jl](https://github.com/ranocha/HyperbolicDiffEq.jl) 

Active with updates for Julia 1.3.

### [https://github.com/gridap/Gridap.jl](https://github.com/gridap/Gridap.jl)

Gridap provides a rich set of tools for the grid-based approximation of PDEs, mainly finite element methods, written in the Julia programming language. Some features of the library are:

 - **Discretization mentods:** Continuous and discontinuous Galerkin methods with Lagrangian, Raviart-Thomas, and Nédélec interpolations of arbitrary order and dimension.
 - **Problem types:** Linear, and non-linear, single-field, and multi-physics problems, both volume-coupled and surface-coupled.
 - **Mesh generation:** Built-in Cartesian mesh generator in arbitrary dimensions; interface with GMSH for unstructured grids via the plugin [GridapGmsh](https://github.com/gridap/GridapGmsh.jl).
 - **Linear and non-linear solvers**: Interfaces with Pardiso and PETSc via the plugins [GridapPardiso](https://github.com/gridap/GridapPardiso.jl) and [GridapPETSc](https://github.com/gridap/GridapPETSc.jl).

### [https://github.com/PetrKryslUCSD/FinEtools.jl.git](https://github.com/PetrKryslUCSD/FinEtools.jl.git)

FinEtools is a package for basic operations on finite element meshes. It also supports a number of packages for applications to heat conduction, acoustics, static and dynamic linear and nonlinear stress analysis, vibrations and fluids,  model reduction, and flexible beams.

### https://github.com/rigetti/DiscretePDEs.jl

DiscretePDEs.jl is a package for discretizing partial differential equations using DiscreteExteriorCalculus.jl. 
In addition to functionality for discretizing arbitrary PDEs, DiscretePDEs.jl also has functionality specifically for modeling electromagnetism.

### [https://github.com/OptimalDesignLab/PDESolver.jl](https://github.com/OptimalDesignLab/PDESolver.jl)

PDESolver is a multi-physics solver primarily focused on Computational Fluid Dynamics. 

### https://github.com/ZenanH/juSFEM
Smoothed FEM.
[Paper ](https://www.sciencedirect.com/science/article/pii/S0898122120300523)

### [https://github.com/dpeschka/jPDE](https://github.com/dpeschka/jPDE)

Partial Differential Equations with Julia, with FEM.

### [www.researchgate.net/..._CFD_Julia_A_Learning_Module_Structuring_an_Introductory_Course_on_Computational_Fluid_Dynamics](https://www.researchgate.net/publication/335398490_CFD_Julia_A_Learning_Module_Structuring_an_Introductory_Course_on_Computational_Fluid_Dynamics)

CFD Julia is a programming module developed for senior undergraduate or graduate-level coursework which teaches the foundations of computational fluid dynamics (CFD). The paper explains various concepts related to spatial and temporal discretization, explicit and implicit numerical schemes, multi-step numerical schemes, higher-order shock-capturing numerical methods, and iterative solvers in CFD. 

### [https://github.com/Paulms/HDiscontinuousGalerkin.jl](https://github.com/Paulms/HDiscontinuousGalerkin.jl)

A finite element toolbox, with focus on Hybridizable Discontinuous Galerkin (at the moment it only works in 2D).

### [https://github.com/ABAtanasov/GalerkinSparseGrids.jl](https://github.com/ABAtanasov/GalerkinSparseGrids.jl)

For solving hyperbolic partial differential equations in higher dimensions, where the curse of dimensionality restricts the computational feasibility of discretization of space using regular grid methods. Instead, the employ  sparse grid construction is employed.

## <a name="fvm"></a>Finite  volume methods

### [https://github.com/madsjulia/FiniteVolume.jl](https://github.com/madsjulia/FiniteVolume.jl)

Finite Volume code. Not much information on the site.

### [https://github.com/j-fu/VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl)

Solver for coupled nonlinear partial differential equations based on the Voronoi finite volume method.

### [https://github.com/simulkade/JFVM.jl](https://github.com/simulkade/JFVM.jl)

Finite volume tool for the transport phenomena in chemical and petroleum engineering and similar fields (linear transient advection-diffusion PDE).
Updated for Julia 1.0.

### [https://github.com/climate-machine/Oceananigans.jl](https://github.com/climate-machine/Oceananigans.jl)

Incompressible fluid flow solver written in Julia that can be run in 1-3 dimensions on CPUs and GPUs. It is designed to solve the rotating Boussinesq equations used in non-hydrostatic ocean modeling but can be used to solve for any incompressible flow.

## <a name="sem"></a>Spectral  element methods

### [https://github.com/pjabardo/SpectralElements.jl](https://github.com/pjabardo/SpectralElements.jl)

Not much activity.

### [https://github.com/ranocha/PolynomialBases.jl](https://github.com/ranocha/PolynomialBases.jl)

A library of functions for polynomial bases used in spectral element methods.

## <a name="bie"></a>Boundary element, Boundary integral methods

### [https://github.com/krcools/BEAST.jl](https://github.com/krcools/BEAST.jl)

This package contains common basis functions and assembly routines for the implementation of boundary element methods. Examples are included for the 2D and 3D Helmholtz equations and for the 3D Maxwell equations.

### NESSie.jl 

Paper: Efficient and intuitive finite element and boundary element methods for nonlocal protein electrostatics in the Julia language [Link](https://www.sciencedirect.com/science/article/pii/S187775031730738X) Also listed as a finite element toolkit.

## <a name="mfe">Mesh free methods


### [www.researchgate.net..._Programming_the_material_point_method_in_Julia](https://www.researchgate.net/publication/312610697_Programming_the_material_point_method_in_Julia)

### [researchgate.net/.../GaLerKia-a-unified-Julia-implementation-of-mesh-and-meshfree-based-Galerkin-methods](https://www.researchgate.net/project/GaLerKia-a-unified-Julia-implementation-of-mesh-and-meshfree-based-Galerkin-methods)

See information above for the finite element methods.

## <a name="vem">Virtual element methods

### [https://github.com/Paulms/jFEMTools.jl](https://github.com/Paulms/jFEMTools.jl)

Tools for FEM and VEM (Virtual Element) methods.

### [researchgate.net/.../GaLerKia-a-unified-Julia-implementation-of-mesh-and-meshfree-based-Galerkin-methods](https://www.researchgate.net/project/GaLerKia-a-unified-Julia-implementation-of-mesh-and-meshfree-based-Galerkin-methods)

See information above for the finite element methods.

## <a name="mm">Multi-method packages

### [https://github.com/climate-machine/CLIMA](https://github.com/climate-machine/CLIMA)

The Climate Machine is a new Earth system model that leverages recent advances in the computational and data sciences to learn directly from a wealth of Earth observations from space and the ground.


## <a name="nonclassical">Non-classical methods

### [https://github.com/JuliaDiffEq/NeuralNetDiffEq.jl](https://github.com/JuliaDiffEq/NeuralNetDiffEq.jl)

A library for solving (partial) differential equations with neural networks. Currently supports parabolic differential equations, though a generic NN-based PDE solver is in progress. Can solve very high dimensional (hundred or thousand) partial differential equations through [universal differential equation](https://arxiv.org/abs/2001.04385) approaches.

## <a name="solvers">Solvers, sparse and hierarchical matrix libraries
  
### [https://github.com/JuliaDiffEq/DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl)

A package for solving time-stepping of differential equations which result from PDE discretizations. Heavy emphasis on large-scale stiff differential equations with sparse Jacobians, i.e. DEs from PDEs. See [the stiff ODE tutorial](https://docs.juliadiffeq.org/dev/tutorials/advanced_ode_example/) for more details.

### [https://github.com/JuliaParallel/PETSc.jl](https://github.com/JuliaParallel/PETSc.jl)

This package provides a high level interface for PETSc.

### [https://github.com/OptimalDesignLab/PETSc2.jl](https://github.com/OptimalDesignLab/PETSc2.jl)

This package provides thin wrappers for PETSc, as well as a few convenience functions that take advantage of multiple dispatch.

### [https://github.com/timholy/PositiveFactorizations.jl](https://github.com/timholy/PositiveFactorizations.jl)

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

### [https://github.com/krcools/ClusterTrees.jl](https://github.com/krcools/ClusterTrees.jl)

Tree data structures for fast multipole methods and H-matrices.


### [https://github.com/JuliaLinearAlgebra/SuiteSparse.jl](https://github.com/JuliaLinearAlgebra/SuiteSparse.jl)


SuiteSparse wrappers in Julia.

### [https://github.com/JuliaLinearAlgebra/Arpack.jl](https://github.com/JuliaLinearAlgebra/Arpack.jl)
Julia wrapper for the arpack library designed to solve large scale eigenvalue problems. ARPACK is a collection of Fortran77 subroutines designed to solve large scale eigenvalue problems.

### [https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl)
Solve sparse linear systems using Algebraic Multigrid (AMG). This works especially well for symmetric positive definite matrices.

### [https://github.com/nep-pack/NonlinearEigenproblems.jl](https://github.com/nep-pack/NonlinearEigenproblems.jl)

This package aims to provide state-of-the-art algorithms to solve the nonlinear eigenvalue problem. This currently includes (but is not restricted to) Newton-type methods, Subspace methods, Krylov methods, contour integral methods, block methods, companion matrix approaches. Problem transformation techniques such as scaling, shifting, deflating are also natively supported by the package.

### [https://github.com/rveltz/PseudoArcLengthContinuation.jl](https://github.com/rveltz/PseudoArcLengthContinuation.jl)

This  package aims at solving equations F(u,λ)=0 where λ∈ℝ starting from an initial guess (u0,λ0). It relies on the pseudo arclength continuation algorithm which provides a predictor (u1,λ1) from (u0,λ0). A Newton method is then used to correct this predictor. The focus is on large scale nonlinear problems and multiple hardwares, and the goal is to use Matrix Free methods on a GPU (see PDE example and Periodic orbit example) or on a cluster to solve non linear PDE, nonlocal problems, compute sub-manifolds and so on.

### [https://github.com/Jutho/KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl)
A Julia package collecting a number of Krylov-based algorithms for linear problems, singular value and eigenvalue problems and the application of functions of linear maps or operators to vectors.

## <a name="geo">Geometry and topology

### [https://github.com/JuliaNLSolvers/Manifolds.jl](https://github.com/JuliaNLSolvers/Manifolds.jl)

Manifolds.jl aims to provide both a unified interface to define and use manifolds as well as a library of manifolds to use for your projects.

### [https://github.com/chakravala/Grassmann.jl](https://github.com/chakravala/Grassmann.jl)

The Grassmann.jl package provides tools for doing computations based on multi-linear algebra, differential geometry, and spin groups using the extended tensor algebra known as Leibniz-Grassmann-Clifford-Hestenes geometric algebra.

## <a name="grids">Mesh and Grid Generation

### [https://github.com/JuliaFEM/Gmsh.jl](https://github.com/JuliaFEM/Gmsh.jl)

Gmsh.jl contains API for Gmsh: a three-dimensional finite element mesh generator. With the help of Gmsh.jl, it is possible add parametric model construction and/or automatic mesh generation to a FEM/FVM pipeline.

### [https://github.com/JuliaGeometry/Triangulate.jl](https://github.com/JuliaGeometry/Triangulate.jl)

Julia wrapper for Johnathan Richard Shewchuk's Triangle mesh generator. 

### [https://github.com/JuliaGeometry/TetGen.jl](https://github.com/JuliaGeometry/TetGen.jl)

The TetGen.jl package is a Julia wrapper for the C++ project TetGen. This wrapper enables TetGen based tetrahedral meshing, and (constrained) 3D Delaunay and Voronoi tesselation.

### [https://github.com/krcools/CompScienceMeshes.jl](https://github.com/krcools/CompScienceMeshes.jl)

Geometry types and algorithms for computational science. Meshes, charts, and neighborhoods.

### https://github.com/JuliaGeometry/VoronoiDelaunay.jl

### [https://github.com/krcools/FMMTrees.jl](https://github.com/krcools/FMMTrees.jl)

Tree data structures for H(2), hierarchical matrices, and FMM-like algorithms (fast multiple methods).

### [https://github.com/JuliaGeometry/DistMesh.jl](https://github.com/JuliaGeometry/DistMesh.jl)

Tetrahedral mesh refinement on signed distance/implicit functions or level sets using TetGen.

## <a name="post"></a>Postprocessing, visualization

### [https://github.com/jipolanco/WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl).

This module allows to write VTK XML files, that can be visualised for example with ParaView. Seems pretty complete, writes compressed files.

## <a name="hpc"></a>HPC, Parallel processing

### [https://github.com/JuliaParallel/MPI.jl](https://github.com/JuliaParallel/MPI.jl)

This provides Julia interface to the Message Passing Interface (MPI), roughly inspired by mpi4py.

### [https://github.com/OptimalDesignLab/PumiInterface.jl](https://github.com/OptimalDesignLab/PumiInterface.jl)

This code provides a way to use PUMI from Julia by wrapping functions in PUMIs APF API. 


## <a name="misc"></a>Miscellanea

### [https://github.com/avigliotti/AD4SM.jl](https://github.com/avigliotti/AD4SM.jl)

Automatic Differentiation for Solid Mechanics in Julia.

### [https://github.com/JuliaRheology/RHEOS.jl](https://github.com/JuliaRheology/RHEOS.jl)

RHEOS, an abbreviation of Rheology Open Source, is a software package written in the Julia programming language that provides tools for analyzing rheological data.

### [https://github.com/KristofferC/Tensors.jl](https://github.com/KristofferC/Tensors.jl)

Efficient computations with symmetric and non-symmetric tensors with support for automatic differentiation.

### [https://github.com/Jutho/TensorOperations.jl](https://github.com/Jutho/TensorOperations.jl)

Fast tensor operations using a convenient Einstein index notation.

## Not categorized yet

https://github.com/JerryLingjieMei/NonlinearPDE
https://github.com/DanPSilva/Partial-Differential-Equations

https://github.com/FourierFlows/FourierFlows.jl
https://github.com/PtFEM/PtFEM.jl
https://github.com/samuelpowell/LibTOAST.jl
