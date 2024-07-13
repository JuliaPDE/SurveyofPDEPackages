# Survey of PDE Packages

State of the ecosystem as of: 10/02/2023

This is a brief list of packages relevant when solving partial differential equations with Julia. The information is mostly gleaned from repositories of packages or from published reports or articles. If

- you think something was missed, 
- you’d like to amend or complement the information, or 
- you wish your software not to be included, 

file an issue, or even better, make it a PR.

For some projects the actual software is not available (or it is not clear how to get it). In that case the document would at least include a link to the paper or other source of information.

## Table of contents

- [General PDE approximation methods](#general)
- [Transform methods](#tm)
- [Finite difference methods](#fdm)
- [Finite element methods](#fem)
- [Finite volume methods](#fvm)
- [Spectral element methods](#sem)
- [Boundary element, Boundary integral methods](#bie)
- [Mesh free methods and particle methods](#mfe)
- [Virtual element methods](#vem)
- [Multi-method packages](#mm)
- [Non-classical methods](#nonclassical)
- [Solvers, sparse and hierarchical matrix libraries](#solvers)
- [Geometry and topology](#geo)
- [Mesh and Grid Generation](#grids)
- [Postprocessing, visualization](#post)
- [HPC, Parallel processing](#hpc)
- [Optimization](#opt)
- [Miscellanea](#misc)

## <a name="general"></a>General PDE approximation methods

### [ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaApproximation/ApproxFun.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaApproximation/ApproxFun.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaApproximation/ApproxFun.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaApproximation/ApproxFun.jl)
![GitHub license](https://img.shields.io/github/license/JuliaApproximation/ApproxFun.jl)

ApproxFun is a package for approximating functions. It is in a similar vein to the Matlab package Chebfun and the Mathematica package RHPackage. Active and high quality project.

### [DiffEqDocs](https://docs.sciml.ai/DiffEqDocs/stable/)

![GitHub contributors](https://img.shields.io/github/contributors/SciML/DifferentialEquations.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/SciML/DifferentialEquations.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/SciML/DifferentialEquations.jl)
![GitHub stars](https://img.shields.io/github/stars/SciML/DifferentialEquations.jl)
![GitHub license](https://img.shields.io/github/license/SciML/DifferentialEquations.jl)

DifferentialEquations.jl: Efficient Differential Equation Solving in Julia. This is a suite for numerically solving differential equations written in Julia and available for use in Julia, Python, and R. The purpose of this package is to supply efficient Julia implementations of solvers for various differential equations. Equations within the realm of this package include:

- Discrete equations (function maps, discrete stochastic (Gillespie/Markov) simulations)
- Ordinary differential equations (ODEs)
- Split and Partitioned ODEs (Symplectic integrators, IMEX Methods)
- Stochastic ordinary differential equations (SODEs or SDEs)
- Stochastic differential-algebraic equations (SDAEs)
- Random differential equations (RODEs or RDEs)
- Differential algebraic equations (DAEs)
- Delay differential equations (DDEs)
- Neutral, retarded, and algebraic delay differential equations (NDDEs, RDDEs, and DDAEs)
- Stochastic delay differential equations (SDDEs)
- Experimental support for stochastic neutral, retarded, and algebraic delay differential equations (SNDDEs, SRDDEs, and SDDAEs)
- Mixed discrete and continuous equations (Hybrid Equations, Jump Diffusions)
- (Stochastic) partial differential equations ((S)PDEs) (with both finite difference and finite element methods)


## <a name="tm"></a>Transform methods

### [julia-pde-benchmark](https://github.com/johnfgibson/julia-pde-benchmark)

![GitHub contributors](https://img.shields.io/github/contributors/johnfgibson/julia-pde-benchmark)
![GitHub closed issues](https://img.shields.io/github/issues-closed/johnfgibson/julia-pde-benchmark)
![GitHub last commit](https://img.shields.io/github/last-commit/johnfgibson/julia-pde-benchmark)
![GitHub stars](https://img.shields.io/github/stars/johnfgibson/julia-pde-benchmark)
![GitHub license](https://img.shields.io/github/license/johnfgibson/julia-pde-benchmark)

Benchmarking a simple PDE integration algorithm in Julia and other languages.
Fourier approach.

### [DFTK.jl](https://github.com//JuliaMolSim/DFTK.jl/)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaMolSim/DFTK.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaMolSim/DFTK.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaMolSim/DFTK.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaMolSim/DFTK.jl)
![GitHub license](https://img.shields.io/github/license/JuliaMolSim/DFTK.jl)

The density-functional toolkit is a library of Julia routines for experimentation with plane-wave-based density-functional theory (DFT): it is an engine to solve nonlinear eigenvector equations discretized in a Fourier basis, applied to the Kohn-Sham equations of electronic structure theory (as well as a couple of others).

### [FourierFlows.jl](https://github.com/FourierFlows/FourierFlows.jl)

![GitHub contributors](https://img.shields.io/github/contributors/FourierFlows/FourierFlows.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/FourierFlows/FourierFlows.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/FourierFlows/FourierFlows.jl)
![GitHub stars](https://img.shields.io/github/stars/FourierFlows/FourierFlows.jl)
![GitHub license](https://img.shields.io/github/license/FourierFlows/FourierFlows.jl)

Tools for solving partial differential equations on periodic domains using Fourier-based pseudospectral methods.

### [PencilFFTs.jl](https://github.com/jipolanco/PencilFFTs.jl) 

![GitHub contributors](https://img.shields.io/github/contributors/jipolanco/PencilFFTs.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/jipolanco/PencilFFTs.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/jipolanco/PencilFFTs.jl)
![GitHub stars](https://img.shields.io/github/stars/jipolanco/PencilFFTs.jl)
![GitHub license](https://img.shields.io/github/license/jipolanco/PencilFFTs.jl)

Fast Fourier transforms of MPI-distributed Julia arrays that can be used for pseudospectral partial differential equation solvers.


## <a name="fdm"></a>Finite difference methods

### [DiffEqOperators.jl](https://github.com/JuliaDiffEq/DiffEqOperators.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaDiffEq/DiffEqOperators.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaDiffEq/DiffEqOperators.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaDiffEq/DiffEqOperators.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaDiffEq/DiffEqOperators.jl)
![GitHub license](https://img.shields.io/github/license/JuliaDiffEq/DiffEqOperators.jl)

Automatic construction of arbitrary order finite difference stencils on regular and irregular grids. Utilizes stencil compilers and matrix-free implementations for low memory high efficiency implementation.

### [sbp.jl](https://github.com/ooreilly/sbp.jl)

![GitHub contributors](https://img.shields.io/github/contributors/ooreilly/sbp.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/ooreilly/sbp.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/ooreilly/sbp.jl)
![GitHub stars](https://img.shields.io/github/stars/ooreilly/sbp.jl)
![GitHub license](https://img.shields.io/github/license/ooreilly/sbp.jl)

Finite difference method, unregistered.

### [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)

![GitHub contributors](https://img.shields.io/github/contributors/ranocha/SummationByPartsOperators.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/ranocha/SummationByPartsOperators.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/ranocha/SummationByPartsOperators.jl)
![GitHub stars](https://img.shields.io/github/stars/ranocha/SummationByPartsOperators.jl)
![GitHub license](https://img.shields.io/github/license/ranocha/SummationByPartsOperators.jl)

A library of classical summation-by-parts (SBP) operators used in finite difference methods to get provably stable semidiscretisations, paying special attention to boundary conditions.

### [EconPDEs.jl](https://github.com/matthieugomez/EconPDEs.jl)

![GitHub contributors](https://img.shields.io/github/contributors/matthieugomez/EconPDEs.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/matthieugomez/EconPDEs.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/matthieugomez/EconPDEs.jl)
![GitHub stars](https://img.shields.io/github/stars/matthieugomez/EconPDEs.jl)
![GitHub license](https://img.shields.io/github/license/matthieugomez/EconPDEs.jl)

This package  solves (systems of) nonlinear ODEs/PDEs arising in economic models (i.e. parabolic/elliptic PDEs arising from HJB equations)
The underlying algorithm is based on a combination of upwinding and fully implicit time stepping, using sparse Jacobians.

### [Partial-Differential-Equations](https://github.com/DanPSilva/Partial-Differential-Equations)

![GitHub contributors](https://img.shields.io/github/contributors/DanPSilva/Partial-Differential-Equations)
![GitHub closed issues](https://img.shields.io/github/issues-closed/DanPSilva/Partial-Differential-Equations)
![GitHub last commit](https://img.shields.io/github/last-commit/DanPSilva/Partial-Differential-Equations)
![GitHub stars](https://img.shields.io/github/stars/DanPSilva/Partial-Differential-Equations)
![GitHub license](https://img.shields.io/github/license/DanPSilva/Partial-Differential-Equations)

Solving partial differential equations using finite difference methods on Julia.

### [hyperbolic_PDE.jl](https://github.com/alsam/hyperbolic_PDE.jl)

![GitHub contributors](https://img.shields.io/github/contributors/alsam/hyperbolic_PDE.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/alsam/hyperbolic_PDE.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/alsam/hyperbolic_PDE.jl)
![GitHub stars](https://img.shields.io/github/stars/alsam/hyperbolic_PDE.jl)
![GitHub license](https://img.shields.io/github/license/alsam/hyperbolic_PDE.jl)

The examples are given from the book Numerical Solution of Hyperbolic Partial Differential Equations by John A. Trangenstein.

### [JuliaIBPM](https://github.com/JuliaIBPM)

An ecosystem of repositories for solving PDEs with the immersed boundary projection method on staggered Cartesian grids.

### [ParallelStencil.jl](https://github.com/omlins/ParallelStencil.jl) 

![GitHub contributors](https://img.shields.io/github/contributors/omlins/ParallelStencil.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/omlins/ParallelStencil)
![GitHub last commit](https://img.shields.io/github/last-commit/omlins/ParallelStencil)
![GitHub stars](https://img.shields.io/github/stars/omlins/ParallelStencil)
![GitHub license](https://img.shields.io/github/license/omlins/ParallelStencil)

Package for writing high-level code for parallel high-performance stencil computations that can be deployed on both GPUs and CPUs

### [ImplicitGlobalGrid.jl](https://github.com/eth-cscs/ImplicitGlobalGrid.jl) 

![GitHub contributors](https://img.shields.io/github/contributors/eth-cscs/ImplicitGlobalGrid.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/eth-cscs/ImplicitGlobalGrid.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/eth-cscs/ImplicitGlobalGrid.jl)
![GitHub stars](https://img.shields.io/github/stars/eth-cscs/ImplicitGlobalGrid.jl)
![GitHub license](https://img.shields.io/github/license/eth-cscs/ImplicitGlobalGrid.jl)

Almost trivial distributed parallelization of stencil-based GPU and CPU applications on a regular staggered grid

## <a name="fem"></a>Finite  element methods

### [FEniCS.jl](https://github.com/JuliaDiffEq/FEniCS.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaDiffEq/FEniCS.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaDiffEq/FEniCS.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaDiffEq/FEniCS.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaDiffEq/FEniCS.jl)
![GitHub license](https://img.shields.io/github/license/JuliaDiffEq/FEniCS.jl)

FEniCS.jl is a wrapper for the FEniCS library for finite element discretizations of PDEs.

### [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl)

![GitHub contributors](https://img.shields.io/github/contributors/Ferrite-FEM/Ferrite.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/Ferrite-FEM/Ferrite.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/Ferrite-FEM/Ferrite.jl)
![GitHub stars](https://img.shields.io/github/stars/Ferrite-FEM/Ferrite.jl)
![GitHub license](https://img.shields.io/github/license/Ferrite-FEM/Ferrite.jl)

A simple finite element toolbox written in Julia. It is actually quite powerful, and it is being actively updated. Some parallels with deal.II might help with the learning curve.

### [JuliaFEM](https://github.com/JuliaFEM/)

The JuliaFEM project develops open-source software for reliable, scalable, distributed Finite Element Method.

Maintains an overview page with documentation at
[http://www.juliafem.org/](http://www.juliafem.org/).

### NESSie.jl

Paper: Efficient and intuitive finite element and boundary element methods for nonlocal protein electrostatics in the Julia language [Link](https://www.sciencedirect.com/science/article/pii/S187775031730738X)

### [GaLerKia-a-unified-Julia-implementation-of-mesh-and-meshfree-based-Galerkin-methods](https://www.researchgate.net/project/GaLerKia-a-unified-Julia-implementation-of-mesh-and-meshfree-based-Galerkin-methods)

Provide a modern and unified simulation platform for finite element methods, meshfree Galerkin methods and virtual element methods. The platform is being built over three Julia independent libraries: FEMLia (for finite elements), MEMLia (for meshfree methods) and VEMLia (for virtual elements).

### [Makhno.jl](https://github.com/pjabardo/Makhno.jl)

![GitHub contributors](https://img.shields.io/github/contributors/pjabardo/Makhno.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/pjabardo/Makhno.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/pjabardo/Makhno.jl)
![GitHub stars](https://img.shields.io/github/stars/pjabardo/Makhno.jl)
![GitHub license](https://img.shields.io/github/license/pjabardo/Makhno.jl)

### [HPFEM.jl](https://github.com/pjabardo/HPFEM.jl)

![GitHub contributors](https://img.shields.io/github/contributors/pjabardo/HPFEM.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/pjabardo/HPFEM.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/pjabardo/HPFEM.jl)
![GitHub stars](https://img.shields.io/github/stars/pjabardo/HPFEM.jl)
![GitHub license](https://img.shields.io/github/license/pjabardo/HPFEM.jl)

HP Finite elements in Julia. One-dimensional. Might have been abandoned.

### [EllipticFEM.jl](https://github.com/gerhardtulzer/EllipticFEM.jl/)

![GitHub contributors](https://img.shields.io/github/contributors/gerhardtulzer/EllipticFEM.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/gerhardtulzer/EllipticFEM.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/gerhardtulzer/EllipticFEM.jl)
![GitHub stars](https://img.shields.io/github/stars/gerhardtulzer/EllipticFEM.jl)
![GitHub license](https://img.shields.io/github/license/gerhardtulzer/EllipticFEM.jl)

FEM Solver for Elliptic, Parabolic and Hyperbolic PDEs Written in Julia. No update in the past three years.

### [BEAST.jl](https://github.com/krcools/BEAST.jl)

![GitHub contributors](https://img.shields.io/github/contributors/krcools/BEAST.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/krcools/BEAST.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/krcools/BEAST.jl)
![GitHub stars](https://img.shields.io/github/stars/krcools/BEAST.jl)
![GitHub license](https://img.shields.io/github/license/krcools/BEAST.jl)

The package can support finite element discretization.
Also listed as a boundary-element code. Also see [TeaTalk.jl](https://github.com/krcools/TeaTalk.jl).

### [SauterSchwabQuadrature.jl](https://github.com/ga96tik/SauterSchwabQuadrature.jl)

![GitHub contributors](https://img.shields.io/github/contributors/ga96tik/SauterSchwabQuadrature.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/ga96tik/SauterSchwabQuadrature.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/ga96tik/SauterSchwabQuadrature.jl)
![GitHub stars](https://img.shields.io/github/stars/ga96tik/SauterSchwabQuadrature.jl)
![GitHub license](https://img.shields.io/github/license/ga96tik/SauterSchwabQuadrature.jl)

### [HyperbolicDiffEq.jl](https://github.com/ranocha/HyperbolicDiffEq.jl)

![GitHub contributors](https://img.shields.io/github/contributors/ranocha/HyperbolicDiffEq.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/ranocha/HyperbolicDiffEq.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/ranocha/HyperbolicDiffEq.jl)
![GitHub stars](https://img.shields.io/github/stars/ranocha/HyperbolicDiffEq.jl)
![GitHub license](https://img.shields.io/github/license/ranocha/HyperbolicDiffEq.jl)

Active with updates for Julia 1.3.

### [Gridap.jl](https://github.com/gridap/Gridap.jl)

![GitHub contributors](https://img.shields.io/github/contributors/gridap/Gridap.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/gridap/Gridap.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/gridap/Gridap.jl)
![GitHub stars](https://img.shields.io/github/stars/gridap/Gridap.jl)
![GitHub license](https://img.shields.io/github/license/gridap/Gridap.jl)

Gridap provides a rich set of tools for the grid-based approximation of PDEs, mainly finite element methods, written in the Julia programming language. Some features of the library are:

 - **Discretization mentods:** Continuous and discontinuous Galerkin methods with Lagrangian, Raviart-Thomas, and Nédélec interpolations of arbitrary order and dimension.
 - **Problem types:** Linear, and non-linear, single-field, and multi-physics problems, both volume-coupled and surface-coupled.
 - **Mesh generation:** Built-in Cartesian mesh generator in arbitrary dimensions; interface with GMSH for unstructured grids via the plugin [GridapGmsh](https://github.com/gridap/GridapGmsh.jl).
 - **Linear and non-linear solvers**: Interfaces with Pardiso and PETSc via the plugins [GridapPardiso](https://github.com/gridap/GridapPardiso.jl) and [GridapPETSc](https://github.com/gridap/GridapPETSc.jl).

### [FinEtools.jl](https://github.com/PetrKryslUCSD/FinEtools.jl)

![GitHub contributors](https://img.shields.io/github/contributors/PetrKryslUCSD/FinEtools.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/PetrKryslUCSD/FinEtools.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/PetrKryslUCSD/FinEtools.jl)
![GitHub stars](https://img.shields.io/github/stars/PetrKryslUCSD/FinEtools.jl)
![GitHub license](https://img.shields.io/github/license/PetrKryslUCSD/FinEtools.jl)

FinEtools is a package for basic operations on finite element meshes. It also supports a number of packages for applications to heat conduction, acoustics, static and dynamic linear and nonlinear stress analysis, vibrations and fluids,  model reduction, and flexible beams.

### [DiscretePDEs.jl](https://github.com/rigetti/DiscretePDEs.jl)

![GitHub contributors](https://img.shields.io/github/contributors/rigetti/DiscretePDEs.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/rigetti/DiscretePDEs.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/rigetti/DiscretePDEs.jl)
![GitHub stars](https://img.shields.io/github/stars/rigetti/DiscretePDEs.jl)
![GitHub license](https://img.shields.io/github/license/rigetti/DiscretePDEs.jl)

DiscretePDEs.jl is a package for discretizing partial differential equations using DiscreteExteriorCalculus.jl.
In addition to functionality for discretizing arbitrary PDEs, DiscretePDEs.jl also has functionality specifically for modeling electromagnetism.

### [PDESolver.jl](https://github.com/OptimalDesignLab/PDESolver.jl)

![GitHub contributors](https://img.shields.io/github/contributors/OptimalDesignLab/PDESolver.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/OptimalDesignLab/PDESolver.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/OptimalDesignLab/PDESolver.jl)
![GitHub stars](https://img.shields.io/github/stars/OptimalDesignLab/PDESolver.jl)
![GitHub license](https://img.shields.io/github/license/OptimalDesignLab/PDESolver.jl)

PDESolver is a multi-physics solver primarily focused on Computational Fluid Dynamics.

### [juSFEM](https://github.com/ZenanH/juSFEM)

![GitHub contributors](https://img.shields.io/github/contributors/ZenanH/juSFEM)
![GitHub closed issues](https://img.shields.io/github/issues-closed/ZenanH/juSFEM)
![GitHub last commit](https://img.shields.io/github/last-commit/ZenanH/juSFEM)
![GitHub stars](https://img.shields.io/github/stars/ZenanH/juSFEM)
![GitHub license](https://img.shields.io/github/license/ZenanH/juSFEM)

Smoothed FEM.
[Paper](https://www.sciencedirect.com/science/article/pii/S0898122120300523)

### [jPDE](https://github.com/dpeschka/jPDE)

![GitHub contributors](https://img.shields.io/github/contributors/dpeschka/jPDE)
![GitHub closed issues](https://img.shields.io/github/issues-closed/dpeschka/jPDE)
![GitHub last commit](https://img.shields.io/github/last-commit/dpeschka/jPDE)
![GitHub stars](https://img.shields.io/github/stars/dpeschka/jPDE)
![GitHub license](https://img.shields.io/github/license/dpeschka/jPDE)

Partial Differential Equations with Julia, with FEM.

### [CFD_Julia_A_Learning_Module_Structuring_an_Introductory_Course_on_Computational_Fluid_Dynamics](https://www.researchgate.net/publication/335398490_CFD_Julia_A_Learning_Module_Structuring_an_Introductory_Course_on_Computational_Fluid_Dynamics)

CFD Julia is a programming module developed for senior undergraduate or graduate-level coursework which teaches the foundations of computational fluid dynamics (CFD). The paper explains various concepts related to spatial and temporal discretization, explicit and implicit numerical schemes, multi-step numerical schemes, higher-order shock-capturing numerical methods, and iterative solvers in CFD.

### [HDiscontinuousGalerkin.jl](https://github.com/Paulms/HDiscontinuousGalerkin.jl)

![GitHub contributors](https://img.shields.io/github/contributors/Paulms/HDiscontinuousGalerkin.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/Paulms/HDiscontinuousGalerkin.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/Paulms/HDiscontinuousGalerkin.jl)
![GitHub stars](https://img.shields.io/github/stars/Paulms/HDiscontinuousGalerkin.jl)
![GitHub license](https://img.shields.io/github/license/Paulms/HDiscontinuousGalerkin.jl)

A finite element toolbox, with focus on Hybridizable Discontinuous Galerkin (at the moment it only works in 2D).

### [GalerkinSparseGrids.jl](https://github.com/ABAtanasov/GalerkinSparseGrids.jl)

![GitHub contributors](https://img.shields.io/github/contributors/ABAtanasov/GalerkinSparseGrids.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/ABAtanasov/GalerkinSparseGrids.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/ABAtanasov/GalerkinSparseGrids.jl)
![GitHub stars](https://img.shields.io/github/stars/ABAtanasov/GalerkinSparseGrids.jl)
![GitHub license](https://img.shields.io/github/license/ABAtanasov/GalerkinSparseGrids.jl)

For solving hyperbolic partial differential equations in higher dimensions, where the curse of dimensionality restricts the computational feasibility of discretization of space using regular grid methods. Instead, the employ  sparse grid construction is employed.

### [LibTOAST.jl](https://github.com/samuelpowell/LibTOAST.jl)

![GitHub contributors](https://img.shields.io/github/contributors/samuelpowell/LibTOAST.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/samuelpowell/LibTOAST.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/samuelpowell/LibTOAST.jl)
![GitHub stars](https://img.shields.io/github/stars/samuelpowell/LibTOAST.jl)
![GitHub license](https://img.shields.io/github/license/samuelpowell/LibTOAST.jl)

LibTOAST.jl is a low-level interface to the TOAST++ library. (TOAST++ is an end-to-end solution for forward modelling and image reconstruction in Diffuse Optical Tomography).

### [PtFEM.jl](https://github.com/PtFEM/PtFEM.jl)

![GitHub contributors](https://img.shields.io/github/contributors/PtFEM/PtFEM.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/PtFEM/PtFEM.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/PtFEM/PtFEM.jl)
![GitHub stars](https://img.shields.io/github/stars/PtFEM/PtFEM.jl)
![GitHub license](https://img.shields.io/github/license/PtFEM/PtFEM.jl)

Programs in chapters 4, 5 and early sections of 6 as described in "Programming the Finite Element Method" by Smith, Griffiths and Margetts.

### [NodesAndModes.jl](https://github.com/jlchan/NodesAndModes.jl)

![GitHub contributors](https://img.shields.io/github/contributors/jlchan/NodesAndModes.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/jlchan/NodesAndModes.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/jlchan/NodesAndModes.jl)
![GitHub stars](https://img.shields.io/github/stars/jlchan/NodesAndModes.jl)
![GitHub license](https://img.shields.io/github/license/jlchan/NodesAndModes.jl)

A collection of routines to compute high order "modal" and "nodal" basis functions and their derivatives on 1D, 2D (triangle, quadrilateral), and 3D (hexahedral, tetrahedral, wedge, pyramid) elements. Routines are also provided for quadrature rules and for computing optimized interpolation points. 

### [StartUpDG.jl](https://github.com/jlchan/StartUpDG.jl)

![GitHub contributors](https://img.shields.io/github/contributors/jlchan/StartUpDG.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/jlchan/StartUpDG.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/jlchan/StartUpDG.jl)
![GitHub stars](https://img.shields.io/github/stars/jlchan/StartUpDG.jl)
![GitHub license](https://img.shields.io/github/license/jlchan/StartUpDG.jl)

A translation of the core Matlab codes from [Nodal Discontinuous Galerkin methods](https://www.springer.com/gp/book/9780387720654) to Julia.

### [eFEMpart](https://github.com/pseastham/eFEMpart)

![GitHub contributors](https://img.shields.io/github/contributors/pseastham/eFEMpart)
![GitHub closed issues](https://img.shields.io/github/issues-closed/pseastham/eFEMpart)
![GitHub last commit](https://img.shields.io/github/last-commit/pseastham/eFEMpart)
![GitHub stars](https://img.shields.io/github/stars/pseastham/eFEMpart)
![GitHub license](https://img.shields.io/github/license/pseastham/eFEMpart)

Finite Element code in the Julia language focused on complex fluid-dynamic and porous-media applications, with the possibility of including a particle simulator in the framework of the discrete element method. The 'eFEM' component allows the use of Finite Elements discretizations to solve common problems in fluid dynamics, and the 'part' refers to mesh-free particle methods (discrete element method) primarily aimed at granular-media simulations where continuum constitutive laws are unavailable.

### [Elfel.jl](https://github.com/PetrKryslUCSD/Elfel.jl)

![GitHub contributors](https://img.shields.io/github/contributors/PetrKryslUCSD/Elfel.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/PetrKryslUCSD/Elfel.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/PetrKryslUCSD/Elfel.jl)
![GitHub stars](https://img.shields.io/github/stars/PetrKryslUCSD/Elfel.jl)
![GitHub license](https://img.shields.io/github/license/PetrKryslUCSD/Elfel.jl)

`Elfel.jl` is an Extensible library for Finite Element methods. It  provides support for the development of Finite Element Method applications, especially in continuum mechanics. Mixed methods with cooperating finite element spaces are supported. High performance is one of the focus points.

### [NNFEM.jl](https://github.com/kailaix/NNFEM.jl)

![GitHub contributors](https://img.shields.io/github/contributors/kailaix/NNFEM.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/kailaix/NNFEM.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/kailaix/NNFEM.jl)
![GitHub stars](https://img.shields.io/github/stars/kailaix/NNFEM.jl)
![GitHub license](https://img.shields.io/github/license/kailaix/NNFEM.jl)

NNFEM is a lightweight educational 2D finite element library with truss and 2D quadrilateral elements. Different constitutive relations are supported, including plane stress/strain, hyperelasticity, elasto-plasticity, etc. It supports unstructured grid. It  supports learning a neural network-based constitutive relations with both direct data (i.e, strain-stress pairs) and indirect data (i.e. full displacement field) via automatic differentiation, and solving finite element problems with network-based constitutive relations.

### [Trixi.jl](https://github.com/trixi-framework/Trixi.jl)

![GitHub contributors](https://img.shields.io/github/contributors/trixi-framework/Trixi.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/trixi-framework/Trixi.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/trixi-framework/Trixi.jl)
![GitHub stars](https://img.shields.io/github/stars/trixi-framework/Trixi.jl)
![GitHub license](https://img.shields.io/github/license/trixi-framework/Trixi.jl)

**Trixi.jl** is a numerical simulation framework for hyperbolic conservation
laws written in [Julia](https://julialang.org). A key objective for the
framework is to be useful to both scientists and students. Therefore, next to
having an extensible design with a fast implementation, Trixi is
focused on being easy to use for new or inexperienced users, including the
installation and postprocessing procedures. Its features include:

* Hierarchical quadtree/octree grid with adaptive mesh refinement
* Native support for 1D, 2D, and 3D simulations
* High-order accuracy in space in time
* Nodal discontinuous Galerkin spectral element methods
  * Kinetic energy-preserving and entropy-stable split forms
  * Entropy-stable shock capturing
  * Positivity-preserving limiting
* Compatible with the [SciML ecosystem for ordinary differential equations](https://diffeq.sciml.ai/latest/)
  * [Explicit low-storage Runge-Kutta time integration](https://diffeq.sciml.ai/latest/solvers/ode_solve/#Low-Storage-Methods)
  * [Strong stability preserving methods](https://diffeq.sciml.ai/latest/solvers/ode_solve/#Explicit-Strong-Stability-Preserving-Runge-Kutta-Methods-for-Hyperbolic-PDEs-(Conservation-Laws))
  * CFL-based and error-based time step control
* Square/cubic domains with periodic and weakly-enforced boundary conditions
* Multiple governing equations:
  * Compressible Euler equations
  * Magnetohydrodynamics equations
  * Hyperbolic diffusion equations for elliptic problems
  * Scalar advection
* Multi-physics simulations
  * [Self-gravitating gas dynamics](https://github.com/trixi-framework/paper-self-gravitating-gas-dynamics)
* Shared-memory parallelization via multithreading
* Visualization of results with Julia-only tools ([Trixi2Img](https://github.com/trixi-framework/Trixi2Img.jl))
  or ParaView/VisIt ([Trixi2Vtk](https://github.com/trixi-framework/Trixi2Vtk.jl))


## <a name="fvm"></a>Finite  volume methods

### [FiniteVolume.jl](https://github.com/madsjulia/FiniteVolume.jl)

![GitHub contributors](https://img.shields.io/github/contributors/madsjulia/FiniteVolume.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/madsjulia/FiniteVolume.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/madsjulia/FiniteVolume.jl)
![GitHub stars](https://img.shields.io/github/stars/madsjulia/FiniteVolume.jl)
![GitHub license](https://img.shields.io/github/license/madsjulia/FiniteVolume.jl)

Finite Volume code. Not much information on the site.

### [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl)

![GitHub contributors](https://img.shields.io/github/contributors/j-fu/VoronoiFVM.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/j-fu/VoronoiFVM.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/j-fu/VoronoiFVM.jl)
![GitHub stars](https://img.shields.io/github/stars/j-fu/VoronoiFVM.jl)
![GitHub license](https://img.shields.io/github/license/j-fu/VoronoiFVM.jl)

Solver for coupled nonlinear partial differential equations based on the Voronoi finite volume method.

### [JFVM.jl](https://github.com/simulkade/JFVM.jl)

![GitHub contributors](https://img.shields.io/github/contributors/simulkade/JFVM.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/simulkade/JFVM.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/simulkade/JFVM.jl)
![GitHub stars](https://img.shields.io/github/stars/simulkade/JFVM.jl)
![GitHub license](https://img.shields.io/github/license/simulkade/JFVM.jl)

Finite volume tool for the transport phenomena in chemical and petroleum engineering and similar fields (linear transient advection-diffusion PDE).
Updated for Julia 1.0.

### [Oceananigans.jl](https://github.com/climate-machine/Oceananigans.jl)

![GitHub contributors](https://img.shields.io/github/contributors/climate-machine/Oceananigans.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/climate-machine/Oceananigans.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/climate-machine/Oceananigans.jl)
![GitHub stars](https://img.shields.io/github/stars/climate-machine/Oceananigans.jl)
![GitHub license](https://img.shields.io/github/license/climate-machine/Oceananigans.jl)

Incompressible fluid flow solver written in Julia that can be run in 1-3 dimensions on CPUs and GPUs. It is designed to solve the rotating Boussinesq equations used in non-hydrostatic ocean modeling but can be used to solve for any incompressible flow.

### [Kinetic.jl](https://github.com/vavrines/Kinetic.jl)

![GitHub contributors](https://img.shields.io/github/contributors/vavrines/Kinetic.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/vavrines/Kinetic.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/vavrines/Kinetic.jl)
![GitHub stars](https://img.shields.io/github/stars/vavrines/Kinetic.jl)
![GitHub license](https://img.shields.io/github/license/vavrines/Kinetic.jl)

A lightweight finite volume toolbox for solving Boltzmann equation and its moments systems (Euler, Navier-Stokes, etc).

### [FiniteVolumeMethod.jl](https://github.com/DanielVandH/FiniteVolumeMethod.jl)

![GitHub contributors](https://img.shields.io/github/contributors/DanielVandH/FiniteVolumeMethod.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/DanielVandH/FiniteVolumeMethod.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/DanielVandH/FiniteVolumeMethod.jl)
![GitHub stars](https://img.shields.io/github/stars/DanielVandH/FiniteVolumeMethod.jl)
![GitHub license](https://img.shields.io/github/license/DanielVandH/FiniteVolumeMethod.jl)

This is a package for solving two-dimensional PDEs of type $\partial_t u + \boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{q}(x, y, t, u) = R(x, y, t, u)$ with the finite volume method with a user-provided triangular mesh. Builds on the DifferentialEquations.jl interface to allow for different algorithms, linear solvers, etc. to be used with ease.

### [HighVoronoi.jl](https://github.com/martinheida/HighVoronoi.jl)

![GitHub contributors](https://img.shields.io/github/contributors/martinheida/HighVoronoi.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/martinheida/HighVoronoi.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/martinheida/HighVoronoi.jl)
![GitHub stars](https://img.shields.io/github/stars/martinheida/HighVoronoi.jl)
![GitHub license](https://img.shields.io/github/license/martinheida/HighVoronoi.jl)

A Julia Package for setting up high dimensional (i.e. any dimension >= 2) Finite Volume problems on Voronoi Meshes.

## <a name="sem"></a>Spectral  element methods

### [SpectralElements.jl](https://github.com/pjabardo/SpectralElements.jl)

![GitHub contributors](https://img.shields.io/github/contributors/pjabardo/SpectralElements.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/pjabardo/SpectralElements.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/pjabardo/SpectralElements.jl)
![GitHub stars](https://img.shields.io/github/stars/pjabardo/SpectralElements.jl)
![GitHub license](https://img.shields.io/github/license/pjabardo/SpectralElements.jl)

Not much activity.

### [PolynomialBases.jl](https://github.com/ranocha/PolynomialBases.jl)

![GitHub contributors](https://img.shields.io/github/contributors/ranocha/PolynomialBases.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/ranocha/PolynomialBases.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/ranocha/PolynomialBases.jl)
![GitHub stars](https://img.shields.io/github/stars/ranocha/PolynomialBases.jl)
![GitHub license](https://img.shields.io/github/license/ranocha/PolynomialBases.jl)

### [CFD_Julia_MAE5093](https://github.com/omersan/CFD_Julia_MAE5093)

![GitHub contributors](https://img.shields.io/github/contributors/omersan/CFD_Julia_MAE5093)
![GitHub closed issues](https://img.shields.io/github/issues-closed/omersan/CFD_Julia_MAE5093)
![GitHub last commit](https://img.shields.io/github/last-commit/omersan/CFD_Julia_MAE5093)
![GitHub stars](https://img.shields.io/github/stars/omersan/CFD_Julia_MAE5093)
![GitHub license](https://img.shields.io/github/license/omersan/CFD_Julia_MAE5093)

This repository contains fundamental codes related to CFD that can be included in any graduate level CFD coursework. A number of codes are included in CFD_Julia module, for instance 1D inviscid Burgers equation, 1D Euler equation, 2D Poisson equation, and 2D incompressible Navier-Stokes equations.

A library of functions for polynomial bases used in spectral element methods.

## <a name="bie"></a>Boundary element, boundary integral methods

### [BEAST.jl](https://github.com/krcools/BEAST.jl)

![GitHub contributors](https://img.shields.io/github/contributors/krcools/BEAST.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/krcools/BEAST.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/krcools/BEAST.jl)
![GitHub stars](https://img.shields.io/github/stars/krcools/BEAST.jl)
![GitHub license](https://img.shields.io/github/license/krcools/BEAST.jl)

This package contains common basis functions and assembly routines for the implementation of boundary element methods. Examples are included for the 2D and 3D Helmholtz equations and for the 3D Maxwell equations.


### NESSie.jl

Paper: Efficient and intuitive finite element and boundary element methods for nonlocal protein electrostatics in the Julia language [Link](https://www.sciencedirect.com/science/article/pii/S187775031730738X) Also listed as a finite element toolkit.

## <a name="mfe">Mesh free methods and particle methods

### [Peridynamics.jl](https://github.com/kaipartmann/Peridynamics.jl)

![GitHub contributors](https://img.shields.io/github/contributors/kaipartmann/Peridynamics.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/kaipartmann/Peridynamics.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/kaipartmann/Peridynamics.jl)
![GitHub stars](https://img.shields.io/github/stars/kaipartmann/Peridynamics.jl)
![GitHub license](https://img.shields.io/github/license/kaipartmann/Peridynamics.jl)

Peridynamics is a non-local formulation of continuum mechanics in which material points represent the continuum. It is particularly well suited for dynamic fracture simulations.

### [Programming_the_material_point_method_in_Julia](https://www.researchgate.net/publication/312610697_Programming_the_material_point_method_in_Julia)

### [GaLerKia-a-unified-Julia-implementation-of-mesh-and-meshfree-based-Galerkin-methods](https://www.researchgate.net/project/GaLerKia-a-unified-Julia-implementation-of-mesh-and-meshfree-based-Galerkin-methods)

See information above for the finite element methods.

## <a name="vem">Virtual element methods

### [jFEMTools.jl](https://github.com/Paulms/jFEMTools.jl)

![GitHub contributors](https://img.shields.io/github/contributors/Paulms/jFEMTools.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/Paulms/jFEMTools.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/Paulms/jFEMTools.jl)
![GitHub stars](https://img.shields.io/github/stars/Paulms/jFEMTools.jl)
![GitHub license](https://img.shields.io/github/license/Paulms/jFEMTools.jl)

Tools for FEM and VEM (Virtual Element) methods.

### [GaLerKia-a-unified-Julia-implementation-of-mesh-and-meshfree-based-Galerkin-methods](https://www.researchgate.net/project/GaLerKia-a-unified-Julia-implementation-of-mesh-and-meshfree-based-Galerkin-methods)

See information above for the finite element methods.

### [eFEMpart](https://github.com/pseastham/eFEMpart)

![GitHub contributors](https://img.shields.io/github/contributors/pseastham/eFEMpart)
![GitHub closed issues](https://img.shields.io/github/issues-closed/pseastham/eFEMpart)
![GitHub last commit](https://img.shields.io/github/last-commit/pseastham/eFEMpart)
![GitHub stars](https://img.shields.io/github/stars/pseastham/eFEMpart)
![GitHub license](https://img.shields.io/github/license/pseastham/eFEMpart)

Listed above under finite elements as well. Finite Element code in the Julia language focused on complex fluid-dynamic and porous-media applications, with the possibility of including a particle simulator in the framework of the discrete element method. The 'eFEM' component allows the use of Finite Elements discretizations to solve common problems in fluid dynamics, and the 'part' refers to mesh-free particle methods (discrete element method) primarily aimed at granular-media simulations where continuum constitutive laws are unavailable.

## <a name="mm">Multi-method packages

### [CLIMA](https://github.com/climate-machine/CLIMA)

![GitHub contributors](https://img.shields.io/github/contributors/climate-machine/CLIMA)
![GitHub closed issues](https://img.shields.io/github/issues-closed/climate-machine/CLIMA)
![GitHub last commit](https://img.shields.io/github/last-commit/climate-machine/CLIMA)
![GitHub stars](https://img.shields.io/github/stars/climate-machine/CLIMA)
![GitHub license](https://img.shields.io/github/license/climate-machine/CLIMA)

The Climate Machine is a new Earth system model that leverages recent advances in the computational and data sciences to learn directly from a wealth of Earth observations from space and the ground.

## <a name="nonclassical">Non-classical methods

### [NeuralNetDiffEq.jl](https://github.com/JuliaDiffEq/NeuralNetDiffEq.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaDiffEq/NeuralNetDiffEq.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaDiffEq/NeuralNetDiffEq.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaDiffEq/NeuralNetDiffEq.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaDiffEq/NeuralNetDiffEq.jl)
![GitHub license](https://img.shields.io/github/license/JuliaDiffEq/NeuralNetDiffEq.jl)

A library for solving (partial) differential equations with neural networks. Currently supports parabolic differential equations, though a generic NN-based PDE solver is in progress. Can solve very high dimensional (hundred or thousand) partial differential equations through [universal differential equation](https://arxiv.org/abs/2001.04385) approaches.

## <a name="solvers">Solvers, sparse and hierarchical matrix libraries

### [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaDiffEq/DifferentialEquations.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaDiffEq/DifferentialEquations.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaDiffEq/DifferentialEquations.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaDiffEq/DifferentialEquations.jl)
![GitHub license](https://img.shields.io/github/license/JuliaDiffEq/DifferentialEquations.jl)

A package for solving time-stepping of differential equations which result from PDE discretizations. Heavy emphasis on large-scale stiff differential equations with sparse Jacobians, i.e. DEs from PDEs. See [the stiff ODE tutorial](https://docs.juliadiffeq.org/dev/tutorials/advanced_ode_example/) for more details.

### [PETSc.jl](https://github.com/JuliaParallel/PETSc.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaParallel/PETSc.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaParallel/PETSc.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaParallel/PETSc.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaParallel/PETSc.jl)
![GitHub license](https://img.shields.io/github/license/JuliaParallel/PETSc.jl)

This package provides a high level interface for PETSc.

### [PETSc2.jl](https://github.com/OptimalDesignLab/PETSc2.jl)

![GitHub contributors](https://img.shields.io/github/contributors/OptimalDesignLab/PETSc2.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/OptimalDesignLab/PETSc2.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/OptimalDesignLab/PETSc2.jl)
![GitHub stars](https://img.shields.io/github/stars/OptimalDesignLab/PETSc2.jl)
![GitHub license](https://img.shields.io/github/license/OptimalDesignLab/PETSc2.jl)

This package provides thin wrappers for PETSc, as well as a few convenience functions that take advantage of multiple dispatch.

### [PositiveFactorizations.jl](https://github.com/timholy/PositiveFactorizations.jl)

![GitHub contributors](https://img.shields.io/github/contributors/timholy/PositiveFactorizations.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/timholy/PositiveFactorizations.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/timholy/PositiveFactorizations.jl)
![GitHub stars](https://img.shields.io/github/stars/timholy/PositiveFactorizations.jl)
![GitHub license](https://img.shields.io/github/license/timholy/PositiveFactorizations.jl)

PositiveFactorizations is a package for computing a positive definite matrix decomposition (factorization) from an arbitrary symmetric input. The motivating application is optimization.

### [FiniteDiff.jl](https://github.com/JuliaDiff/FiniteDiff.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaDiff/FiniteDiff.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaDiff/FiniteDiff.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaDiff/FiniteDiff.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaDiff/FiniteDiff.jl)
![GitHub license](https://img.shields.io/github/license/JuliaDiff/FiniteDiff.jl)

This package is for calculating derivatives, gradients, Jacobians, Hessians, etc. numerically.

### [SparseDiffTools.jl](https://github.com/JuliaDiff/SparseDiffTools.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaDiff/SparseDiffTools.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaDiff/SparseDiffTools.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaDiff/SparseDiffTools.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaDiff/SparseDiffTools.jl)
![GitHub license](https://img.shields.io/github/license/JuliaDiff/SparseDiffTools.jl)

This package is for exploiting sparsity in Jacobians and Hessians to accelerate computations. Matrix-free Jacobian-vector product and Hessian-vector product operators are provided that are compatible with AbstractMatrix-based libraries like IterativeSolvers.jl for easy and efficient Newton-Krylov implementation. It is possible to perform matrix coloring, and utilize coloring in Jacobian and Hessian construction.

### [HierarchicalMatrices.jl](https://github.com/JuliaMatrices/HierarchicalMatrices.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaMatrices/HierarchicalMatrices.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaMatrices/HierarchicalMatrices.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaMatrices/HierarchicalMatrices.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaMatrices/HierarchicalMatrices.jl)
![GitHub license](https://img.shields.io/github/license/JuliaMatrices/HierarchicalMatrices.jl)

This package provides a flexible framework for hierarchical matrices in Julia.

### [LowRankApprox.jl](https://github.com/JuliaMatrices/LowRankApprox.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaMatrices/LowRankApprox.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaMatrices/LowRankApprox.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaMatrices/LowRankApprox.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaMatrices/LowRankApprox.jl)
![GitHub license](https://img.shields.io/github/license/JuliaMatrices/LowRankApprox.jl)

This Julia package provides fast low-rank approximation algorithms for BLAS/LAPACK-compatible matrices based on some of the latest technology in adaptive randomized matrix sketching.

### [LazyBandedMatrices.jl](https://github.com/JuliaMatrices/LazyBandedMatrices.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaMatrices/LazyBandedMatrices.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaMatrices/LazyBandedMatrices.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaMatrices/LazyBandedMatrices.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaMatrices/LazyBandedMatrices.jl)
![GitHub license](https://img.shields.io/github/license/JuliaMatrices/LazyBandedMatrices.jl)

This package supports lazy banded and block-banded matrices.

### [BlockBandedMatrices.jl](https://github.com/JuliaMatrices/BlockBandedMatrices.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaMatrices/BlockBandedMatrices.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaMatrices/BlockBandedMatrices.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaMatrices/BlockBandedMatrices.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaMatrices/BlockBandedMatrices.jl)
![GitHub license](https://img.shields.io/github/license/JuliaMatrices/BlockBandedMatrices.jl)

A Julia package for representing block-block-banded matrices and banded-block-banded matrices.

### [kernelmatrices.jl](https://bitbucket.org/cgeoga/kernelmatrices.jl)

This software suite is a companion to the manuscript Scalable Gaussian Process Computations using Hierarchical Matrices.

### [ClusterTrees.jl](https://github.com/krcools/ClusterTrees.jl)

![GitHub contributors](https://img.shields.io/github/contributors/krcools/ClusterTrees.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/krcools/ClusterTrees.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/krcools/ClusterTrees.jl)
![GitHub stars](https://img.shields.io/github/stars/krcools/ClusterTrees.jl)
![GitHub license](https://img.shields.io/github/license/krcools/ClusterTrees.jl)

Tree data structures for fast multipole methods and H-matrices.

### [Arpack.jl](https://github.com/JuliaLinearAlgebra/Arpack.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaLinearAlgebra/Arpack.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaLinearAlgebra/Arpack.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaLinearAlgebra/Arpack.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaLinearAlgebra/Arpack.jl)
![GitHub license](https://img.shields.io/github/license/JuliaLinearAlgebra/Arpack.jl)

Julia wrapper for the arpack library designed to solve large scale eigenvalue problems. ARPACK is a collection of Fortran77 subroutines designed to solve large scale eigenvalue problems.

### [AlgebraicMultigrid.jl](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaLinearAlgebra/AlgebraicMultigrid.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaLinearAlgebra/AlgebraicMultigrid.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaLinearAlgebra/AlgebraicMultigrid.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaLinearAlgebra/AlgebraicMultigrid.jl)
![GitHub license](https://img.shields.io/github/license/JuliaLinearAlgebra/AlgebraicMultigrid.jl)

Solve sparse linear systems using Algebraic Multigrid (AMG). This works especially well for symmetric positive definite matrices.

### [NonlinearEigenproblems.jl](https://github.com/nep-pack/NonlinearEigenproblems.jl)

![GitHub contributors](https://img.shields.io/github/contributors/nep-pack/NonlinearEigenproblems.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/nep-pack/NonlinearEigenproblems.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/nep-pack/NonlinearEigenproblems.jl)
![GitHub stars](https://img.shields.io/github/stars/nep-pack/NonlinearEigenproblems.jl)
![GitHub license](https://img.shields.io/github/license/nep-pack/NonlinearEigenproblems.jl)

This package aims to provide state-of-the-art algorithms to solve the nonlinear eigenvalue problem. This currently includes (but is not restricted to) Newton-type methods, Subspace methods, Krylov methods, contour integral methods, block methods, companion matrix approaches. Problem transformation techniques such as scaling, shifting, deflating are also natively supported by the package.

### [BifurcationKit.jl](https://github.com/rveltz/BifurcationKit.jl)

![GitHub contributors](https://img.shields.io/github/contributors/rveltz/BifurcationKit.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/rveltz/BifurcationKit.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/rveltz/BifurcationKit.jl)
![GitHub stars](https://img.shields.io/github/stars/rveltz/BifurcationKit.jl)
![GitHub license](https://img.shields.io/github/license/rveltz/BifurcationKit.jl)

This organisation provides a set of tools to compute (automatically) bifurcation diagrams of ODE, DDE, PDE, nonlocal problems, etc. It can thus tackle large dimensional problems, by for example, using a GPU. The computation of periodic orbits and their bifurcations is also possible.

### [KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl)

![GitHub contributors](https://img.shields.io/github/contributors/Jutho/KrylovKit.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/Jutho/KrylovKit.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/Jutho/KrylovKit.jl)
![GitHub stars](https://img.shields.io/github/stars/Jutho/KrylovKit.jl)
![GitHub license](https://img.shields.io/github/license/Jutho/KrylovKit.jl)

A Julia package collecting a number of Krylov-based algorithms for linear problems, singular value and eigenvalue problems and the application of functions of linear maps or operators to vectors.

### [LinearOperators.jl](https://github.com/JuliaSmoothOptimizers/LinearOperators.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaSmoothOptimizers/LinearOperators.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaSmoothOptimizers/LinearOperators.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaSmoothOptimizers/LinearOperators.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaSmoothOptimizers/LinearOperators.jl)
![GitHub license](https://img.shields.io/github/license/JuliaSmoothOptimizers/LinearOperators.jl)

Operators behave like matrices, but are defined by their effect when applied to a vector. They can be transposed, conjugated, or combined with other operators cheaply. The costly operation is deferred until multiplied with a vector.

### [LinearMaps.jl](https://github.com/Jutho/LinearMaps.jl)

![GitHub contributors](https://img.shields.io/github/contributors/Jutho/LinearMaps.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/Jutho/LinearMaps.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/Jutho/LinearMaps.jl)
![GitHub stars](https://img.shields.io/github/stars/Jutho/LinearMaps.jl)
![GitHub license](https://img.shields.io/github/license/Jutho/LinearMaps.jl)

A Julia package for defining and working with linear maps, also known as linear transformations or linear operators acting on vectors. The only requirement for a LinearMap is that it can act on a vector (by multiplication) efficiently.

### [AbstractOperators.jl](https://github.com/kul-forbes/AbstractOperators.jl)

![GitHub contributors](https://img.shields.io/github/contributors/kul-forbes/AbstractOperators.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/kul-forbes/AbstractOperators.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/kul-forbes/AbstractOperators.jl)
![GitHub stars](https://img.shields.io/github/stars/kul-forbes/AbstractOperators.jl)
![GitHub license](https://img.shields.io/github/license/kul-forbes/AbstractOperators.jl)

Abstract operators extend the syntax typically used for matrices to linear mappings of arbitrary dimensions and nonlinear functions. Unlike matrices however, abstract operators apply the mappings with specific efficient algorithms that minimize memory requirements. This is particularly useful in iterative algorithms and in first order large-scale optimization algorithms.

### [MUMPS.jl](https://github.com/JuliaSmoothOptimizers/MUMPS.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaSmoothOptimizers/MUMPS.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaSmoothOptimizers/MUMPS.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaSmoothOptimizers/MUMPS.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaSmoothOptimizers/MUMPS.jl)
![GitHub license](https://img.shields.io/github/license/JuliaSmoothOptimizers/MUMPS.jl)

MUMPS is a library for the solution of large linear systems using a factorization. Structure can be exploited, such as symmetry, or symmetry and definiteness. The factorization and solve phases can be performed in parallel.

### [Pardiso.jl](https://github.com/JuliaSparse/Pardiso.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaSparse/Pardiso.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaSparse/Pardiso.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaSparse/Pardiso.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaSparse/Pardiso.jl)
![GitHub license](https://img.shields.io/github/license/JuliaSparse/Pardiso.jl)

The Pardiso.jl package provides an interface for using PARDISO 6.0 and Intel MKL PARDISO from the Julia language.

### [Metis.jl](https://github.com/JuliaSparse/Metis.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaSparse/Metis.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaSparse/Metis.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaSparse/Metis.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaSparse/Metis.jl)
![GitHub license](https://img.shields.io/github/license/JuliaSparse/Metis.jl)

Metis.jl is a Julia wrapper to the Metis library which is a library for partitioning unstructured graphs, partitioning meshes, and computing fill-reducing orderings of sparse matrices.

### [IterativeSolvers.jl](https://github.com/JuliaMath/IterativeSolvers.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaMath/IterativeSolvers.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaMath/IterativeSolvers.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaMath/IterativeSolvers.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaMath/IterativeSolvers.jl)
![GitHub license](https://img.shields.io/github/license/JuliaMath/IterativeSolvers.jl)

IterativeSolvers is a Julia package that provides iterative algorithms for solving linear systems, eigensystems, and singular value problems.

### [Preconditioners.jl](https://github.com/mohamed82008/Preconditioners.jl)

![GitHub contributors](https://img.shields.io/github/contributors/mohamed82008/Preconditioners.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/mohamed82008/Preconditioners.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/mohamed82008/Preconditioners.jl)
![GitHub stars](https://img.shields.io/github/stars/mohamed82008/Preconditioners.jl)
![GitHub license](https://img.shields.io/github/license/mohamed82008/Preconditioners.jl)

Selected pre-conditioners for iterative  solvers of coupled linear algebraic equations. Including AMG.

## <a name="geo">Geometry and topology

### [Manifolds.jl](https://github.com/JuliaNLSolvers/Manifolds.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaNLSolvers/Manifolds.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaNLSolvers/Manifolds.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaNLSolvers/Manifolds.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaNLSolvers/Manifolds.jl)
![GitHub license](https://img.shields.io/github/license/JuliaNLSolvers/Manifolds.jl)

Manifolds.jl aims to provide both a unified interface to define and use manifolds as well as a library of manifolds to use for your projects.

### [Grassmann.jl](https://github.com/chakravala/Grassmann.jl)

![GitHub contributors](https://img.shields.io/github/contributors/chakravala/Grassmann.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/chakravala/Grassmann.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/chakravala/Grassmann.jl)
![GitHub stars](https://img.shields.io/github/stars/chakravala/Grassmann.jl)
![GitHub license](https://img.shields.io/github/license/chakravala/Grassmann.jl)

The Grassmann.jl package provides tools for doing computations based on multi-linear algebra, differential geometry, and spin groups using the extended tensor algebra known as Leibniz-Grassmann-Clifford-Hestenes geometric algebra.

### [DiscreteDifferentialGeometry.jl](https://github.com/mewertd2/DiscreteDifferentialGeometry.jl)

![GitHub contributors](https://img.shields.io/github/contributors/mewertd2/DiscreteDifferentialGeometry.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/mewertd2/DiscreteDifferentialGeometry.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/mewertd2/DiscreteDifferentialGeometry.jl)
![GitHub stars](https://img.shields.io/github/stars/mewertd2/DiscreteDifferentialGeometry.jl)
![GitHub license](https://img.shields.io/github/license/mewertd2/DiscreteDifferentialGeometry.jl)

The DiscreteDifferentialGeometry Julia package defines types and methods to implement Discrete Differential Geometry.

## <a name="grids">Mesh and Grid Generation

### [Gmsh.jl](https://github.com/JuliaFEM/Gmsh.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaFEM/Gmsh.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaFEM/Gmsh.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaFEM/Gmsh.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaFEM/Gmsh.jl)
![GitHub license](https://img.shields.io/github/license/JuliaFEM/Gmsh.jl)

Gmsh.jl contains API for Gmsh: a three-dimensional finite element mesh generator. With the help of Gmsh.jl, it is possible add parametric model construction and/or automatic mesh generation to a FEM/FVM pipeline.

### [Triangulate.jl](https://github.com/JuliaGeometry/Triangulate.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaGeometry/Triangulate.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaGeometry/Triangulate.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaGeometry/Triangulate.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaGeometry/Triangulate.jl)
![GitHub license](https://img.shields.io/github/license/JuliaGeometry/Triangulate.jl)

Julia wrapper for Johnathan Richard Shewchuk's Triangle mesh generator.

### [TetGen.jl](https://github.com/JuliaGeometry/TetGen.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaGeometry/TetGen.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaGeometry/TetGen.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaGeometry/TetGen.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaGeometry/TetGen.jl)
![GitHub license](https://img.shields.io/github/license/JuliaGeometry/TetGen.jl)

The TetGen.jl package is a Julia wrapper for the C++ project TetGen. This wrapper enables TetGen based tetrahedral meshing, and (constrained) 3D Delaunay and Voronoi tesselation.

### [CompScienceMeshes.jl](https://github.com/krcools/CompScienceMeshes.jl)

![GitHub contributors](https://img.shields.io/github/contributors/krcools/CompScienceMeshes.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/krcools/CompScienceMeshes.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/krcools/CompScienceMeshes.jl)
![GitHub stars](https://img.shields.io/github/stars/krcools/CompScienceMeshes.jl)
![GitHub license](https://img.shields.io/github/license/krcools/CompScienceMeshes.jl)

Geometry types and algorithms for computational science. Meshes, charts, and neighborhoods.

### [VoronoiDelaunay.jl](https://github.com/JuliaGeometry/VoronoiDelaunay.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaGeometry/VoronoiDelaunay.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaGeometry/VoronoiDelaunay.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaGeometry/VoronoiDelaunay.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaGeometry/VoronoiDelaunay.jl)
![GitHub license](https://img.shields.io/github/license/JuliaGeometry/VoronoiDelaunay.jl)

Fast, robust construction of 2D Delaunay and Voronoi tessellations on generic point types.

### [MiniQhull.jl](https://github.com/gridap/MiniQhull.jl)

![GitHub contributors](https://img.shields.io/github/contributors/gridap/MiniQhull.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/gridap/MiniQhull.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/gridap/MiniQhull.jl)
![GitHub stars](https://img.shields.io/github/stars/gridap/MiniQhull.jl)
![GitHub license](https://img.shields.io/github/license/gridap/MiniQhull.jl)

A simple Julia wrapper around Qhull to compute Delaunay triangulations in arbitrary dimensions. This is a direct wrapper using `ccall` (in constrast to `QHull.jl` that is a wrapper around a python wrapper of Qhull.)

### [FMMTrees.jl](https://github.com/krcools/FMMTrees.jl)

![GitHub contributors](https://img.shields.io/github/contributors/krcools/FMMTrees.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/krcools/FMMTrees.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/krcools/FMMTrees.jl)
![GitHub stars](https://img.shields.io/github/stars/krcools/FMMTrees.jl)
![GitHub license](https://img.shields.io/github/license/krcools/FMMTrees.jl)

Tree data structures for H(2), hierarchical matrices, and FMM-like algorithms (fast multiple methods).

### [DistMesh.jl](https://github.com/JuliaGeometry/DistMesh.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaGeometry/DistMesh.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaGeometry/DistMesh.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaGeometry/DistMesh.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaGeometry/DistMesh.jl)
![GitHub license](https://img.shields.io/github/license/JuliaGeometry/DistMesh.jl)

Tetrahedral mesh refinement on signed distance/implicit functions or level sets using TetGen.

### [HalfEdges.jl](https://github.com/digitaldomain/HalfEdges.jl)

![GitHub contributors](https://img.shields.io/github/contributors/digitaldomain/HalfEdges.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/digitaldomain/HalfEdges.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/digitaldomain/HalfEdges.jl)
![GitHub stars](https://img.shields.io/github/stars/digitaldomain/HalfEdges.jl)
![GitHub license](https://img.shields.io/github/license/digitaldomain/HalfEdges.jl)

Halfedge data structure for navigating and querying polygonal meshes.

### [MeshCore.jl](https://github.com/PetrKryslUCSD/MeshCore.jl)

![GitHub contributors](https://img.shields.io/github/contributors/PetrKryslUCSD/MeshCore.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/PetrKryslUCSD/MeshCore.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/PetrKryslUCSD/MeshCore.jl)
![GitHub stars](https://img.shields.io/github/stars/PetrKryslUCSD/MeshCore.jl)
![GitHub license](https://img.shields.io/github/license/PetrKryslUCSD/MeshCore.jl)

Small package for the topology of meshes for the Finite Element Methods (FEM). All essential topological incidence relations are provided: see the guide. The library provides efficient storage of coordinates and connectivities in static arrays for speed of access.

### [MeshSteward.jl](https://github.com/PetrKryslUCSD/MeshSteward.jl)

![GitHub contributors](https://img.shields.io/github/contributors/PetrKryslUCSD/MeshSteward.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/PetrKryslUCSD/MeshSteward.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/PetrKryslUCSD/MeshSteward.jl)
![GitHub stars](https://img.shields.io/github/stars/PetrKryslUCSD/MeshSteward.jl)
![GitHub license](https://img.shields.io/github/license/PetrKryslUCSD/MeshSteward.jl)

Small package for the management of meshes for the Finite Element Methods (FEM). A convenience layer on top of [MeshCore.jl](https://github.com/PetrKryslUCSD/MeshCore.jl).
 
 ### [DelaunayTriangulation.jl](https://github.com/DanielVandH/DelaunayTriangulation.jl)

![GitHub contributors](https://img.shields.io/github/contributors/DanielVandH/DelaunayTriangulation.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/DanielVandH/DelaunayTriangulation.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/DanielVandH/DelaunayTriangulation.jl)
![GitHub stars](https://img.shields.io/github/stars/DanielVandH/DelaunayTriangulation.jl)
![GitHub license](https://img.shields.io/github/license/DanielVandH/DelaunayTriangulation.jl)
 
Package for constructing unconstrained or constrained Delaunay triangulations of planar point sets, with support for mesh refinement with constant or custom area and angle constraints, with a generic interface for geometric primitives to allow for easy customisation. Exact predicates are used to make the procedures robust.

## <a name="post"></a>Postprocessing, visualization

### [WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl).

![GitHub contributors](https://img.shields.io/github/contributors/jipolanco/WriteVTK.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/jipolanco/WriteVTK.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/jipolanco/WriteVTK.jl)
![GitHub stars](https://img.shields.io/github/stars/jipolanco/WriteVTK.jl)
![GitHub license](https://img.shields.io/github/license/jipolanco/WriteVTK.jl)

This module allows to write VTK XML files, that can be visualised for example with ParaView. Seems pretty complete, writes compressed files.

### [VTKDataTypes.jl](https://github.com/mohamed82008/VTKDataTypes.jl)

![GitHub contributors](https://img.shields.io/github/contributors/mohamed82008/VTKDataTypes.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/mohamed82008/VTKDataTypes.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/mohamed82008/VTKDataTypes.jl)
![GitHub stars](https://img.shields.io/github/stars/mohamed82008/VTKDataTypes.jl)
![GitHub license](https://img.shields.io/github/license/mohamed82008/VTKDataTypes.jl)

VTKDataTypes.jl presents a Julia type system for representing and manipulating VTK data natively in Julia.

## <a name="hpc"></a>HPC, Parallel processing

### [MPI.jl](https://github.com/JuliaParallel/MPI.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaParallel/MPI.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaParallel/MPI.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaParallel/MPI.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaParallel/MPI.jl)
![GitHub license](https://img.shields.io/github/license/JuliaParallel/MPI.jl)

This provides Julia interface to the Message Passing Interface (MPI), roughly inspired by mpi4py.

### [PumiInterface.jl](https://github.com/OptimalDesignLab/PumiInterface.jl)

![GitHub contributors](https://img.shields.io/github/contributors/OptimalDesignLab/PumiInterface.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/OptimalDesignLab/PumiInterface.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/OptimalDesignLab/PumiInterface.jl)
![GitHub stars](https://img.shields.io/github/stars/OptimalDesignLab/PumiInterface.jl)
![GitHub license](https://img.shields.io/github/license/OptimalDesignLab/PumiInterface.jl)

This code provides a way to use PUMI from Julia by wrapping functions in PUMIs APF API.


## <a name="opt"></a>Optimization

### [TopOpt.jl](https://github.com/mohamed82008/TopOpt.jl)

![GitHub contributors](https://img.shields.io/github/contributors/mohamed82008/TopOpt.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/mohamed82008/TopOpt.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/mohamed82008/TopOpt.jl)
![GitHub stars](https://img.shields.io/github/stars/mohamed82008/TopOpt.jl)
![GitHub license](https://img.shields.io/github/license/mohamed82008/TopOpt.jl)

A WIP topology optimization package in pure Julia.

## <a name="misc"></a>Miscellanea

### [AD4SM.jl](https://github.com/avigliotti/AD4SM.jl)

![GitHub contributors](https://img.shields.io/github/contributors/avigliotti/AD4SM.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/avigliotti/AD4SM.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/avigliotti/AD4SM.jl)
![GitHub stars](https://img.shields.io/github/stars/avigliotti/AD4SM.jl)
![GitHub license](https://img.shields.io/github/license/avigliotti/AD4SM.jl)

Automatic Differentiation for Solid Mechanics in Julia.

### [RHEOS.jl](https://github.com/JuliaRheology/RHEOS.jl)

![GitHub contributors](https://img.shields.io/github/contributors/JuliaRheology/RHEOS.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JuliaRheology/RHEOS.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/JuliaRheology/RHEOS.jl)
![GitHub stars](https://img.shields.io/github/stars/JuliaRheology/RHEOS.jl)
![GitHub license](https://img.shields.io/github/license/JuliaRheology/RHEOS.jl)

RHEOS, an abbreviation of Rheology Open Source, is a software package written in the Julia programming language that provides tools for analyzing rheological data.

### [Tensors.jl](https://github.com/KristofferC/Tensors.jl)

![GitHub contributors](https://img.shields.io/github/contributors/KristofferC/Tensors.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/KristofferC/Tensors.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/KristofferC/Tensors.jl)
![GitHub stars](https://img.shields.io/github/stars/KristofferC/Tensors.jl)
![GitHub license](https://img.shields.io/github/license/KristofferC/Tensors.jl)

Efficient computations with symmetric and non-symmetric tensors with support for automatic differentiation.

### [TensorOperations.jl](https://github.com/Jutho/TensorOperations.jl)

![GitHub contributors](https://img.shields.io/github/contributors/Jutho/TensorOperations.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/Jutho/TensorOperations.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/Jutho/TensorOperations.jl)
![GitHub stars](https://img.shields.io/github/stars/Jutho/TensorOperations.jl)
![GitHub license](https://img.shields.io/github/license/Jutho/TensorOperations.jl)

Fast tensor operations using a convenient Einstein index notation.

### [TensorToolbox.jl](https://github.com/lanaperisa/TensorToolbox.jl)

![GitHub contributors](https://img.shields.io/github/contributors/lanaperisa/TensorToolbox.jl)
![GitHub closed issues](https://img.shields.io/github/issues-closed/lanaperisa/TensorToolbox.jl)
![GitHub last commit](https://img.shields.io/github/last-commit/lanaperisa/TensorToolbox.jl)
![GitHub stars](https://img.shields.io/github/stars/lanaperisa/TensorToolbox.jl)
![GitHub license](https://img.shields.io/github/license/lanaperisa/TensorToolbox.jl)

Julia package for tensors (dense tensors,
Tucker format,
Kruskal (CP) format,
Hierarchical Tucker format,
Tensor Train format).

### [NonlinearPDE](https://github.com/JerryLingjieMei/NonlinearPDE)

![GitHub contributors](https://img.shields.io/github/contributors/JerryLingjieMei/NonlinearPDE)
![GitHub closed issues](https://img.shields.io/github/issues-closed/JerryLingjieMei/NonlinearPDE)
![GitHub last commit](https://img.shields.io/github/last-commit/JerryLingjieMei/NonlinearPDE)
![GitHub stars](https://img.shields.io/github/stars/JerryLingjieMei/NonlinearPDE)
![GitHub license](https://img.shields.io/github/license/JerryLingjieMei/NonlinearPDE)

Solution of a nonlinear elliptic PDE with multi-grid.

## Waiting for classification

