module Solvers

abstract type AbstractElementVolumeMatrices end
abstract type AbstractBoundaryFluxMatrices end

export AbstractElementVolumeMatrices
export AbstractBoundaryFluxMatrices

using Gmsh: gmsh

using ..Meshes: AbstractMesh, SimpleMesh, SimpleMeshConnectivity, numNodes, numElements, createFaces, createEdges, Line, Triangle, Tetrahedron, AbstractElement

using ..FunctionSpaces: AbstractFiniteElement, SimpleFiniteElementSpace, HOneElement, HOneSpace, HOneTriangle

using ..Physics: AbstractProblemType, EMTMProblemType, EMTEProblemType, EM3DProblemType, 
                    DirichletBoundaryCondition, NeumannBoundaryCondition, RobinBoundaryCondition, 
                    ComplexEMElementConstitutives, 
                    AbstractEMBoundaryCondition, PECEMBoundaryCondition, PMCEMBoundaryCondition, ABCEMBoundaryCondition, EMElementBoundaryCondition

include("AbstractSolvers.jl")
include("EMSolverUtilities.jl")

export setupElementConstitutives, setupBoundaryConditions

include("EMTMFEMSolver.jl")

export EMTMElementVolumeMatrices, EMTMBoundaryFluxMatrices, elementVolumeMatrices


include("EMDGMSolver.jl")

end