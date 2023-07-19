module FunctionSpaces

using Gmsh: gmsh

# import ..Meshes
using ..Meshes: AbstractMesh, SimpleMesh, numNodes, numElements, createFaces, createEdges, Line, Triangle, Tetrahedron, AbstractElement

#using ..Physics: DirichletBoundaryCondition, NeumannBoundaryCondition, RobinBoundaryCondition


abstract type AbstractFiniteElement{T<:AbstractElement} end
abstract type AbstractFunctionSpace{T<:AbstractFiniteElement} end

# commented out because we don't need to necessarily plan for future support of non HCurl or H1 elements
# for vector functions with scalar coefficients
# abstract type VectorFiniteElement{T} <: FiniteElement{T} end

# for scalar function with scalar coefficients
# abstract type ScalarFiniteElement{T} <: FiniteElement{T} end


include("HOneElement.jl")
export HOneElement, HOneLIne, HOneTriangle, HOneTetrahedron

include("HCurlElement.jl")
export HCurlElement, HCurlLine, HCurlTriangle, HCurlTetrahedron, numEdgeFunctions, numFaceFunctions, numBubbleFunctions, numEdges, numFaces, numVolumes

struct SimpleFiniteElementSpace{T} <: AbstractFunctionSpace{T}
    mesh::SimpleMesh
    local_space::T
    dof_map::Array{Int64,2} # (nFunctions, nElements)
    physical_group_mapping::Dict{Tuple{Int32,Int32},Array{Int32,2}}
end
numDOF(space::SimpleFiniteElementSpace) = maximum(space.dof_map)
export SimpleFiniteElementSpace
export numDOF

function getPhysicsDOF(space::SimpleFiniteElementSpace, dim::Integer, physics_tag::Integer)
    key = (dim, physics_tag)
    return vec(space.physical_group_mapping[key])
end
export getPhysicsDOF


include("HCurlSpace.jl")
export HCurlSpace, numDOF

export massMatrix, curlMatrix, integrateRobinBoundary, integratePointSource

include("HOneSpace.jl")
export HOneSpace, integrate

end