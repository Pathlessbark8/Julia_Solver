#=-------------------------------------------------------------------
 Dependencies
-------------------------------------------------------------------=#

using Revise

includet("Mesh/Mesh.jl")
using .Meshes
Revise.track(Meshes, "Mesh/AbstractMesh.jl")
Revise.track(Meshes, "Mesh/GmshMeshUtils.jl")
Revise.track(Meshes, "Mesh/SimpleMesh.jl")

includet("FunctionSpaces/FunctionSpaces.jl")
using .FunctionSpaces
Revise.track(FunctionSpaces, "FunctionSpaces/HOneElement.jl")
Revise.track(FunctionSpaces, "FunctionSpaces/HOneSpace.jl")
Revise.track(FunctionSpaces, "FunctionSpaces/HCurlElement.jl")
Revise.track(FunctionSpaces, "FunctionSpaces/HCurlSpace.jl")

includet("Physics/Physics.jl")
using .Physics
Revise.track(Physics, "Physics/AbstractPhysics.jl")
Revise.track(Physics, "Physics/EMPhysics.jl")


includet("Sources/Sources.jl")
using .Sources
Revise.track(Sources, "Sources/AbstractSources.jl")
Revise.track(Sources, "Sources/EMSources.jl")

includet("Solvers/Solvers.jl")
using .Solvers
Revise.track(Solvers, "Solvers/AbstractSolvers.jl")
Revise.track(Solvers, "Solvers/EMSolverUtilities.jl")
Revise.track(Solvers, "Solvers/EMTMFEMSolver.jl")

includet("Extras/Visualization.jl")
using .Visualization
Revise.track(Visualization, "Extras/MakieViz.jl")

using Gmsh: gmsh
using SparseArrays
using StaticArrays
using LinearAlgebra

include("Extras/Analytic.jl")

using CairoMakie

include("src/EIL.jl")

#=-------------------------------------------------------------------
 Main 
-------------------------------------------------------------------=#

#=------------------------------------
 Initialize Gmsh
-------------------------------------=#

if gmsh.isInitialized() == 1
    gmsh.finalize()
end
gmsh.initialize()

#=------------------------------------
 Problem definition (user input proxy)
-------------------------------------=#

# basic formulation
problem_type = EMTMProblemType()
solution_order = 1

# problem geometry
mesh_filename = "./CircleMeshFine.msh"

#Frequency Input
frequency = 3e8

#Constitutives Input
constitutives_table = Dict{Int64,Any}()
constitutives_table[3000] = SimpleComplexEMConstitutives("Freespace", 3000, 1.0 + 0.0im, 1.0 + 0.0im)
constitutives_table[3001] = SimpleComplexEMConstitutives("Target", 3001, 1.0 + 0.0im, 1.0 + 0.0im)
constitutives_table[3002] = SimpleComplexEMConstitutives("Other", 3002, 1.0 + 0.0im, 1.0 + 0.0im)

boundary_condition_table = Dict{Int64, Physics.AbstractEMBoundaryCondition}()
boundary_condition_table[2000] = ABCEMBoundaryCondition()

sources = Vector{Sources.AbstractEMSource}(undef,0)
push!(sources, Sources.EMTMGaussianCurrentDensity(0.5, 0.5, 0.001, 0.001, 1.0))

#=------------------------------------
# Mesh 
------------------------------------=#

mesh = loadMesh(mesh_filename)
mesh_connectivity = setupMeshConnectivity(mesh)

#=------------------------------------
 Constitutives
------------------------------------=#
element_constitutives = setupElementConstitutives(problem_type, mesh, constitutives_table)

#=------------------------------------
 Boundary Conditions
------------------------------------=#

element_boundary_conditions = setupBoundaryConditions(problem_type, mesh, mesh_connectivity, boundary_condition_table)
 
#=------------------------------------
 Basis / Finite Element Space
------------------------------------=#

function_space = HOneSpace(mesh, HOneElement(mesh.properties.dim, solution_order)) #change or add constructor to auto generate element based on order
element_matrices = elementVolumeMatrices(problem_type, function_space)

gmsh.finalize()

