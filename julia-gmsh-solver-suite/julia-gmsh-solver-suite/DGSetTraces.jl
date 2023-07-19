#=-------------------------------------------------------------------
 Dependencies
-------------------------------------------------------------------=#


using LinearAlgebra
using StaticArrays
using Plots

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

using GLMakie

if gmsh.isInitialized() == 1
    gmsh.finalize()
end
gmsh.initialize()


# problem geometry
mesh_filename = "./CircleMeshInitialDG.msh"
# mesh_filename = "./CircleMeshVeryFine.msh"

setMeshOrder = 2
@time mesh = setTracesForEdges(mesh_filename,setMeshOrder)
# print(mesh)
# mesh_connectivity = setupMeshConnectivity(mesh)
