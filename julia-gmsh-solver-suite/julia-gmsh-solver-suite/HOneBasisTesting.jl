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

using GLMakie

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


# problem geometry
mesh_filename = "./CircleMesh.msh"


mesh = loadMesh(mesh_filename)
mesh_connectivity = setupMeshConnectivity(mesh)

#=------------------------------------
 Basis / Finite Element Space
------------------------------------=#

solution_order = 3
element_type = mesh.properties.type
println(element_type)

local_points_u = range(0.0,stop=1.0,length=100)
local_points_v = range(0.0,stop=1.0,length=100)
local_points = Array{Float64}(undef, 0)

for u in local_points_u
    for v in local_points_v
        if (u + v <= 1) 
            append!(local_points, [u, v, 0])
        end
    end
end

#local_points = [0.2 0.8 0 0.4 0.6 0]

num_local_points = length(local_points) ÷ 3

local_points

# generate global points in each element so we can determine unique dofs
num_comp_legendre, basis_functions_legendre, num_orient_legendre = gmsh.model.mesh.getBasisFunctions(element_type, local_points[:], "H1Legendre$solution_order")
num_comp_lagrange, basis_functions_lagrange, num_orient_lagrange = gmsh.model.mesh.getBasisFunctions(element_type, local_points[:], "Lagrange$solution_order")

num_basis_legendre = (((length(basis_functions_legendre) ÷ num_comp_legendre) ÷ num_local_points) ÷ num_orient_legendre)
num_basis_lagrange = (((length(basis_functions_lagrange) ÷ num_comp_lagrange) ÷ num_local_points) ÷ num_orient_lagrange)
basis_functions_legendre = reshape(basis_functions_legendre, num_comp_legendre, num_basis_legendre, num_local_points, num_orient_legendre)
basis_functions_lagrange = reshape(basis_functions_lagrange, num_comp_lagrange, num_basis_lagrange, num_local_points, num_orient_lagrange)

println(size(basis_functions_lagrange))
# layout = Layout(
#     title="Mt Bruno Elevation", autosize=false,
#     width=500, height=500,
#     margin=attr(l=65, r=50, b=65, t=90)
# )
component_index = 1
orientation_index = 1
println(length(local_points)÷ 3)
println(length(basis_functions_lagrange))
println(string(num_comp_legendre) * " " * string(num_basis_legendre) * " " * string(num_local_points) * " " * string(num_orient_legendre))

fig1 = GLMakie.Figure();
ax = [Axis(fig1[1,i]) for i in 1:num_basis_lagrange]
s=[];
s = [GLMakie.heatmap!(ax[basis_index],local_points[1:3:end], local_points[2:3:end], basis_functions_lagrange[component_index, basis_index, :, orientation_index]) for basis_index=1:num_basis_lagrange]
display(fig1)

typeof(basis_functions_lagrange)

sc=[]
fig3 = GLMakie.Figure();

for basis_index in 1:num_basis_lagrange
    axs = Axis(fig3[1, 1]); 
    aspect=(1, 1, 1)
    display(plot(local_points[1:3:end],local_points[2:3:end],basis_functions_lagrange[component_index, basis_index, :, orientation_index]))
end


# fig2 = Figure()
# s=[]
# s = [GLMakie.heatmap(local_points[1:3:end], local_points[2:3:end], basis_functions_legendre[component_index, basis_index, :, orientation_index]) for basis_index=1:num_basis_legendre]
# display(fig2)

# for basis_index in 1:num_basis_legendre
#     fig3=GLMakie.Figure()
#     display(GLMakie.plot(local_points[1:3:end],local_points[2:3:end],basis_functions_legendre[component_index, basis_index, :, orientation_index],st=:surface,c=my_cg))# sc = display(fig);
# end

#display(plot(s[1], s[2]))