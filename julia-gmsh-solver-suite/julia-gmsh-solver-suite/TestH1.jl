using Revise
includet("Mesh/Mesh.jl")
using .Meshes
Revise.track(Meshes, "Mesh/GmshMeshUtils.jl")
Revise.track(Meshes, "Mesh/SimpleMesh.jl")

includet("Physics/Physics.jl")
using .Physics

includet("FunctionSpaces/FunctionSpaces.jl")
using .FunctionSpaces
Revise.track(FunctionSpaces, "FunctionSpaces/H1Element.jl")

includet("Extras/Visualization.jl")
using .Visualization
Revise.track(Visualization, "Extras/MakieViz.jl")

using Gmsh: gmsh
using SparseArrays
using StaticArrays
using LinearAlgebra

using CairoMakie

# reset gmsh
if gmsh.isInitialized() == 1
    gmsh.finalize()
end

gmsh.initialize()

basis_order = 3

local_element = H1Element(2, basis_order)

points = local_element.quad_points
funcs = local_element.functions

x = points[1:3:end]
y = points[2:3:end]

f = Figure()

for i =1:3
    z = funcs[i, :, 1]
    ax = Axis(f[1, i], aspect=1)
    tr = tricontourf!(ax, x, y, z)
    # scatter!(x, y, color=z, strokewidth=1, strokecolor = :black)
end

nEdge = 3*basis_order-3
for i = 1:nEdge
    z = funcs[3+i, :, 1]
    ax = Axis(f[2, i], aspect=1)
    tr = tricontourf!(ax, x, y, z)
    # scatter!(x, y, color=z, strokewidth=1, strokecolor = :black)
end

nFace = convert(Int64, (basis_order - 1) * (basis_order - 2) / 2)
for i = 1:nFace
    z = funcs[3+nEdge+i,:,1]
    ax = Axis(f[3,i], aspect=1)
    tr = tricontourf!(ax, x, y, z)
end
display(f)

gmsh.finalize()
