import Gmsh:gmsh #Julia package available for API that was installed
import LinearAlgebra

#using PyPlot
using Plots
using Base.Iterators
#using GLMakie

#if the code crashes and gmsh doesn't get finalized, problems ensue on next run
try
    gmsh.finalize()
catch

end

#start up the gmsh api and set a few of my own GUI display preferences
gmsh.initialize()
gmsh.option.setNumber("Geometry.Volumes", 1)
gmsh.option.setNumber("Geometry.Surfaces", 1)
gmsh.option.setNumber("Geometry.Normals", 30)

gmsh.open("./PECCircle.msh")

problem_dimension = 2
geometric_order = 1;

tetrahedralType = gmsh.model.mesh.getElementType("tetrahedron", geometric_order)
triangleType = gmsh.model.mesh.getElementType("triangle",geometric_order)
elementType = -1




if problem_dimension == 3
    elementType = tetrahedralType
    nodesPerElement = 4
    facesPerElement = 4
    edgesPerElement = 6
    nodesPerFace = 3
elseif problem_dimension == 2
    elementType = triangleType
    nodesPerElement = 3
    facesPerElement = 3
    edgesPerElement = 3
    nodesPerFace = 2 
    nodesPerEdge = 1
else
    error("Unknown Problem Dimension")
end



nodeTags, nodeCoords = gmsh.model.mesh.getNodes()
nodeCoords = reshape(nodeCoords, (3, length(nodeTags))) #column major so each point is in a column.
nodeTags2Index = Dict{Int64,Int64}() #should be unsigned int but hard to read
for i in eachindex(nodeTags)
    nodeTags2Index[nodeTags[i]] = i
end

elements = gmsh.model.mesh.getElementsByType(elementType)
elementTags = elements[1]
elementNodes = elements[2] 
numElements = length(elementTags)
elementNodes = reshape(elementNodes, (nodesPerElement, numElements)) #questions about this operation in relation to memory

elementTags2Index = Dict{Int64,Int64}() #should be unsigned int but hard to read
for i in eachindex(elementTags)
    elementTags2Index[elementTags[i]] = i
end


#testing vector basis functions on triangles

element_index = 1
element_tag = elementTags[element_index]
local_points_u = range(0.0,1.0,8)
local_points_v = range(0.0,1.0,8)
local_points = Array{Float64}(undef, 0)
dummy_data = Array{Float64}(undef, 0)

for u in local_points_u
    for v in local_points_v
        if (u + v <= 1.0 && u + v >= 0.0) 
            #println("Adding point (", u, ", ", v, ")")
            append!(local_points, [u, v, 0])
            append!(dummy_data, u + v)
        end
    end
end
n_local_points = length(local_points) ÷ 3

local_points_u = local_points[1:3:end]
local_points_v = local_points[2:3:end]

basis_order = 1
basis_index = 1
orientation_indexes = []

n_L_orientations = gmsh.model.mesh.getNumberOfOrientations(elementType,"Lagrange")
n_GradL_orientations = gmsh.model.mesh.getNumberOfOrientations(elementType,"GradLagrange")
n_HcurlL_orientations = gmsh.model.mesh.getNumberOfOrientations(elementType,"HcurlLegendre")

_ , L, _ = gmsh.model.mesh.getBasisFunctions(elementType, local_points, "Lagrange" * string(basis_order), orientation_indexes)
_ , GradL, _ = gmsh.model.mesh.getBasisFunctions(elementType, local_points, "GradLagrange" * string(basis_order), orientation_indexes)
_ , HcurlL, _ = gmsh.model.mesh.getBasisFunctions(elementType, local_points, "HcurlLegendre" * string(basis_order), orientation_indexes)

n_basis_L = length(L) ÷ n_local_points ÷ n_GradL_orientations
n_basis_GradL = length(GradL) ÷ n_local_points ÷ n_GradL_orientations ÷ 3
n_basis_HcurlL = length(HcurlL) ÷ n_local_points ÷ n_HcurlL_orientations ÷ 3





t = gmsh.view.add("Test")
#gmsh.view.addHomogeneousModelData(t, 0, gmsh.model.getCurrent(), "ElementData", elementTags, elementTags)

VectorData = [[1.0 1.0 1.0] for i in 1:length(elementTags)]
VectorData = collect(Iterators.Flatten(VectorData))
gmsh.view.addHomogeneousModelData(t, 0, gmsh.model.getCurrent(), "ElementData", elementTags, VectorData)

gmsh.view.write(t, "./Q.msh")

gmsh.open("./Q.msh")
gmsh.fltk.run()


gmsh.finalize()
