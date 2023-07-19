#=============================================================
A first attempt at Julia-based GMsh for FEM and DGM solvers
Ian Jeffrey
December 21, 2023

Much of the work here follows the x7.py example here:
https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_11_1/tutorials/python/x7.py#L44
==============================================================#

import Gmsh:gmsh #Julia package available for API that was installed
import LinearAlgebra

#using PyPlot
using Plots
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

make_3d_geo = false
make_2d_geo = false
read_3d_geo = false
read_2d_geo = true
if make_3d_geo
    # create a box with PEC boundaries using the OpenCascade (OCC) interface
    box_id = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1);
    gmsh.model.occ.synchronize();
    gmsh.model.addPhysicalGroup(3, [box_id], 3000, "freespace")
    gmsh.model.geo.synchronize();
    box_boundary = gmsh.model.getBoundary((3,box_id))
    box_boundary = abs.(getfield.(box_boundary, 2))
    gmsh.model.addPhysicalGroup(2,box_boundary, 2000, "pec")
    gmsh.model.mesh.generate(3)
    gmsh.write("PECCube.msh")
    gmsh.fltk.run()
elseif make_2d_geo
    circle_id = gmsh.model.occ.addCircle(0,0,0, 1, -1, 0, 2*pi);
    gmsh.model.occ.synchronize()
    curve_loop_id = gmsh.model.geo.addCurveLoop([circle_id])
    surface_id = gmsh.model.geo.addPlaneSurface([curve_loop_id]);
    gmsh.model.geo.synchronize();
    gmsh.model.addPhysicalGroup(2, [surface_id], 3000, "freespace")
    circle_boundary = gmsh.model.getBoundary((2,circle_id))
    circle_boundary = abs.(getfield.(circle_boundary, 2))
    gmsh.model.addPhysicalGroup(1, circle_boundary, 2000, "pec")
    gmsh.model.mesh.generate(2)
    gmsh.write("PECCircle.msh")
    gmsh.fltk.run()
elseif read_3d_geo
    gmsh.open("./PECCube.msh")
elseif read_2d_geo
    gmsh.open("./PECCircle.msh")
end

#==================================
 Access to Nodes

 Notes: 
    -nodeTags may or may not be sequential (not sure that Gmsh guarantees anything)
    -nodeTags essentially serves as a local-index to tag lookup array
    -we can create a tag to local-index lookup (not sure if this is useful yet)
===================================#

nodeTags, nodeCoords = gmsh.model.mesh.getNodes()
nodeCoords = reshape(nodeCoords, (3, length(nodeTags))) #column major so each point is in a column.
nodeTags2Index = Dict{Int64,Int64}() #should be unsigned int but hard to read
for i in eachindex(nodeTags)
    nodeTags2Index[nodeTags[i]] = i
end

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

#==================================
 Access to Elements

 Notes: 
    -We will assume for now that we are dealing with tet meshes
    -we can create a tag to local-index lookup (not sure if this is useful yet)
    -elements contains two arrays, tags and nodes for elements
===================================#

elements = gmsh.model.mesh.getElementsByType(elementType)
elementTags = elements[1]
elementNodes = elements[2] 
numElements = length(elementTags)
elementNodes = reshape(elementNodes, (nodesPerElement, numElements)) #questions about this operation in relation to memory

elementTags2Index = Dict{Int64,Int64}() #should be unsigned int but hard to read
for i in eachindex(elementTags)
    elementTags2Index[elementTags[i]] = i
end

#we can create edges for the elements so we have access to them
gmsh.model.mesh.createEdges() #should be for all dimensions and all tags?
gmsh.model.mesh.createFaces() #should be for all dimensions and all tags?

#get nodes for edges - comes as a single list
edgeNodes = gmsh.model.mesh.getElementEdgeNodes(elementType)
@assert length(edgeNodes)==2*edgesPerElement*numElements 

#get edge tags and orientations - comes as a single list
#edge tags provides a mapping to unique edge indexes
edgeTags, edgeOrientations = gmsh.model.mesh.getEdges(edgeNodes)
maxEdgeTag = maximum(edgeTags)
numEdges = length(edgeTags)

faceNodes = gmsh.model.mesh.getElementFaceNodes(elementType,3) #3 indicates triangular faces4@assert length(faceNodes)==nodesPerFace*facesPerElement*numElements #4 is specific to tets
faceTags, faceOrientations = gmsh.model.mesh.getFaces(3, faceNodes) #3 here says that face nodes are triples
maxFaceTag = maximum(faceTags)
numFaces = length(faceTags)

#----------------------------------
# create edge-to-element maps
# each edge can belong to more than one element
# edgeTags provide a unique identifier
#----------------------------------

edges2Elements = [Vector{Int64}() for i in 1:maxEdgeTag]

for i in eachindex(edgeTags)
    #1 based indexing of tets -- edge indexes go 1, 2, 3, 4, 5, 6, 7 with corresponding elements being 1, 1, 1, 1, 1, 1, 2
    tet_index = div(i-1,6) + 1 #this is local indexing from 1 now, to change make tetTags[tet_index]
    if ~(tet_index in edges2Elements[edgeTags[i]])
        append!(edges2Elements[edgeTags[i]], tet_index)  # 6 edges per tetrahedron
    end
end

active_edge_count = 0
average_elements_per_edge = 0
for i in eachindex(edges2Elements)
    global average_elements_per_edge += length(edges2Elements[i])
    global active_edge_count += 1
end
average_elements_per_edge/=active_edge_count
println("Average elements per edge = ", average_elements_per_edge)

#----------------------------------
# create face-to-element maps
# each face can belong to more than one element
# faceTags provide a unique identifier
#----------------------------------

faces2Elements = [Vector{Int64}() for i in 1:maxFaceTag]

for i in eachindex(faceTags) 
    tet_index = div(i-1,4) + 1 #this is local indexing from 1 now, to change make tetTags[tet_index]
    if ~(tet_index in faces2Elements[faceTags[i]])
        append!(faces2Elements[faceTags[i]], tet_index)
    end
end

active_face_count = 0
average_elements_per_face = 0
for i in eachindex(faces2Elements)
    global average_elements_per_face += length(faces2Elements[i])
    global active_face_count += 1
end
average_elements_per_face/=active_face_count
println("Average elements per face = ", average_elements_per_face)

#now want to set up element neighbours and connectivity to other dimensions (for boundary conditions)

#----------------------------------
# create element-to-element (neighbour) maps
# each element has neighbours through its faces
# if a neighbour doesn't exist, set to -1
#----------------------------------

Elements2Elements = [Vector{Int64}() for i in 1:numElements]
for i in eachindex(Elements2Elements)
    Elements2Elements[i] = [-1, -1 , -1, -1]
end

for i in eachindex(faceTags) 
    tet_index = div(i-1,4) + 1 #this is local indexing from 1 now, to change make tetTags[tet_index]
    face_index = mod(i, 4) + 1 #1-based
    for j in faces2Elements[faceTags[i]] 
        if j != tet_index 
            Elements2Elements[tet_index][face_index] = j
        end
    end
end

# start playing around with information available for an element

#parv - parametric coordinates (local) for nodes
#numpriv - number of primary nodes
#numv - number of nodes
  name, dim, order, numv, parv, numpriv = gmsh.model.mesh.getElementProperties(elementType)



# try some quadrature points etc.

#these points are for integrating functions of the given order (so there are generally fewer points than required)
qpoints, qweights = gmsh.model.mesh.getIntegrationPoints(elementType, "Gauss3")
W = LinearAlgebra.Diagonal(qweights)

ele_jacobians, ele_J, ele_globalcoords = gmsh.model.mesh.getJacobian(elementTags[1],qpoints)
#gmsh/model/mesh/getElementProperties
J = LinearAlgebra.Diagonal(ele_J)
_ , basis_at_qpoints, _ = gmsh.model.mesh.getBasisFunctions(elementType, qpoints, "Lagrange1",)

#calculate a mass matrix by integrating products of basis functions over the element
# nBasis = 4
# nQpoints = length(qweights)
# B = reshape(basis_at_qpoints, (nBasis, nQpoints)) 
# M = B*(W*J)*B'


function Nrect(p1, p2, p3, p4, x, y)
    
    #note these basis functions are only for a true rectangle (not a quad)

    centroid = (p1 + p2 + p3 + p4)/4.0
    
    e1 = p2 - p1;
    e2 = p3 - p4;
    e3 = p4 - p1;
    e4 = p3 - p2;
    
    l1 = LinearAlgebra.norm(e1,2)
    l2 = LinearAlgebra.norm(e2,2)
    l3 = LinearAlgebra.norm(e3,3)
    l4 = LinearAlgebra.norm(e4,4)

    N1x = 1/l2*((centroid[2] + l2/2.0) .- Y)
    N1y = 0*N1x
    N1 = [N1x N1y]

    N2x = 1/l2*(Y .- (centroid[2] + l2/2.0))
    N2y = 0*N2x
    N2 = [N2x N2y]
    
    N3y = 1/l1*((centroid[1] + l1/2.0) .- X)
    N3x = 0*N3y
    N3 = [N3x N3y]

    N4y = 1/l1*(X .- (centroid[1] + l1/2.0))
    N4x = 0.0*N4y
    N4 = [N4x N4y]

    return N1, N2, N3, N4
    
end

function Ntri(p1, p2, p3, x, y)
    
    # so we could write this function from scratch, or we could use the gmsh api
    # and go from there.
    

    return N1, N2, N3
    
end




# x = range(start=0, stop=1, length=11)
# y = range(start=0, stop=1, length=11)
# #assume x and y are columns
# X = x' .* ones(length(y))
# Y = y .* ones(length(x))'

# Y = reshape(Y, (length(x)*length(y), 1))
# X = reshape(X, (length(x)*length(y), 1))

# p = plot(NodesX, NodesY, seriestype=:scatter, label="vertices")
# #p = quiver!(X, Y, quiver2d=(N1[:,1], N1[:,2]), label="N1")
# #p = quiver!(X, Y, quiver2d=(N2[:,1], N2[:,2]), label="N2")
# #p = quiver!(X, Y, quiver2d=(N3[:,1], N3[:,2]), label="N3")
# p = quiver!(X, Y, quiver2d=(N4[:,1], N4[:,2]), label="N4")
# display(p)

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

# condition is that u + v <= 1

ele_jacobians, ele_J, ele_global_coords = gmsh.model.mesh.getJacobian(elementTags[1],local_points)

ele_global_coords_x = ele_global_coords[1:3:end]
ele_global_coords_y = ele_global_coords[2:3:end]

local_points_u = local_points[1:3:end]
local_points_v = local_points[2:3:end]

# p2 = plot(ele_global_coords_x, ele_global_coords_y, seriestype=:scatter, label="global_points")
# display(p2)

# p3 = plot(local_points_u, local_points_v, seriestype=:scatter, label="local_points")
# display(p3)


#local_points = local_points[:] 
_ , L1, _ = gmsh.model.mesh.getBasisFunctions(elementType, local_points, "Lagrange1")
_ , L2, _ = gmsh.model.mesh.getBasisFunctions(elementType, local_points, "Lagrange2")
_ , L6, _ = gmsh.model.mesh.getBasisFunctions(elementType, local_points, "Lagrange6")

nbasis_L1 = Int64(length(L1)/(length(local_points)/3))
nbasis_L2 = Int64(length(L2)/(length(local_points)/3))
nbasis_L6 = Int64(length(L6)/(length(local_points)/3))

L2index = 1;
#L2plot = plot(local_points_u, local_points_v, L2[L2index:nbasis_L2:end], seriestype=:scatter, label="L2 basis")
#display(L2plot)

L6index = 6;
#L6plot = plot(local_points_u, local_points_v, L6[L6index:nbasis_L6:end], seriestype=:scatter, label="L6 basis")
#display(L6plot)

_ , GradL1, _ = gmsh.model.mesh.getBasisFunctions(elementType, local_points, "GradLagrange1")
_ , GradL2, _ = gmsh.model.mesh.getBasisFunctions(elementType, local_points, "GradLagrange2")
_ , GradL6, _ = gmsh.model.mesh.getBasisFunctions(elementType, local_points, "GradLagrange6")

nbasis_GradL1 = Int64(length(GradL1)/(length(local_points)/3)/3)
nbasis_GradL2 = Int64(length(GradL2)/(length(local_points)/3)/3)
nbasis_GradL6 = Int64(length(GradL6)/(length(local_points)/3)/3)

GradL1index = 1
GradL1x = GradL1[(GradL1index-1)*3 + 1:3*nbasis_GradL1:end]
GradL1y = GradL1[(GradL1index-1)*3 + 2:3*nbasis_GradL1:end]
# GradL1plot = quiver!(local_points_u, local_points_v, quiver2d=(0.1*GradL1x, 0.1*GradL1y), label="Grad Basis")
# display(GradL1plot)

ncomp_Hcurl1, Hcurl1, norient_Hcurl1 = gmsh.model.mesh.getBasisFunctions(elementType, local_points, "HcurlLegendre1", [1])
_, Hcurl2, _ = gmsh.model.mesh.getBasisFunctions(elementType, local_points, "HcurlLegendre2")
_, Hcurl6, _ = gmsh.model.mesh.getBasisFunctions(elementType, local_points, "HcurlLegendre6")

norient_Hcurl1 = 1

nbasis_Hcurl1 = Int64(length(Hcurl1)/(length(local_points) ÷ 3)/ncomp_Hcurl1/norient_Hcurl1)
nbasis_Hcurl2 = Int64(length(Hcurl2)/(length(local_points) ÷ 3)/ncomp_Hcurl1/norient_Hcurl1)
nbasis_Hcurl6 = Int64(length(Hcurl6)/(length(local_points)/3)/ncomp_Hcurl1/norient_Hcurl1)

Hcurlindex = 6
Hcurl1x = Hcurl1[(Hcurlindex-1)*3 + 1:ncomp_Hcurl1*nbasis_Hcurl1:Int64(end/norient_Hcurl1)]
Hcurl1y = Hcurl1[(Hcurlindex-1)*3 + 2:ncomp_Hcurl1*nbasis_Hcurl1:Int64(end/norient_Hcurl1)]
Hcurl1plot = quiver(local_points_u, local_points_v, quiver2d=(0.05*Hcurl1x, 0.05*Hcurl1y), label="Curl Basis")
display(Hcurl1plot)




#s = scatter(local_points_u, local_points_v, L13)
#display(s)

# outfile = "./basisprofile.txt"
# open(outfile, "w") do f
#     for i in eachindex(local_points_u)
#         println(f, local_points_u[i], " ", local_points_v[i], " ", L21[i], " ", L22[i], " ", L23[i])
#     end
# end


gmsh.finalize()
