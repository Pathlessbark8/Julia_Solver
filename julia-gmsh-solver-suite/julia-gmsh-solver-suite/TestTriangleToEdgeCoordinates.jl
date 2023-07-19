using Revise

if isdefined(@__MODULE__, :EIL)
    
else
    include("EIL.jl")
end

using LinearAlgebra
using Plots
using Parameters
using NodesAndModes
using .EIL

using Gmsh:gmsh

gmsh.initialize()
gmsh.option.setNumber("Geometry.Volumes", 1)
gmsh.option.setNumber("Geometry.Surfaces", 1)
gmsh.option.setNumber("Geometry.Normals", 30)

#TM Example
mesh_filename = "./PECCircle.msh"
problem_type = EIL.EMTMProblemType()
basis_type = EIL.NodalBasisType()
solution_order = 1
integration_order = solution_order*solution_order + 1

gmsh.open(mesh_filename)

# ==================================================================================
# The purpose of this test script is to figure out:
# 1) the edge ordering and node ordering in gmsh triangles
# 2) the local triangle coordinates along each edge
# 3) how to integrate triangle basis functions along the edges of the triangle
#
# The answers to these queations are:
# 1) E1 = [v1 v2], E2 = [v2 v3], E3 = [v3 v1] where the triangle nodes are v1, v2, v3
# 2) E1 = (u, 0, 0), E2 = (u, 1-u, 0) or (1-v, v, 0), E3 = (0, v, 0)
# 3) Code is provided to integrate around the contour of a triangle and compare to a 
#    direct evaluation of permiter length.
# ==================================================================================

#basic types for this unit test
tri_type = gmsh.model.mesh.getElementType("triangle", 1)
line_type = gmsh.model.mesh.getElementType("line", 1)

#-----------------------------------------------------------------------------------
# To answer question #1 we will get the edges of each element and examine the nodes
# then compare the edge nodes to the node ordering in each element and confirm 
# the scheme used to number the edges.
# 
# The following code confirms that triangles with nodes [N1, N2, N3] have edges
# E1 = (N1, N2)
# E2 = (N2, N3)
# E3 = (N3, N1)
#-----------------------------------------------------------------------------------

# get triangle nodes
tri_element_tags, tri_element_nodes = gmsh.model.mesh.getElementsByType(tri_type) 
tri_element_nodes = convert.(Int64, tri_element_nodes)
tri_element_nodes = reshape(tri_element_nodes, 3, length(tri_element_nodes) ÷ 3)
n_tri = length(tri_element_tags)

# get triangle edge nodes
# order is [e1a1n1 e1a1n2 e1a2n1 ...] where e = element, a = edge, n = node, so the length is 6 * n_tri
tri_edge_nodes = convert.(Int64,gmsh.model.mesh.getElementEdgeNodes(tri_type))   
gmsh.model.mesh.createEdges(tri_type)
tri_element_edge_tags, tri_element_edge_orientations = gmsh.model.mesh.getEdges(tri_edge_nodes)   # 3 * number of elements
tri_element_edge_nodes = reshape(tri_edge_nodes, 2, length(tri_edge_nodes) ÷ 2)

#confirm node ordering for edges on triangles to be E1: (1,2)  E2: (2,3)  E3: (3, 1)
for iele in eachindex(tri_element_tags)
    edge_index_offset = 3*(iele-1)
    #println("iele = ", iele, " of ",length(tri_elements), " edge index offeset = ", edge_index_offset)
    edges = tri_element_edge_nodes[:, edge_index_offset .+ (1:3)]
    nodes = tri_element_nodes[:, iele]
    @assert(edges[:,1] == [nodes[1]; nodes[2]])
    @assert(edges[:,2] == [nodes[2]; nodes[3]])
    @assert(edges[:,3] == [nodes[3]; nodes[1]])
end

#-----------------------------------------------------------------------------------
# To answer question #2 we will take global points along edges and 
# produce the local coordinates in the triangle. For simplicity  we will just use
# the primary nodes corresponding to the edges. The following code confirms that 
# the coordinates (u,v,w) = (u,v,0) in 2D triangles are:
# On E1 v = 0 (coordinates are (u, 0, 0))
# On E2 1 - u - v = 0 (coordinates are (u, 1-u) or (1-v,v)) (note let t = 1 - u - v)
# On E3 u = 0 (coordinates are (0, v, 0))
#-----------------------------------------------------------------------------------

# we will need the node global coordinates
node_tags, node_coords, _ = gmsh.model.mesh.getNodesByElementType(tri_type)
node_tags = convert.(Int64, node_tags)
node_coords = reshape(node_coords, (3, length(node_coords) ÷ 3))

# verify edge coordinates for each element in the mesh
for iele in eachindex(tri_element_tags)
    
    ele_node_coords = node_coords[:, 3*(iele-1) .+ (1:3)]

    # ===== Edge 1 =======
    # node 1
    e1n1x = ele_node_coords[1,1]  
    e1n1y = ele_node_coords[2,1]
    e1n1z = ele_node_coords[3,1]
    e1n1u, e1n1v, e1n1w = gmsh.model.mesh.getLocalCoordinatesInElement(tri_element_tags[iele], e1n1x, e1n1y, e1n1z)
    e1n1t = 1 - e1n1u - e1n1v #additional barycentric coordinate t

    e1n2x = ele_node_coords[1,2]
    e1n2y = ele_node_coords[2,2]
    e1n2z = ele_node_coords[3,2]
    e1n2u, e1n2v, e1n2w = gmsh.model.mesh.getLocalCoordinatesInElement(tri_element_tags[iele], e1n2x, e1n2y, e1n2z)
    e1n2t = 1 - e1n2u - e1n2v #additional barycentric coordinate t
    
    @assert(abs(e1n2v - e1n1v - 0.0) < 1e-12) #edge 1 should have v = 0
    
    # ===== Edge 2 =======
    e2n1x = ele_node_coords[1,2]
    e2n1y = ele_node_coords[2,2]
    e2n1z = ele_node_coords[3,2]
    e2n1u, e2n1v, e2n1w = gmsh.model.mesh.getLocalCoordinatesInElement(tri_element_tags[iele], e2n1x, e2n1y, e2n1z)
    e2n1t = 1 - e2n1u - e2n1v #additional barycentric coordinate t

    e2n2x = ele_node_coords[1,3]
    e2n2y = ele_node_coords[2,3]
    e2n2z = ele_node_coords[3,3]
    e2n2u, e2n2v, e2n2w = gmsh.model.mesh.getLocalCoordinatesInElement(tri_element_tags[iele], e2n2x, e2n2y, e2n2z)
    e2n2t = 1 - e2n2u - e2n2v #additional barycentric coordinate t

    @assert(abs(e2n2t - e2n1t - 0.0) < 1e-12) #edge 2 should have 1 - u - v = 0
    
    # ===== Edge 3 =======
    e3n1x = ele_node_coords[1,3]
    e3n1y = ele_node_coords[2,3]
    e3n1z = ele_node_coords[3,3]
    e3n1u, e3n1v, e3n1w = gmsh.model.mesh.getLocalCoordinatesInElement(tri_element_tags[iele], e3n1x, e3n1y, e3n1z)
    e3n1t = 1 - e3n1u - e3n1v #additional barycentric coordinate t

    e3n2x = ele_node_coords[1,1]
    e3n2y = ele_node_coords[2,1]
    e3n2z = ele_node_coords[3,1]
    e3n2u, e3n2v, e3n2w = gmsh.model.mesh.getLocalCoordinatesInElement(tri_element_tags[iele], e3n2x, e3n2y, e3n2z)
    e3n2t = 1 - e3n2u - e3n2v #additional barycentric coordinate t

    @assert(abs(e3n2u - e3n1u - 0.0) < 1e-12) #edge 3 should have u = 0
    
end

#-----------------------------------------------------------------------------------
# To answer question #3 we need to first create line (1d) elements corresponding to 
# the boundaries of (edges) of each element. 
#-----------------------------------------------------------------------------------

new_mesh_entity = gmsh.model.addDiscreteEntity(1) #add a 1D mesh entity for the 1D edge elements

max_element_tag = gmsh.model.mesh.getMaxElementTag() # new element tags (which are unique in the mesh) must be offset from existing tags
line_element_tag_offset = max_element_tag

unique_line_tags = Set()    # we will use a set to keep track of new line elements we have already added
tags_for_lines = []         # list of new tags for line elements to be created
nodes_for_lines = []   # list of nodes defining the line elements to be created

for iedge in eachindex(tri_element_edge_tags)  # this loop over edges (non-unique) there are 3 for each element
    if !(in(tri_element_edge_tags[iedge] , unique_line_tags)) # if we haven't added this line yet
        push!(unique_line_tags, tri_element_edge_tags[iedge])   # add the line tag to the set 
        push!(tags_for_lines, tri_element_edge_tags[iedge] + line_element_tag_offset) # add the line tag to the list, offset appropriately
        append!(nodes_for_lines, tri_edge_nodes[2*(iedge-1) .+ (1:2)]) # add the node indexes in the mesh that make up the line    
    end
end

gmsh.model.mesh.addElementsByType(new_mesh_entity, line_type, tags_for_lines, nodes_for_lines)

# Now that we've created the elements, let's get them from the mesh
# Note that if we had lines in the mesh before, we will still have them in addition to the element edges

line_element_tags, line_element_nodes = gmsh.model.mesh.getElementsByType(line_type)

# We can now confirm that the nodes making up these line elements properly represent the edges
# Note that in order to do the comparison we need to sort the list of nodes as one element's edge nodes
# may oppose the edge nodes of an edge previously created 
for iedge in eachindex(tri_element_edge_tags)
    
    edge_tag = tri_element_edge_tags[iedge]
    edge_nodes_from_tri = tri_element_edge_nodes[:, iedge]
    line_tag = edge_tag + line_element_tag_offset
    _, line_nodes, _, _ = gmsh.model.mesh.getElement(line_tag)

    @assert(sort!(line_nodes) == sort!(edge_nodes_from_tri))   #the fact that we have to sort these is a problem
end

#-----------------------------------------------------------------------------------
# To finish answering question #3 we next need to generate points along edges/lines
# and calculate the 1-D jacobians along those lines in order to be able to perform
# integration. 
#
# This following code evaluates integrals around the contours. It isn't clear to me
# how to handle direction along the edges.
#-----------------------------------------------------------------------------------

# get a local integration rule on line elements for evaluating integrals
line_quad_rule = gmsh.model.mesh.getIntegrationPoints(line_type,"Gauss$integration_order")
line_n_quad_points = length(line_quad_rule[1]) ÷ 3
line_quad_points = reshape(line_quad_rule[1], 3, line_n_quad_points) #store as 3 x n_quad_points
line_quad_weights = line_quad_rule[2]

for iele in eachindex(tri_element_tags)
    #iele = 1
    ele_node_coords = node_coords[:, 3*(iele-1) .+ (1:3)]

    # ===== Edge 1 =======
    e1n1x = ele_node_coords[1,1]
    e1n1y = ele_node_coords[2,1]
    e1n1z = ele_node_coords[3,1]

    e1n2x = ele_node_coords[1,2]
    e1n2y = ele_node_coords[2,2]
    e1n2z = ele_node_coords[3,2]
    
    e1_length = sqrt((e1n2x - e1n1x)^2 + (e1n2y - e1n1y)^2 + (e1n2z - e1n1z)^2)

    tri_e1_edge_tag = tri_element_edge_tags[3*(iele-1) + 1]
    tri_e1_line_tag = tri_e1_edge_tag + line_element_tag_offset
    line_e1_jacobian, line_e1_J, line_e1_global_quad_coords = gmsh.model.mesh.getJacobian(tri_e1_line_tag, line_quad_points[:])
    e1_length_integrated = sum(line_quad_weights .* line_e1_J)

    #println("e1 direct = ", e1_length, " integrated = ", e1_length_integrated)

    # ===== Edge 2 =======
    e2n1x = ele_node_coords[1,2]
    e2n1y = ele_node_coords[2,2]
    e2n1z = ele_node_coords[3,2]

    e2n2x = ele_node_coords[1,3]
    e2n2y = ele_node_coords[2,3]
    e2n2z = ele_node_coords[3,3]

    e2_length = sqrt((e2n2x - e2n1x)^2 + (e2n2y - e2n1y)^2 + (e2n2z - e2n1z)^2)

    tri_e2_edge_tag = tri_element_edge_tags[3*(iele-1) + 2]
    tri_e2_line_tag = tri_e2_edge_tag + line_element_tag_offset
    line_e2_jacobian, line_e2_J, line_e2_global_quad_coords = gmsh.model.mesh.getJacobian(tri_e2_line_tag, line_quad_points[:])
    e2_length_integrated = sum(line_quad_weights .* line_e2_J)

    #println("e2 direct = ", e2_length, " integrated = ", e2_length_integrated)

    # ===== Edge 3 =======
    e3n1x = ele_node_coords[1,3]
    e3n1y = ele_node_coords[2,3]
    e3n1z = ele_node_coords[3,3]
    
    e3n2x = ele_node_coords[1,1]
    e3n2y = ele_node_coords[2,1]
    e3n2z = ele_node_coords[3,1]
    
    e3_length = sqrt((e3n2x - e3n1x)^2 + (e3n2y - e3n1y)^2 + (e3n2z - e3n1z)^2)

    tri_e3_edge_tag = tri_element_edge_tags[3*(iele-1) + 3]
    tri_e3_line_tag = tri_e3_edge_tag + line_element_tag_offset
    line_e3_jacobian, line_e3_J, line_e1_global_quad_coords = gmsh.model.mesh.getJacobian(tri_e3_line_tag, line_quad_points[:])
    e3_length_integrated = sum(line_quad_weights .* line_e3_J)

    #println("e3 direct = ", e3_length, " integrated = ", e3_length_integrated)

    triangle_circumference_direct = e1_length + e2_length + e3_length
    triangle_circumference_integrated = e1_length_integrated + e2_length_integrated + e3_length_integrated

    println("Direct Length = ", triangle_circumference_direct, " Integrated = ", triangle_circumference_integrated)


end


