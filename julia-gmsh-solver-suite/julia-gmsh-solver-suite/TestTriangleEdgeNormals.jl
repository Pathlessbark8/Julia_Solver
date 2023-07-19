using Revise

if isdefined(@__MODULE__, :EIL)
    
else
    include("EIL.jl")
end

using LinearAlgebra
using Plots
using Parameters
using NodesAndModes
using StaticArrays
using .EIL

using Gmsh:gmsh

gmsh.finalize()
setupGmsh() 

# ==================================================================================
# The purpose of this test script is to figure out:
# 1) how to compute the outer normal to each edge of a triangular element from the jacobians
# ==================================================================================
# The outer normal to the edge (surface) of a triangle is equal to the edge tangent 
# crossed onto the element surface normal. 
# 
# We can get the edge tangent by making 1d elements corresponding to edges, 
# computing the jacobian matrix for the element, and taking a1 (first column) 
# of that jacobian matrix.
#
# We can get the surface normal by taking a3 from the triangle jacobian.
#
# It turns out that orienting the 1d element properly requires confirming its orientation
# with the original triangle (the element edge nodes don't necessarily reflect the original triangle)
#
#-----------------------------------------------------------------------------------


#basic types for this unit test
tri_type = gmsh.model.mesh.getElementType("triangle", 1)
line_type = gmsh.model.mesh.getElementType("line", 1)

#-----------------------------------------------------------------------------------
# Start by creating a circular mesh with a normal zaxis = [alpha, beta, gamma]
#-----------------------------------------------------------------------------------

zaxis = [3, 2, 2]

circle_id = gmsh.model.occ.addCircle(0,0,0, 1, -1, 0, 2*pi, zaxis);
lc1 = 0.3
gmsh.model.occ.mesh.setSize(gmsh.model.occ.getEntities(0), lc1)
gmsh.model.occ.synchronize()
curve_loop_id = gmsh.model.geo.addCurveLoop([circle_id])
surface_id = gmsh.model.geo.addPlaneSurface([curve_loop_id]);

gmsh.model.geo.synchronize();
gmsh.model.addPhysicalGroup(2, [surface_id], 3000, "freespace")
circle_boundary = gmsh.model.getBoundary((2,circle_id))
circle_boundary = abs.(getfield.(circle_boundary, 2))
gmsh.model.addPhysicalGroup(1, circle_boundary, 2000, "pec")
gmsh.model.mesh.generate(2)
gmsh.write("PECCircleAngle.msh")
# gmsh.fltk.run()

#-----------------------------------------------------------------------------------
# Pull the triangle tags and nodes
#-----------------------------------------------------------------------------------

tri_element_tags, tri_element_nodes = gmsh.model.mesh.getElementsByType(tri_type) 
tri_element_nodes = convert.(Int64, tri_element_nodes)
tri_element_nodes = reshape(tri_element_nodes, 3, length(tri_element_nodes) ÷ 3)
n_tri = length(tri_element_tags)

#-----------------------------------------------------------------------------------
# Get the triangle edge nodes
# order is [e1a1n1 e1a1n2 e1a2n1 ...] where e = element, a = edge, n = node, so the length is 6 * n_tri
# 
# *** note that these edge nodes aren't necessarily in the order corresponding to the way the nodes 
# show up in the element
# -----------------------------------------------------------------------------------
tri_edge_nodes = convert.(Int64,gmsh.model.mesh.getElementEdgeNodes(tri_type))   
gmsh.model.mesh.createEdges(tri_type)
tri_element_edge_tags, tri_element_edge_orientations = gmsh.model.mesh.getEdges(tri_edge_nodes)   # 3 * number of elements
tri_element_edge_nodes = reshape(tri_edge_nodes, 2, length(tri_edge_nodes) ÷ 2)

#-----------------------------------------------------------------------------------
# Create elements for the triangle edges so that we can compute jacobians along them
#
# *** note the resulting elements add another layer of "reordering" as multiple edges
# for elements correspond to the same physical line and we only want to make the line
# into an element once. if the node ordering is off, the line will be misaligned to the
# edge for some element(s).
#-----------------------------------------------------------------------------------
new_mesh_entity = gmsh.model.addDiscreteEntity(1) #add a 1D mesh entity for the 1D edge elements

max_element_tag = gmsh.model.mesh.getMaxElementTag() # new element tags (which are unique in the mesh) must be offset from existing tags
line_element_tag_offset = max_element_tag

unique_line_tags = Set()    # we will use a set to keep track of new line elements we have already added
tags_for_lines = []         # list of new tags for line elements to be created
nodes_for_lines = []   # list of nodes defining the line elements to be created

# here we create gmsh line elements corresponding to the edges
for iedge in eachindex(tri_element_edge_tags)  # this loop over edges (non-unique) there are 3 for each element
    if !(in(tri_element_edge_tags[iedge] , unique_line_tags)) # if we haven't added this line yet
        push!(unique_line_tags, tri_element_edge_tags[iedge])   # add the line tag to the set 
        push!(tags_for_lines, tri_element_edge_tags[iedge] + line_element_tag_offset) # add the line tag to the list, offset appropriately
        append!(nodes_for_lines, tri_edge_nodes[2*(iedge-1) .+ (1:2)]) # add the node indexes in the mesh that make up the line    
    end
end
gmsh.model.mesh.addElementsByType(new_mesh_entity, line_type, tags_for_lines, nodes_for_lines)

#-----------------------------------------------------------------------------------
# Get the newly created line elements from the mesh
#-----------------------------------------------------------------------------------

line_element_tags, line_element_nodes = gmsh.model.mesh.getElementsByType(line_type) 
line_element_nodes = convert.(Int64, line_element_nodes)
line_element_nodes = reshape(line_element_nodes, 2, length(line_element_nodes) ÷ 2)
n_line = length(line_element_tags)

#-----------------------------------------------------------------------------------
# Get the node tags and coordinates
# - here the tags are required for confirming ordering
# - the coordinates are used as a secondary check on edge orientation but aren't 
# strictly required to compute the edge normals
#-----------------------------------------------------------------------------------

line_node_tags, line_node_coords, _ = gmsh.model.mesh.getNodesByElementType(line_type)
line_node_tags = convert.(Int64, line_node_tags)
line_node_coords = reshape(line_node_coords, 3, length(line_node_coords) ÷ 3)

tri_node_tags, tri_node_coords, _ = gmsh.model.mesh.getNodesByElementType(tri_type)
tri_node_tags = convert.(Int64, tri_node_tags)
tri_node_coords = reshape(tri_node_coords, 3, length(tri_node_coords) ÷ 3)


#-----------------------------------------------------------------------------------
# Next we will compute the outward edge normal to each edge of each triangle in the
# mesh. This is evaluated at the barycentre of each edge but could easily be 
# extended to evaluation at an arbitrary set of quadrature points along the edge
#
# Of particular interest is the fact that the edge_orientation being +/- does not
# affect the geometry of the problem. We suspect this is only related to the definition
# of growing local coordinate along the edge for the purpose of integration.
#-----------------------------------------------------------------------------------

verbose = true #if you want output set true

#loop over triangles
for tri_index in 1:3 #eachindex(tri_element_tags)
 
    tri_tag = tri_element_tags[tri_index]
    
    if verbose println("Triangle ", tri_index, " with tag ", tri_tag) end

    ele_node_coords = tri_node_coords[:, 3*(tri_index-1) .+ (1:3)]

    #loop over edges
    for tri_edge_index in [1,2,3]

        edge_tag = tri_element_edge_tags[3*(tri_index-1) + tri_edge_index]
        edge_orientation = tri_element_edge_orientations[3*(tri_index-1) + tri_edge_index]
        
        if verbose println("\tEdge # ", tri_edge_index, " with tag ", edge_tag, " and orientation ", edge_orientation) end

        # these node tags come from the getElementEdgeNodes function call
        edge_nodes_from_tri_element_edge_nodes = tri_element_edge_nodes[:, 3*(tri_index-1) + tri_edge_index]
        
        if verbose println("\t\tEdge nodes from triangle edge_nodes = ", edge_nodes_from_tri_element_edge_nodes) end

        # get the nodes of the edge directly from the element (as opposed to from the element edges)
        # including the coordinates for checking tangents
        if tri_edge_index == 1
            edge_p1 = ele_node_coords[:,1]
            edge_p2 = ele_node_coords[:,2]
            edge_nodes_from_tri = tri_element_nodes[1:2,tri_index]
        elseif tri_edge_index == 2
            edge_p1 = ele_node_coords[:,2]
            edge_p2 = ele_node_coords[:,3]
            edge_nodes_from_tri = tri_element_nodes[2:3,tri_index]
        elseif tri_edge_index == 3
            edge_p1 = ele_node_coords[:,3]
            edge_p2 = ele_node_coords[:,1]
            edge_nodes_from_tri = [tri_element_nodes[3,tri_index], tri_element_nodes[1,tri_index]]
        end
        edge_tangent = edge_p2 - edge_p1
        edge_tangent_normalized = edge_tangent/norm(edge_tangent)

        if verbose
            println("\t\tEdge nodes from tri element = ", edge_nodes_from_tri) 
            println("\t\tEdge p1 = ", edge_p1)
            println("\t\tEdge p2 = ", edge_p2)
            println("\t\tEdge tangent = ", edge_tangent)
            println("\t\tEdge tangent normalized = ", edge_tangent_normalized)
        end

        # we will grab the nodes from the line as well as this is where the jacobian will be coming from
        # and we should check alignment of the line with the triangle edge
        line_tag = edge_tag + line_element_tag_offset
        _, line_nodes, _, _ = gmsh.model.mesh.getElement(line_tag)
        line_nodes = convert.(Int64, line_nodes)

        if verbose
            println("\t\tCorresponding Line has tag ", line_tag)
            println("\t\tNodes from line = ", line_nodes)
        end
        
        # check alignment of line nodes and tri edge nodes
        if line_nodes == edge_nodes_from_tri_element_edge_nodes
            line_to_edge_multiplier = 1
        elseif sort(line_nodes) == sort(edge_nodes_from_tri_element_edge_nodes)
            line_to_edge_multiplier = -1
        else
            @assert(0==1)
        end

        # check alignment of tri edge nodes and triangle nodes
        if edge_nodes_from_tri_element_edge_nodes == edge_nodes_from_tri
            edge_to_tri_multiplier = 1
        elseif sort(edge_nodes_from_tri_element_edge_nodes) == sort(edge_nodes_from_tri)
            edge_to_tri_multiplier = -1
        else
            @assert(0==1)
        end

        # check alignment of line nodes and triangle nodes -- only this should be needed 
        # to get the outward normal adjusted properly
        if line_nodes == edge_nodes_from_tri
            line_to_tri_multiplier = 1
        elseif sort(line_nodes) == sort(edge_nodes_from_tri)
            line_to_tri_multiplier = -1
        else
            @assert(0==1)
        end
        
        if verbose
            println("\t\tline_to_edge_multiplier = ", line_to_edge_multiplier)
            println("\t\tedge_to_tri_multiplier = ", edge_to_tri_multiplier)
            println("\t\tline_to_tri_multiplier = ", line_to_tri_multiplier)
        end
        
        # here, if we wanted to compute the outward normal directly we would need to determine how the line nodes 
        # orient with respect to the edge of the element. however, we may just correct the normal to be outwards
        # after-the-fact

        line_coord = 0
        line_local_point = [line_coord, 0, 0]

        tri_coord = (1+line_coord)/2
        if tri_edge_index == 1
            tri_local_point = [tri_coord, 0, 0]
        elseif tri_edge_index == 2
            tri_local_point = [tri_coord, 1-tri_coord, 0]
        elseif tri_edge_index == 3
            tri_local_point = [0,tri_coord, 0]
        else
            @assert(0==1)
        end

        # get the geometric info from the line element
        line_jacobian, line_determinant, line_global_point = gmsh.model.mesh.getJacobian(line_tag, line_local_point[:])
        line_jacobian = reshape(line_jacobian, 3, 3)

        line_a1 = line_jacobian[:,1]
        line_a2 = line_jacobian[:,2]
        line_a3 = line_jacobian[:,3]

        # get the geometric info from the triangle
        tri_jacobian, tri_determinant, tri_global_point = gmsh.model.mesh.getJacobian(tri_tag, tri_local_point[:])
        tri_jacobian = reshape(tri_jacobian, 3, 3)

        tri_a1 = tri_jacobian[:,1]
        tri_a2 = tri_jacobian[:,2]        
        tri_a3 = tri_jacobian[:,3]

        # confirm we have the same point from both
        @assert(norm(line_global_point - tri_global_point) < 1e-12)

        # grab the surface normal - comes normalised
        surface_normal = tri_a3 

        # grab the edge vector in the direction of the triangle edge definition
        edge_tangent_from_line = line_to_tri_multiplier .* line_a1
        edge_tangent_from_line = edge_tangent_from_line/norm(edge_tangent_from_line) # normalize just for comparison with element

        if verbose println("\t\tedge_tangent_from_line = ", edge_tangent_from_line) end
    
        edge_normal = cross(edge_tangent_from_line, surface_normal)
        edge_normal = edge_normal / norm(edge_normal)
        if verbose println("\t\tedge_normal = ", edge_normal) end

        # Confirmation # 1: take dot product with the vector from the center of the element to the midpoint of the edge
        # which should be positive
        _, _, edge_centre_point = gmsh.model.mesh.getJacobian(line_tag, [0,0,0])
        _, _, tri_centre_point = gmsh.model.mesh.getJacobian(tri_tag, [1/3, 1/3, 0])
        tri_to_edge = edge_centre_point - tri_centre_point
        edge_normal_direction = dot(edge_normal, tri_to_edge)
        if verbose println("\t\tedge_normal_direction = ", edge_normal_direction) end
        @assert(edge_normal_direction > 0)

        #second, take the cross product with the edge normal to the line 
        surface_normal_from_edge_normal = cross(edge_normal, edge_tangent_from_line)
        surface_normal_from_edge_normal = surface_normal_from_edge_normal/norm(surface_normal_from_edge_normal)
        @assert(norm(surface_normal_from_edge_normal - surface_normal) < 1e-12)
        
        if verbose 
            println("\t\tsurface normal = ", surface_normal)
            println("\t\tsurface_normal_from_edge_normal = ", surface_normal_from_edge_normal)
        end
        

   end
end
 

