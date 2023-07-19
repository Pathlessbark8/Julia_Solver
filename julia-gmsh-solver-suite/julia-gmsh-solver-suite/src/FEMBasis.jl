using NodesAndModes
using StructArrays

abstract type AbstractBasisType end

struct NodalBasisType <: AbstractBasisType

end

struct EdgeBasisType <: AbstractBasisType

end

abstract type AbstractDOFMap end

struct NodalDOFMap <: AbstractDOFMap

    dof_idx2unique::Array{Int64, 1}            # n_global_points
    ele_dofidx2global::Array{Int64, 2}         # n_basis x n_ele  -- likely not needed?
    ele_dofidx2unique::Array{Int64, 2}         # n_basis x n_ele

end


struct EdgeDOFMap <: AbstractDOFMap
    n_unique_edges::Int64
    unique_edge_tags::Vector{Int64}
    unique_edge_nodes::Array{Int64,2}         #2 x n_edges
    element_edge_tags::Array{Int64,2}         #3 x nelements for triangles (should be made general)
    element_edge_orientations::Array{Int64,2} #3 x n_elements for triangles (should be made general)
end

abstract type AbstractProblemType end
struct EMTMType <: AbstractProblemType end
struct EMTEType <: AbstractProblemType end
struct EM3DType <: AbstractProblemType end


#=------------------------------------
Nodal Basis Struct
-------------------------------------=#
struct LocalBasis
    
    element_type::Int64
    solution_order::Int64
    intergration_order::Int64

    # quad points and geometry at those points
    n_quad_points::Int64
    quad_points::Array{Float64,2}           # 3 x n_quad_points
    quad_weights::Array{Float64,1}          # n_quad_points
    
    # counts, orientations, components
    n_basis::Int64                          # number of basis functions (order dependent)
    n_orient::Int64                         # number of orientations for a given element
    n_comp_basis::Int64                     # number of components in the basis (1 or 3 typically)
    n_comp_grad_basis::Int64                # number of components in the gradient of basis (3)
    n_comp_curl_basis::Int64                # number of components in the curl of basis (3)

    #basis functions and their derivatives at quad points
    basis::Array{Float64,4}                 # n_comp_basis x n_basis x n_quad_points x n_orient
    grad_basis::Array{Float64,4}            # n_comp_grad_basis x n_basis x n_quad_points x n_orient
    curl_basis::Array{Float64,4}            # n_comp_curl_basis x n_basis x n_quad_points x n_orient

end


function setupLocalBasis(mesh::Mesh, element_type, basis_type::String, solution_order::Int64, basis_grad::Bool=false, basis_curl::Bool=false)

    #quadrature rule for integrals of products of basis functions
    integration_order = solution_order^2  
    quad_rule = gmsh.model.mesh.getIntegrationPoints(element_type,"Gauss$integration_order")
    n_quad_points = length(quad_rule[1]) ÷ 3
    quad_points = reshape(quad_rule[1], 3, n_quad_points) #store as 3 x n_quad_points
    quad_weights = quad_rule[2]

    #evaluate basis functions at (local) quad points
    basis_name = string(basis_type, solution_order)
    n_comp_basis::Int64, basis_functions, n_orient::Int64 = gmsh.model.mesh.getBasisFunctions(element_type, quad_points[:], basis_name)
    n_basis = (((length(basis_functions) ÷ n_comp_basis) ÷ n_quad_points) ÷ n_orient)
    basis = reshape(basis_functions, n_comp_basis, n_basis, n_quad_points, n_orient)

    #evaluate gradient of basis function at (local) quad points
    if basis_grad
        grad_basis_name = string("Grad", basis_name)
        n_comp_grad_basis::Int64, grad_basis_functions, n_orient = gmsh.model.mesh.getBasisFunctions(element_type, quad_points[:], grad_basis_name)
        @assert n_basis == (((length(grad_basis_functions) ÷ n_comp_grad_basis) ÷ n_quad_points) ÷ n_orient)
        grad_basis = reshape(grad_basis_functions, n_comp_grad_basis, n_basis, n_quad_points, n_orient)
    else
        n_comp_grad_basis = 0
        grad_basis = Array{Float64,4}(undef,0,0,0,0)
    end

    if basis_curl
        curl_basis_name = string("Curl", basis_name)
        n_comp_curl_basis::Int64, curl_basis_functions, n_orient = gmsh.model.mesh.getBasisFunctions(element_type, quad_points[:], curl_basis_name)
        @assert n_basis == (((length(curl_basis_functions) ÷ n_comp_curl_basis) ÷ n_quad_points) ÷ n_orient)
        curl_basis = reshape(curl_basis_functions, n_comp_curl_basis, n_basis, n_quad_points, n_orient)
    else
        n_comp_curl_basis = 0
        curl_basis = Array{Float64,4}(undef,0,0,0,0)
    end

    local_basis = LocalBasis(element_type, solution_order, integration_order, n_quad_points, quad_points, quad_weights, 
                            n_basis, n_orient, n_comp_basis, 
                            n_comp_grad_basis, n_comp_curl_basis, 
                            basis, grad_basis, curl_basis)

    return local_basis

end


function setupLocalBasis(mesh::Mesh, solution_order::Int64, problem_type::EMTMType)
    geometric_order = 1
    triangle_type = gmsh.model.mesh.getElementType("triangle", geometric_order)
    local_basis = setupLocalBasis(mesh, triangle_type, "Lagrange", solution_order, true, false)
    return local_basis            
end

function setupLocalBasis(mesh::Mesh, solution_order::Int64, problem_type::EMTEType)
    geometric_order = 1    
    triangle_type = gmsh.model.mesh.getElementType("triangle", geometric_order)
    local_basis = setupLocalBasis(mesh, triangle_type, "HcurlLegendre", solution_order, false, true)
    return local_basis
end

function setupLocalBasis(mesh::Mesh, solution_order::Int64, problem_type::EM3DType)
    geometric_order = 1
    tetrahedron_type = gmsh.model.mesh.getElementType("tetrahedron", geometric_order)
    local_basis = setupLocalBasis(mesh, tetrahedron_type, "HcurlLegendre", solution_order, false, true)
    return local_basis
end


#=------------------------------------
setupDOFMap -- Nodal Basis Version 
-------------------------------------=#
function setupDOFMap(element_type, solution_order, basis_type::NodalBasisType)
    
    #local nodes associated with nodal basis
    if element_type == gmsh.model.mesh.getElementType("triangle",1)
        r, s = nodes(Tri(), solution_order) #call to NodesAndModes
        t = 0 .* r
    elseif element_type == gmsh.model.mesh.getElementType("tetrahedron",1)
        r, s, t = nodes(Tet(), solution_order) #call to NodesAndModes
    else
        assert(0==1)
    end
    
    #map from [-1 1] to expected gmsh interval [0 1]
    r = (r .+ 1)./2.0                     
    s = (s .+ 1)./2.0
    t = (t .+ 1)./2.0
    n_local_points = length(r)

    #create local points array
    local_points = Array{Float64, 1}(undef, 3*n_local_points)
    local_points[1:3:end] = r
    local_points[2:3:end] = s
    local_points[3:3:end] = t
    local_points = reshape(local_points, 3, n_local_points) #storing as 3 x n_local_points

    #global points at local points for each element - this comes from gmsh one element at a time
    global_points = Array{Float64,2}(undef,3,0)
    _, _, global_points = gmsh.model.mesh.getJacobians(element_type, local_points[:])
    n_global_points = length(global_points) ÷ 3
    global_points = reshape(global_points, 3, n_global_points) #storing as 3 x n_global_points
   
    #create a list of unique nodes from global points as follows:
    # round the points for comparison
    rounded_global_points = round.(global_points, digits=10)
    rounded_global_points = StructArray((rounded_global_points[1,:], rounded_global_points[2,:], rounded_global_points[3,:]))

    # sort the points
    sort_indexes = sortperm(rounded_global_points)
    sorted_rounded_global_points = rounded_global_points[sort_indexes]

    # find unique points (hence the rounding) and create a idx2unique map - just march through points until
    # some change is found, then you have a new point
    node_idx2unique = Array{Int64, 1}(undef, n_global_points)
    next_unique_index = 1
    for i in eachindex(sorted_rounded_global_points)
        if (i == 1)
            node_idx2unique[sort_indexes[i]] = next_unique_index
        else
            if (sorted_rounded_global_points[i] != sorted_rounded_global_points[i-1])
                next_unique_index+=1
            end
            node_idx2unique[sort_indexes[i]] = next_unique_index
        end
    end
    n_dof = next_unique_index #the number of degrees of freedom in this nodal basis

    #definitely a better way to do this with reshaping 1:n_global_points

    n_elements = n_global_points ÷ n_local_points
    n_basis = n_local_points
    
    ele_idx2unique = Array{Int64, 2}(undef, n_basis, n_elements)
    ele_idx2global = Array{Int64, 2}(undef, n_basis, n_elements)

    for iele in range(1,n_elements)
        global_ids = (1:n_basis) .+ (iele-1)*n_basis
        unique_ids = [node_idx2unique[global_id] for global_id in global_ids]
        ele_idx2global[:,iele] = global_ids
        ele_idx2unique[:,iele] = unique_ids
    end

    nodal_dof_map = NodalDOFMap(node_idx2unique, ele_idx2global, ele_idx2unique)

    return nodal_dof_map
    
end

#=------------------------------------
setupDOFMap -- Edge Basis Version 
-------------------------------------=#
function setupDOFMap(element_type, solution_order, basis_type::EdgeBasisType)
    
    basis_type = "HcurlLegendre"
    geometric_order = 1 #assumed for now
    
    #local edges based on element_type
    local_edge_nodes = gmsh.model.mesh.getElementEdgeNodes(element_type)
    
    #create edges - assume for entire mesh
    gmsh.model.mesh.createEdges()
    
    #get unique edges
    unique_edges = gmsh.model.mesh.getAllEdges()
    unique_edge_tags = unique_edges[1]
    unique_edge_nodes = unique_edges[2]
    n_unique_edges = length(unique_edge_tags)
    unique_edge_nodes = reshape(unique_edge_nodes, 2, n_unique_edges)

    #get element edges and orientations
    element_edge_tags, element_edge_orientations = gmsh.model.mesh.getEdges(local_edge_nodes)

    if element_type == gmsh.model.mesh.getElementType("triangle", geometric_order)
        edges_per_ele = 3
    elseif element_type == gmsh.model.mesh.getElementType("tetrahedron",geometric_order)
        edges_per_ele = 6
    else
        @assert(0==1)
    end

    n_elements = length(element_edge_tags) ÷ edges_per_ele

    element_edge_tags = reshape(element_edge_tags, 3, n_elements)
    element_edge_orientations = reshape(element_edge_tags, 3, n_elements)
    @assert(sort(unique(unique_edge_tags)) == 1:n_unique_edges) #just check to confirm edge tags are inclusive of all values from 1 to number of unique edges. A mapping will be required if this ever breaks

    
    @assert(n_elements == length(element_edge_tags[1,:]))

    #construct edge basis struct
    edge_dof_map = EdgeDOFMap(n_unique_edges, unique_edge_tags, unique_edge_nodes, 
                              element_edge_tags, element_edge_orientations)

    return edge_dof_map
    
end








# global_quad_points::Array{Float64, 3}   # 3 x n_quad_points x n_elements
# # geometric info at quad points
# jacobians::Array{Float64,4}             # 3 x 3 x n_quad_points x nelements
# J::Array{Float64,2}                     # n_quad_points x nelements

# #element to basis indexing helpers
# ele_idx2global::Array{Int64, 2}         # n_basis x n_ele
# ele_idx2unique::Array{Int64, 2}         # n_basis x n_ele

# #indexing helpers for elements (assumes constant order)
# n_elements = n_global_points ÷ n_local_points
# @assert(n_elements == mesh.elements[3].count) #assuming 2D elements -- this assert may have to go
# ele_idx2unique = Array{Int64, 2}(undef, n_basis, n_elements)
# ele_idx2global = Array{Int64, 2}(undef, n_basis, n_elements)

# #global points, jabocians and determinants at quad points -- allocation
# jacobians_at_quad_points, determinants_at_quad_points, global_quad_points = gmsh.model.mesh.getJacobians(element_type, quad_points[:])
# n_elements = (length(global_quad_points) ÷ 3) ÷ n_quad_points
# global_quad_points = reshape(global_quad_points, 3, n_quad_points, n_elements)
# determinants_at_quad_points = reshape(determinants_at_quad_points, n_quad_points, n_elements) #scalar for each point and each element
# jacobians_at_quad_points = reshape(jacobians_at_quad_points, 3, 3, n_quad_points, n_elements) #3x3 matrix for each point and each element

# #it may be of interest to map unique nodal points to gmsh points for plotting purposes
# node_idx2gmsh = Array{Int64, 1}(undef, n_global_points)
# for i in eachindex(node_idx2unique) #total number of nodes
#     element_id = ((i-1) ÷ 3) + 1 #we know we have 3 nodes per element for a 1st order solution
#     node = global_points[:,i]
#     found = false
#     for iele_node = 1:3
#         ele_gmsh_node_tag = mesh.elements[3].nodes[iele_node,element_id] #assuming this ordering stands up
#         ele_gmsh_node_idx = mesh.nodes.tag2idx[ele_gmsh_node_tag]
#         mesh_node = mesh.nodes.coords[:,ele_gmsh_node_idx]
#         if (norm(node - mesh_node) < 1e-12)
#             found = true
#             node_idx2gmsh[i] = ele_gmsh_node_idx
#             break
#         end
#     end
#     @assert(found == true)
# end


