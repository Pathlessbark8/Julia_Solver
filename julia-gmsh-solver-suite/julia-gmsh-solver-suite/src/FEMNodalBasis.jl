using Gmsh:gmsh
using Parameters

#=------------------------------------
Nodal Basis Struct
-------------------------------------=#
@with_kw struct NodalBasis
    
    n_elements::Int64

    #local nodal basis points
    n_local_points::Int64
    local_points::Array{Float64,2}          # 3 x n_local_points
    n_global_points::Int64      
    global_points::Array{Float64,2}         # 3 x n_global_points
    n_unique_global_points::Int64
    node_idx2unique::Array{Int64,1}         # map to unique index for each global point
    node_idx2gmsh::Array{Int64, 1}          # map to gmsh index for each global point (only works for 1st order)
    
    #quad points and geometry at those points
    n_quad_points::Int64
    quad_points::Array{Float64,2}           # 3 x n_quad_points
    quad_weights::Array{Float64,1}          # n_quad_points
    global_quad_points::Array{Float64, 3}   # 3 x n_quad_points x n_elements
    jacobians::Array{Float64,4}             # 3 x 3 x n_quad_points x nelements
    J::Array{Float64,2}                     # n_quad_points x nelements
    
    #basis functions at quad points
    n_basis::Int64
    n_orient::Int64
    n_comp_basis::Int64    
    basis::Array{Float64,4}                 # n_comp_basis x n_basis x n_quad_points x n_orient
    n_comp_grad_basis::Int64
    grad_basis::Array{Float64,4}            # n_comp x n_basis x n_quad_points x n_orient

    #element to basis indexing helpers
    ele_idx2global::Array{Int64, 2}         # n_basis x n_ele
    ele_idx2unique::Array{Int64, 2}         # n_basis x n_ele
end


#=------------------------------------
setupNodalBasis() 
-------------------------------------=#
function setupNodalBasisOnTriangles(mesh::Mesh, solution_order::Int64)
    
    basis_type = "Lagrange"
    geometric_order = 1 #assumed for now
    triangle_element_type = gmsh.model.mesh.getElementType("triangle", geometric_order)
    triangle_element_type_index = findfirst(mesh.element_types .== triangle_element_type)
    @assert(triangle_element_type_index !== nothing)
    
    #local nodes associated with nodal basis
    r, s = nodes(Tri(), solution_order) #call to NodesAndModes
    n_local_points = length(r)

    #map from [-1 1] to expected gmsh interval [0 1]
    r = (r .+ 1)./2.0                     
    s = (s .+ 1)./2.0
    
    #create local points array
    local_points = Array{Float64, 1}(undef, 3*n_local_points)
    local_points[1:3:end] = r
    local_points[2:3:end] = s
    local_points[3:3:end] = 0 .* r #despite being 3D - so we need a 2D tag somewhere to automate this
    local_points = reshape(local_points, 3, n_local_points) #storing as 3 x n_local_points

    #quadrature rule for integrals of products of basis functions
    integration_order = solution_order^2  
    quad_rule = gmsh.model.mesh.getIntegrationPoints(triangle_element_type,"Gauss$integration_order")
    n_quad_points = length(quad_rule[1]) ÷ 3
    quad_points = reshape(quad_rule[1], 3, n_quad_points) #store as 3 x n_quad_points
    quad_weights = quad_rule[2]

    #global points at local points -- allocation
    global_points = Array{Float64,2}(undef,3,0)

    #global points, jabocians and determinants at quad points -- allocation
    global_quad_points = Array{Float64, 3}(undef, 3, n_quad_points, 0)
    jacobians_at_quad_points = Array{Float64, 4}(undef, 3, 3, n_quad_points, 0)
    determinants_at_quad_points = Array{Float64, 2}(undef, n_quad_points, 0)

    #we need to loop over groups and entities in the same way we did
    #when creating elements in order to be consistent with our element list
  
    #loop over all physical groups in the mesh
    for iter in mesh.iterators

        group_dim = iter.group_dim
        group_tag = iter.group_tag
        entity = iter.entity
        element_type = iter.element_type
        element_type_index = iter.element_type_index

        #TODO: expand this to a larger comparison as needed/desired
        if element_type != triangle_element_type
            continue #only interested in these elements for this basis
        end

        #global coords associated with nodal basis
        _, _, entity_global_points = gmsh.model.mesh.getJacobians(element_type, local_points[:], entity)
        n_entity_global_points = length(entity_global_points) ÷ 3
        entity_global_points = reshape(entity_global_points, 3, n_entity_global_points) #storing as 3 x n_global_points
        global_points = [global_points entity_global_points] #concatenate to global list

        entity_jacobians_at_quad_points, entity_determinants_at_quad_points, entity_global_quad_points = gmsh.model.mesh.getJacobians(element_type, quad_points[:], entity)
        n_entity_elements = (length(entity_global_quad_points) ÷ 3) ÷ n_quad_points
        entity_global_quad_points = reshape(entity_global_quad_points, 3, n_quad_points, n_entity_elements)
        entity_determinants_at_quad_points = reshape(entity_determinants_at_quad_points, n_quad_points, n_entity_elements) #scalar for each point and each element
        entity_jacobians_at_quad_points = reshape(entity_jacobians_at_quad_points, 3, 3, n_quad_points, n_entity_elements) #3x3 matrix for each point and each element

        global_quad_points = cat(global_quad_points, entity_global_quad_points, dims=3)
        jacobians_at_quad_points = cat(jacobians_at_quad_points, entity_jacobians_at_quad_points, dims=4)
        determinants_at_quad_points = cat(determinants_at_quad_points, entity_determinants_at_quad_points, dims=2)

    end #iterator

    _, n_global_points = size(global_points)

    #create a list of unique nodes from global points as follows:
    # round the points for comparison
    rounded_global_points = round.(global_points, digits=10)
    rounded_global_points = StructArray((rounded_global_points[1,:], rounded_global_points[2,:], rounded_global_points[3,:]))

    # sort the points
    sort_indexes = sortperm(rounded_global_points)
    sorted_rounded_global_points = rounded_global_points[sort_indexes]

    # find unique points (hence the rounding) and create a idx2unique map
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
    n_unique_global_points = next_unique_index

    # The global points list, 3 x n_global_points, is constructed as follows: for a given element all points corresponding to the local basis points.
    # If we want to use gmsh plotting on nodes (for say first-order elements) we need to map those nodal basis points to the actual gmsh nodes representing
    # the geometry. The easiest way to do this is to loop over the elements, the points in the elements, and check which vertex they are equal to (so it seems).
    # Above, we have a list of unique points that is literally numbered from 1 to # unique nodes. The mapping to these unique nodes is over a similar sequential mapping
    # in the global nodes list (completely uncoupled from gmsh numbering).

    node_idx2gmsh = Array{Int64, 1}(undef, n_global_points)
    for i in eachindex(node_idx2unique) #total number of nodes
        element_id = ((i-1) ÷ 3) + 1 #we know we have 3 nodes per element for a 1st order solution
        node = global_points[:,i]
        found = false
        for iele_node = 1:3
            ele_gmsh_node_tag = mesh.elements[3].nodes[iele_node,element_id] #assuming this ordering stands up
            ele_gmsh_node_idx = mesh.nodes.tag2idx[ele_gmsh_node_tag]
            mesh_node = mesh.nodes.coords[:,ele_gmsh_node_idx]
            if (norm(node - mesh_node) < 1e-12)
                found = true
                node_idx2gmsh[i] = ele_gmsh_node_idx
                break
            end
        end
        @assert(found == true)
    end


    #evaluate basis functions at (local) quad points
    basis_name = string(basis_type, solution_order)
    n_comp_basis::Int64, basis_functions, n_orient::Int64 = gmsh.model.mesh.getBasisFunctions(triangle_element_type, quad_points[:], basis_name)
    n_basis = (((length(basis_functions) ÷ n_comp_basis) ÷ n_quad_points) ÷ n_orient)
    basis_functions = reshape(basis_functions, n_comp_basis, n_basis, n_quad_points, n_orient)

    #evaluate gradient of basis function at (local) quad points
    grad_basis_name = string("Grad", basis_name)
    n_comp_grad_basis::Int64, grad_basis_functions, n_orient = gmsh.model.mesh.getBasisFunctions(triangle_element_type, quad_points[:], grad_basis_name)
    @assert n_basis == (((length(grad_basis_functions) ÷ n_comp_grad_basis) ÷ n_quad_points) ÷ n_orient)
    grad_basis_functions = reshape(grad_basis_functions, n_comp_grad_basis, n_basis, n_quad_points, n_orient)

    #indexing helpers for elements (assumes constant order)
    n_elements = n_global_points ÷ n_local_points
    @assert(n_elements == mesh.elements[3].count) #assuming 2D elements -- this assert may have to go
    ele_idx2unique = Array{Int64, 2}(undef, n_basis, n_elements)
    ele_idx2global = Array{Int64, 2}(undef, n_basis, n_elements)
    
    #definitely a better way to do this with reshaping 1:n_global_points
    for iele in range(1,n_elements)
        global_ids = (1:n_basis) .+ (iele-1)*n_basis
        unique_ids = [node_idx2unique[global_id] for global_id in global_ids]
        ele_idx2global[:,iele] = global_ids
        ele_idx2unique[:,iele] = unique_ids
    end

    #construct nodal basis struct
    nodal_basis = NodalBasis(n_elements, n_local_points, local_points, n_global_points, global_points,
                        n_unique_global_points, node_idx2unique, node_idx2gmsh,
                        n_quad_points, quad_points, quad_weights, global_quad_points, jacobians_at_quad_points, determinants_at_quad_points, 
                        n_basis, n_orient, n_comp_basis, basis_functions,
                        n_comp_grad_basis, grad_basis_functions, ele_idx2global, ele_idx2unique)

    return nodal_basis
    
    
end
