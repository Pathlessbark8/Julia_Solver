using Gmsh:gmsh
using Parameters

#=------------------------------------
Nodal Basis Struct
-------------------------------------=#
@with_kw struct EdgeBasis
    
    n_elements::Int64

    #edges
    n_unique_edges::Int64
    unique_edge_tags::Vector{Int64}
    unique_edge_nodes::Array{Int64,2}         #2 x n_edges
    element_edge_tags::Array{Int64,2}         #3 x nelements for triangles (should be made general)
    element_edge_orientations::Array{Int64,2} #3 x n_elements for triangles (should be made general)
    element_basis_orientations::Vector{Int64} #n_elements

    #quad points and geometry at those points
    n_quad_points::Int64
    quad_points::Array{Float64,2}           # 3 x n_quad_points
    quad_weights::Array{Float64,1}          # n_quad_points
    global_quad_points::Array{Float64, 3}   # 3 x n_quad_points x n_elements
    jacobians::Array{Float64,4}             # 3 x 3 x n_quad_points x nelements
    J::Array{Float64,2}                     # n_quad_points x nelements
    
    #basis functions at quad points
    n_basis::Int64                          # will be a multiple of the number of edges ?
    n_orient::Int64                         # number of possible orientations
    n_comp_basis::Int64                     # of components in the basis 
    basis::Array{Float64,4}                 # n_comp_basis x n_basis x n_quad_points x n_orient
    n_comp_curl_basis::Int64
    curl_basis::Array{Float64,4}            # n_comp x n_basis x n_quad_points x n_orient

end


#=------------------------------------
setupEdgeBasisOnTriangles() 
-------------------------------------=#
function setupEdgeBasisOnTriangles(mesh::Mesh, solution_order::Int64)
    
    basis_type = "HcurlLegendre"
    geometric_order = 1 #assumed for now
    triangle_element_type = gmsh.model.mesh.getElementType("triangle", geometric_order)
    triangle_element_type_index = findfirst(mesh.element_types .== triangle_element_type)
    @assert(triangle_element_type_index !== nothing)
    
    #create edges in the mesh for triangles
    #start by creating (dim, tag) entity pairs
    local_edge_nodes = gmsh.model.mesh.getElementEdgeNodes(triangle_element_type)

    edge_dim_tag_pairs = Vector{Tuple{Int64,Int64}}(undef,0)
    for i in eachindex(mesh.physical_groups)

        group_dim = mesh.physical_groups[i][1]
        group_tag = mesh.physical_groups[i][2]
        group_entities = gmsh.model.getEntitiesForPhysicalGroup(group_dim, group_tag)
    
        #loop over entries in physical group
        for entity in group_entities
            
            entity_element_types = gmsh.model.mesh.getElementTypes(group_dim, entity)
            
            for itype in eachindex(entity_element_types)

                #TODO: expand this to a larger comparison as needed/desired
                if entity_element_types[itype] != triangle_element_type
                    continue #only interested in these elements for this basis
                end

                push!(edge_dim_tag_pairs, (group_dim, entity))
            end
        end
    end

    #create edges
    gmsh.model.mesh.createEdges(edge_dim_tag_pairs)
    
    #get unique edges
    unique_edges = gmsh.model.mesh.getAllEdges()
    unique_edge_tags = unique_edges[1]
    unique_edge_nodes = unique_edges[2]
    n_unique_edges = length(unique_edge_tags)
    unique_edge_nodes = reshape(unique_edge_nodes, 2, n_unique_edges)

    #get element edges and orientations
    element_edge_tags, element_edge_orientations = gmsh.model.mesh.getEdges(local_edge_nodes)
    element_edge_tags = reshape(element_edge_tags, 3, mesh.elements[3].count)
    element_edge_orientations = reshape(element_edge_tags, 3, mesh.elements[3].count)
    @assert(sort(unique(unique_edge_tags)) == 1:n_unique_edges) #just check to confirm edge tags are inclusive of all values from 1 to number of unique edges. A mapping will be required if this ever breaks

    #quadrature rule for integrals of products of basis functions
    integration_order = solution_order^2 + 1 #may need an extra term here depending on how things are defined 
    quad_rule = gmsh.model.mesh.getIntegrationPoints(triangle_element_type,"Gauss$integration_order")
    n_quad_points = length(quad_rule[1]) ÷ 3
    quad_points = reshape(quad_rule[1], 3, n_quad_points) #store as 3 x n_quad_points
    quad_weights = quad_rule[2]

    global_quad_points = Array{Float64, 3}(undef, 3, n_quad_points, 0)
    jacobians_at_quad_points = Array{Float64, 4}(undef, 3, 3, n_quad_points, 0)
    determinants_at_quad_points = Array{Float64, 2}(undef, n_quad_points, 0)
    element_basis_orientations  = Vector{Int64}(undef,0)
    #we need to loop over groups and entities in the same way we did
    #when creating elements in order to be consistent with our element list
  
    basis_name = string(basis_type, solution_order)

    #loop over all physical groups in the mesh -- some of this code is common to nodal and edge basis, duplication should be removed 
    for i in eachindex(mesh.physical_groups)

        group_dim = mesh.physical_groups[i][1]
        group_tag = mesh.physical_groups[i][2]
        group_entities = gmsh.model.getEntitiesForPhysicalGroup(group_dim, group_tag)
    
        #loop over entries in physical group
        for entity in group_entities
            
            entity_element_types = gmsh.model.mesh.getElementTypes(group_dim, entity)
            
            for itype in eachindex(entity_element_types)

                #TODO: expand this to a larger comparison as needed/desired
                if entity_element_types[itype] != triangle_element_type
                    continue #only interested in these elements for this basis
                end

                entity_element_basis_orientations = gmsh.model.mesh.getBasisFunctionsOrientation(triangle_element_type, basis_name)
                element_basis_orientations = cat(element_basis_orientations, entity_element_basis_orientations, dims=1)

                entity_jacobians_at_quad_points, entity_determinants_at_quad_points, entity_global_quad_points = gmsh.model.mesh.getJacobians(entity_element_types[itype], quad_points[:], entity)
                n_entity_elements = (length(entity_global_quad_points) ÷ 3) ÷ n_quad_points
                entity_global_quad_points = reshape(entity_global_quad_points, 3, n_quad_points, n_entity_elements)
                entity_determinants_at_quad_points = reshape(entity_determinants_at_quad_points, n_quad_points, n_entity_elements) #scalar for each point and each element
                entity_jacobians_at_quad_points = reshape(entity_jacobians_at_quad_points, 3, 3, n_quad_points, n_entity_elements) #3x3 matrix for each point and each element

                global_quad_points = cat(global_quad_points, entity_global_quad_points, dims=3)
                jacobians_at_quad_points = cat(jacobians_at_quad_points, entity_jacobians_at_quad_points, dims=4)
                determinants_at_quad_points = cat(determinants_at_quad_points, entity_determinants_at_quad_points, dims=2)

            end #each type in entity
        end # group entities
    end #physical groups

    #evaluate basis functions at (local) quad points
    n_comp_basis::Int64, basis_functions, n_orient::Int64 = gmsh.model.mesh.getBasisFunctions(triangle_element_type, quad_points[:], basis_name)
    n_basis = (((length(basis_functions) ÷ n_comp_basis) ÷ n_quad_points) ÷ n_orient)
    basis_functions = reshape(basis_functions, n_comp_basis, n_basis, n_quad_points, n_orient)

    #evaluate gradient of basis function at (local) quad points
    curl_basis_name = string("Curl", basis_name)
    n_comp_curl_basis::Int64, curl_basis_functions, n_orient = gmsh.model.mesh.getBasisFunctions(triangle_element_type, quad_points[:], curl_basis_name)
    @assert n_basis == (((length(curl_basis_functions) ÷ n_comp_curl_basis) ÷ n_quad_points) ÷ n_orient)
    curl_basis_functions = reshape(curl_basis_functions, n_comp_curl_basis, n_basis, n_quad_points, n_orient)

    

    n_elements = mesh.elements[3].count
    @assert(n_elements == length(element_edge_tags[1,:]))

    #construct edge basis struct
    edge_basis = EdgeBasis(n_elements, n_unique_edges, unique_edge_tags, unique_edge_nodes, 
                            element_edge_tags, element_edge_orientations, element_basis_orientations,
                            n_quad_points, quad_points, quad_weights, global_quad_points, jacobians_at_quad_points, determinants_at_quad_points, 
                            n_basis, n_orient, n_comp_basis, basis_functions,
                            n_comp_curl_basis, curl_basis_functions)

    return edge_basis
    
end


#=------------------------------------
setupEdgeBasisOnTriangles() 
-------------------------------------=#
function setupEdgeBasisOnTetrahedrals(mesh::Mesh, solution_order::Int64)
    
    basis_type = "HcurlLegendre"
    geometric_order = 1 #assumed for now
    tetrahedron_element_type = gmsh.model.mesh.getElementType("tetrahedron", geometric_order)
    tetrahedron_element_type_index = findfirst(mesh.element_types .== tetrahedron_element_type)
    @assert(tetrahedron_element_type_index !== nothing)
    
    #create edges in the mesh for triangles
    #start by creating (dim, tag) entity pairs
    local_edge_nodes = gmsh.model.mesh.getElementEdgeNodes(tetrahedron_element_type)

    edge_dim_tag_pairs = Vector{Tuple{Int64,Int64}}(undef,0)
    for i in eachindex(mesh.physical_groups)

        group_dim = mesh.physical_groups[i][1]
        group_tag = mesh.physical_groups[i][2]
        group_entities = gmsh.model.getEntitiesForPhysicalGroup(group_dim, group_tag)
    
        #loop over entries in physical group
        for entity in group_entities
            
            entity_element_types = gmsh.model.mesh.getElementTypes(group_dim, entity)
            
            for itype in eachindex(entity_element_types)

                #TODO: expand this to a larger comparison as needed/desired
                if entity_element_types[itype] != tetrahedron_element_type
                    continue #only interested in these elements for this basis
                end

                push!(edge_dim_tag_pairs, (group_dim, entity))
            end
        end
    end

    #create edges
    gmsh.model.mesh.createEdges(edge_dim_tag_pairs)
    
    #get unique edges
    unique_edges = gmsh.model.mesh.getAllEdges()
    unique_edge_tags = unique_edges[1]
    unique_edge_nodes = unique_edges[2]
    n_unique_edges = length(unique_edge_tags)
    unique_edge_nodes = reshape(unique_edge_nodes, 2, n_unique_edges)

    #get element edges and orientations
    element_edge_tags, element_edge_orientations = gmsh.model.mesh.getEdges(local_edge_nodes)
    element_edge_tags = reshape(element_edge_tags, 6, mesh.elements[tetrahedron_element_type_index].count)
    element_edge_orientations = reshape(element_edge_tags, 6, mesh.elements[tetrahedron_element_type_index].count)
    @assert(sort(unique(unique_edge_tags)) == 1:n_unique_edges) #just check to confirm edge tags are inclusive of all values from 1 to number of unique edges. A mapping will be required if this ever breaks

    #quadrature rule for integrals of products of basis functions
    integration_order = solution_order^2 + 1 #may need an extra term here depending on how things are defined 
    quad_rule = gmsh.model.mesh.getIntegrationPoints(tetrahedron_element_type,"Gauss$integration_order")
    n_quad_points = length(quad_rule[1]) ÷ 3
    quad_points = reshape(quad_rule[1], 3, n_quad_points) #store as 3 x n_quad_points
    quad_weights = quad_rule[2]

    global_quad_points = Array{Float64, 3}(undef, 3, n_quad_points, 0)
    jacobians_at_quad_points = Array{Float64, 4}(undef, 3, 3, n_quad_points, 0)
    determinants_at_quad_points = Array{Float64, 2}(undef, n_quad_points, 0)
    element_basis_orientations  = Vector{Int64}(undef,0)
    #we need to loop over groups and entities in the same way we did
    #when creating elements in order to be consistent with our element list
  
    basis_name = string(basis_type, solution_order)

    #loop over all physical groups in the mesh -- some of this code is common to nodal and edge basis, duplication should be removed 
    for i in eachindex(mesh.physical_groups)

        group_dim = mesh.physical_groups[i][1]
        group_tag = mesh.physical_groups[i][2]
        group_entities = gmsh.model.getEntitiesForPhysicalGroup(group_dim, group_tag)
    
        #loop over entries in physical group
        for entity in group_entities
            
            entity_element_types = gmsh.model.mesh.getElementTypes(group_dim, entity)
            
            for itype in eachindex(entity_element_types)

                #TODO: expand this to a larger comparison as needed/desired
                if entity_element_types[itype] != tetrahedron_element_type
                    continue #only interested in these elements for this basis
                end

                entity_element_basis_orientations = gmsh.model.mesh.getBasisFunctionsOrientation(tetrahedron_element_type, basis_name)
                element_basis_orientations = cat(element_basis_orientations, entity_element_basis_orientations, dims=1)

                entity_jacobians_at_quad_points, entity_determinants_at_quad_points, entity_global_quad_points = gmsh.model.mesh.getJacobians(entity_element_types[itype], quad_points[:], entity)
                n_entity_elements = (length(entity_global_quad_points) ÷ 3) ÷ n_quad_points
                entity_global_quad_points = reshape(entity_global_quad_points, 3, n_quad_points, n_entity_elements)
                entity_determinants_at_quad_points = reshape(entity_determinants_at_quad_points, n_quad_points, n_entity_elements) #scalar for each point and each element
                entity_jacobians_at_quad_points = reshape(entity_jacobians_at_quad_points, 3, 3, n_quad_points, n_entity_elements) #3x3 matrix for each point and each element

                global_quad_points = cat(global_quad_points, entity_global_quad_points, dims=3)
                jacobians_at_quad_points = cat(jacobians_at_quad_points, entity_jacobians_at_quad_points, dims=4)
                determinants_at_quad_points = cat(determinants_at_quad_points, entity_determinants_at_quad_points, dims=2)

            end #each type in entity
        end # group entities
    end #physical groups

    #evaluate basis functions at (local) quad points
    n_comp_basis::Int64, basis_functions, n_orient::Int64 = gmsh.model.mesh.getBasisFunctions(tetrahedron_element_type, quad_points[:], basis_name)
    n_basis = (((length(basis_functions) ÷ n_comp_basis) ÷ n_quad_points) ÷ n_orient)
    basis_functions = reshape(basis_functions, n_comp_basis, n_basis, n_quad_points, n_orient)

    #evaluate gradient of basis function at (local) quad points
    curl_basis_name = string("Curl", basis_name)
    n_comp_curl_basis::Int64, curl_basis_functions, n_orient = gmsh.model.mesh.getBasisFunctions(tetrahedron_element_type, quad_points[:], curl_basis_name)
    @assert n_basis == (((length(curl_basis_functions) ÷ n_comp_curl_basis) ÷ n_quad_points) ÷ n_orient)
    curl_basis_functions = reshape(curl_basis_functions, n_comp_curl_basis, n_basis, n_quad_points, n_orient)

    n_elements = mesh.elements[tetrahedron_element_type].count
    @assert(n_elements == length(element_edge_tags[1,:]))

    #construct edge basis struct
    edge_basis = EdgeBasis(n_elements, n_unique_edges, unique_edge_tags, unique_edge_nodes, 
                            element_edge_tags, element_edge_orientations, element_basis_orientations,
                            n_quad_points, quad_points, quad_weights, global_quad_points, jacobians_at_quad_points, determinants_at_quad_points, 
                            n_basis, n_orient, n_comp_basis, basis_functions,
                            n_comp_curl_basis, curl_basis_functions)

    return edge_basis
    
end
