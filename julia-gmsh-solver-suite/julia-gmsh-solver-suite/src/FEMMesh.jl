using Gmsh:gmsh
#using Parameters

#=------------------------------------
Nodes Struct (Mesh Vertices)
-------------------------------------=#
struct Nodes
    count::Int64              
    tags::Array{Int64}          
    coords::Array{Float64,2}
    tag2idx::Dict{Int64,Int64}
end

#=------------------------------------
Elements Struct (Mesh Elements)
-------------------------------------=#
struct Elements
    count::Int64
    tags::Array{Int64,1}
    types::Array{Int64,1}
    nodes_per_element::Int64
    nodes::Array{Int64,2}
    tag2idx::Dict{Int64,Int64}
    physics_tags::Array{Int64,1}
end

#=------------------------------------
Mesh Struct
-------------------------------------=#
struct Mesh
    filename::String
    element_types::Vector{Int64}
    nodes::Nodes
    elements::Vector{Elements}
end

struct MeshConnectivity
    
    EdgeToEle1D::Dict{Int64, Vector{Tuple{Int64, Int64}}}
    EdgeToEle2D::Dict{Int64, Vector{Tuple{Int64, Int64}}}
    EdgeToEle3D::Dict{Int64, Vector{Tuple{Int64, Int64}}}
    FaceToEle2D::Dict{Int64, Vector{Tuple{Int64, Int64}}}
    FaceToEle3D::Dict{Int64, Vector{Tuple{Int64, Int64}}}

    Ele2DToEle2D::Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64}}
    Ele3DToEle3D::Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64}}
    Ele2DToEle1D::Dict{Tuple{Int64, Int64}, Int64}
    Ele3DToEle2D::Dict{Tuple{Int64, Int64}, Int64}

end

#=------------------------------------
setupGmsh() -- setup basis options we prefer
-------------------------------------=#
function setupGmsh()
    gmsh.initialize()
    gmsh.option.setNumber("Geometry.Volumes", 1)
    gmsh.option.setNumber("Geometry.Surfaces", 1)
    gmsh.option.setNumber("Geometry.Normals", 30)
end

#=------------------------------------
readMesh() -- reads a gmsh mesh - some options hardcoded
-------------------------------------=#
function readMesh(filename::String)
    
    # assumed mesh type
    geometric_order = 1

    # create list of types of elements we want to read
    element_types_to_read = Vector{Int64}(undef,0)
    push!(element_types_to_read, gmsh.model.mesh.getElementType("point", geometric_order))
    push!(element_types_to_read, gmsh.model.mesh.getElementType("line", geometric_order))
    push!(element_types_to_read, gmsh.model.mesh.getElementType("triangle", geometric_order))
    push!(element_types_to_read, gmsh.model.mesh.getElementType("tetrahedron", geometric_order))
    n_element_types_to_read = length(element_types_to_read)

    # determine number of nodes in each element we want to read
    element_types_to_read_num_nodes = Vector{Int64}(undef,n_element_types_to_read)
    for itype in eachindex(element_types_to_read)
        ele_type = element_types_to_read[itype]
        _, _, _, num_nodes, _ ,_ = gmsh.model.mesh.getElementProperties(ele_type)
        element_types_to_read_num_nodes[itype] = num_nodes        
    end
    
    # file read
    gmsh.open(filename)
    
    # set up nodes - this accounts for all geometric nodes (vertices) in the mesh
    node_tags, node_coords = gmsh.model.mesh.getNodes()
    node_count = length(node_tags)
    node_coords = reshape(node_coords, (3, node_count))
    node_tag2idx = Dict{Int64,Int64}()
    for i in eachindex(node_tags)
        node_tag2idx[node_tags[i]] = i
    end
    nodes = Nodes(node_count, node_tags, node_coords, node_tag2idx)

    # create a dictionary that maps element tags to physical groups
    gmsh_physical_groups = gmsh.model.getPhysicalGroups()
    ele_tag_to_physics_tag = Dict{Int64, Int64}()

    for group in gmsh_physical_groups
        
        group_dim = group[1]
        group_tag = group[2]
        group_entities = gmsh.model.getEntitiesForPhysicalGroup(group_dim, group_tag)
        
        #loop over entries in physical groups
        for entity in group_entities
            
            entity_element_types = gmsh.model.mesh.getElementTypes(group_dim, entity)

            for itype in eachindex(entity_element_types)

                #check to ensure the type in the entity is something we want to read
                element_type_index = findfirst(element_types_to_read .== entity_element_types[itype])
                if element_type_index === nothing
                    continue
                end

                element_type = element_types_to_read[element_type_index]
                entity_elements = gmsh.model.mesh.getElementsByType(element_type, entity)
                entity_element_tags = entity_elements[1]
                
                for tag in entity_element_tags

                    if haskey(ele_tag_to_physics_tag, tag) == false
                        ele_tag_to_physics_tag[tag] = group_tag
                    else
                        assert(0==1)
                    end
                end
            end
        end
    end

    #create data arrays for element structs
    element_counts = Vector{Int64}(undef, n_element_types_to_read)
    element_counts = fill(0,n_element_types_to_read)
    element_tags = [Vector{Int64}(undef, 0) for i in 1:n_element_types_to_read]
    element_types = [Vector{Int64}(undef, 0) for i in 1:n_element_types_to_read]
    element_physical_tags = [Vector{Int64}(undef, 0) for i in 1:n_element_types_to_read]
    element_nodes = [Array{Int64,2}(undef,element_types_to_read_num_nodes[itype],0) for itype in 1:n_element_types_to_read]
        
    elements = Vector{Elements}(undef,0)

    # next, read all elements of the type we are interested in and create elements storage
    
    for itype in eachindex(element_types_to_read)

        element_type = element_types_to_read[itype]
        type_elements = gmsh.model.mesh.getElementsByType(element_type)
        
        element_tags[itype] = type_elements[1]
        type_elements_count = length(element_tags[itype])
        element_counts[itype] = type_elements_count

        if element_counts[itype] > 0
            nodes_per_element = length(type_elements[2]) ÷ type_elements_count
            element_nodes[itype] = reshape(type_elements[2], nodes_per_element, type_elements_count)
            element_types[itype] = fill(element_type, type_elements_count)
        
            element_physical_tags[itype] = [ele_tag_to_physics_tag[element_tags[itype][itag]] for itag in eachindex(element_tags[itype])]
        
            @assert(element_counts[itype] == length(element_tags[itype]))
            @assert(element_counts[itype] == length(element_types[itype]))
            @assert(element_counts[itype] == length(element_physical_tags[itype]))
        end

        #enable direct lookup of element by tag
        ele_tag2idx = Dict{Int64,Int64}() 
        for i in eachindex(element_tags[itype])
            ele_tag2idx[element_tags[itype][i]] = i
        end

        #create element struct for this type
        elements_of_type = Elements(element_counts[itype],element_tags[itype], 
                                    element_types[itype], 
                                    element_types_to_read_num_nodes[itype], 
                                    element_nodes[itype], ele_tag2idx, 
                                    element_physical_tags[itype])
        push!(elements, elements_of_type)
    end

    #create mesh struct
    mesh = Mesh(filename, element_types_to_read, nodes, elements)
    return mesh
end


function setupMeshConnectivity(mesh::Mesh)

    gmsh.model.mesh.createEdges() #creates all unique edges
    gmsh.model.mesh.createFaces() #creates all unique faces
    
    line_element_type = gmsh.model.mesh.getElementType("line", 1)
    tri_element_type = gmsh.model.mesh.getElementType("triangle", 1)
    tet_element_type = gmsh.model.mesh.getElementType("tetrahedron", 1)

    line_edge_nodes = gmsh.model.mesh.getElementEdgeNodes(line_element_type)
    tri_edge_nodes = gmsh.model.mesh.getElementEdgeNodes(tri_element_type)
    tet_edge_nodes = gmsh.model.mesh.getElementEdgeNodes(tet_element_type)

    tri_face_nodes = gmsh.model.mesh.getElementFaceNodes(tri_element_type, 3)
    tet_face_nodes = gmsh.model.mesh.getElementFaceNodes(tet_element_type, 3)

    line_edge_tags, line_edge_orientations = gmsh.model.mesh.getEdges(line_edge_nodes)
    tri_edge_tags, tri_edge_orientations = gmsh.model.mesh.getEdges(tri_edge_nodes)
    tet_edge_tags, tet_edge_orientations = gmsh.model.mesh.getEdges(tet_edge_nodes)

    tri_face_tags, tri_face_orientations = gmsh.model.mesh.getFaces(3, tri_face_nodes)
    tet_face_tags, tet_face_orientations = gmsh.model.mesh.getFaces(3, tet_face_nodes)

    #create 
    line_tags, line_node_tags = gmsh.model.mesh.getElementsByType(line_element_type)
    tri_tags, tri_node_tags = gmsh.model.mesh.getElementsByType(tri_element_type)
    tet_tags, tet_node_tags = gmsh.model.mesh.getElementsByType(tet_element_type)

    #edge to ele dictionaries
    EdgeToEle1D = Dict{Int64, Vector{Tuple{Int64, Int64}}}()
    EdgeToEle2D = Dict{Int64, Vector{Tuple{Int64, Int64}}}()
    EdgeToEle3D = Dict{Int64, Vector{Tuple{Int64, Int64}}}()
    
    #face to ele dictionaries
    FaceToEle2D = Dict{Int64, Vector{Tuple{Int64, Int64}}}()
    FaceToEle3D = Dict{Int64, Vector{Tuple{Int64, Int64}}}()

    #ele to ele dictionaries
    Ele2DToEle2D = Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64}}()
    Ele3DToEle3D = Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64}}()
    Ele2DToEle1D = Dict{Tuple{Int64, Int64}, Int64}()
    Ele3DToEle2D = Dict{Tuple{Int64, Int64}, Int64}()

    #compute Edge to Ele 1D connectivity
    for iedge in eachindex(line_edge_tags)
        line_index = iedge
        local_edge_index = 1
        if haskey(EdgeToEle1D,line_edge_tags[iedge]) == false
            EdgeToEle1D[line_edge_tags[iedge]] = [(line_tags[line_index], local_edge_index)] # connected on edge 1
        else
            push!(EdgeToEle1D[line_edge_tags[iedge]], (line_tags[line_index], local_edge_index)) #each line has one edge connected on edge 1
        end
    end

    #compute Edge to Ele 2D connectivity
    for iedge in eachindex(tri_edge_tags)
        tri_index = ((iedge - 1) ÷ 3) + 1
        local_edge_index = mod(iedge - 1, 3) + 1

        if haskey(EdgeToEle2D, tri_edge_tags[iedge]) == false
            EdgeToEle2D[tri_edge_tags[iedge]] = [(tri_tags[tri_index], local_edge_index)]
        else
            push!(EdgeToEle2D[tri_edge_tags[iedge]], (tri_tags[tri_index], local_edge_index)) #each tri has three edges
        end
    end

    #compute Edge to Ele 3D connectivity
    for iedge in eachindex(tet_edge_tags)
        tet_index = ((iedge - 1) ÷ 6) + 1
        local_edge_index = mod(iedge - 1, 6) + 1
        if haskey(EdgeToEle3D, tet_edge_tags[iedge]) == false
            EdgeToEle3D[tet_edge_tags[iedge]] = [(tet_tags[tet_index], local_edge_index)]
        else
            push!(EdgeToEle3D[tet_edge_tags[iedge]], (tet_tags[tet_index], local_edge_index)) #each tet has three edges)
        end
    end

    #compute Face to Ele 2D connectivity
    for iface in eachindex(tri_face_tags)
        tri_index = iface
        local_face_index = 1
        if haskey(FaceToEle2D, tri_face_tags[iface]) == false
            FaceToEle2D[tri_face_tags[iface]] = [(tri_tags[tri_index], local_face_index)]
        else
            push!(FaceToEle2D[tri_face_tags[iface]], (tri_tags[tri_index], local_face_index)) #each tri can have only one face
        end
    end

    #compute Fdge to Ele 3D connectivity
    for iface in eachindex(tet_face_tags)
        tet_index = ((iface - 1) ÷ 4) + 1
        local_face_index = mod(iface - 1, 4) + 1
        if haskey(FaceToEle3D, tet_face_tags[iface]) == false
            FaceToEle3D[tet_face_tags[iface]] = [(tet_tags[tet_index], local_face_index)]
        else
            push!(FaceToEle3D[tet_face_tags[iface]], (tet_tags[tet_index], local_face_index)) #each tet has four faces
        end
    end

    #compute Ele 2D to Ele 2D connectivity
    for iedge in eachindex(EdgeToEle2D)
        edge_element_connectivity = EdgeToEle2D[iedge]
        for iele in edge_element_connectivity, jele in edge_element_connectivity 
            if iele[1] != jele[1]
                if haskey(Ele2DToEle2D, iele) == false
                    Ele2DToEle2D[iele] = jele
                else
                    #push!(EleToEle2D[(ielement,iface)], (jelement, jface))
                    @assert(0==1)
                end
            end
        end
    end

    #compute Ele 3D to Ele 3D connectivity
    for iface in eachindex(FaceToEle3D)
        face_element_connectivity = FaceToEle3D[iface]        
        for iele in face_element_connectivity, jele in face_element_connectivity 
            if iele[1] != jele[1]
                if haskey(Ele3DToEle3D, iele) == false
                    Ele3DToEle3D[iele] = jele
                else
                    #push!(EleToEle2D[(ielement,iface)], (jelement, jface))
                    @assert(0==1)
                end
            end
        end
    end

    #compute Ele 2D to Ele 1D connectivity
    for iedge in eachindex(EdgeToEle2D) #each index over a dictionary gives us the entries that exist!
        edge_Ele2D_connectivity = EdgeToEle2D[iedge]
        if haskey(EdgeToEle1D, iedge) == false
            continue
        else
            elements_1d = EdgeToEle1D[iedge]
            @assert(length(elements_1d) == 1)
            element_1d = elements_1d[1]

            for iele in edge_Ele2D_connectivity
                if haskey(Ele2DToEle1D, iele) == false
                    Ele2DToEle1D[iele] = element_1d[1]
                else
                    @assert(0==1)
                end
            end
        end
    end
    
    #compute Ele 3D to Ele 2D connectivity
    for iface in eachindex(FaceToEle3D) #each index over a dictionary gives us the entries that exist!
        face_Ele3D_connectivity = FaceToEle3D[iface]
        if haskey(FaceToEle2D, iface) == false
            continue
        else
            elements_2d = FaceToEle2D[iface]
            @assert(length(elements_2d) == 1)
            element_2d = elements_2d[1]

            for iele in face_Ele3D_connectivity
                if haskey(Ele3DToEle2D, iele) == false
                    Ele3DToEle2D[iele] = element_2d[1]
                else
                    @assert(0==1)
                end
            end
        end
    end

    connectivity = MeshConnectivity(EdgeToEle1D, EdgeToEle2D, EdgeToEle3D, FaceToEle2D, FaceToEle3D, Ele2DToEle2D, Ele3DToEle3D, Ele2DToEle1D, Ele3DToEle2D)

    return connectivity
end




