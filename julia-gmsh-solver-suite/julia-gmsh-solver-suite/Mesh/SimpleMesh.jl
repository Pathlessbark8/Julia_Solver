using Gmsh: gmsh

"""
A simple gmsh mesh. 
Supports only homogenous mesh element types and simplex elements.
"""
struct SimpleMesh <: AbstractMesh
    properties::ElementProperties
    nodes::Nodes
    elements::Elements
    physics_map::Dict{Tuple{Int32, Int32}, Vector{Int64}} # (dim, tag) => element tags
end

mutable struct ElementEdgeProperties 
    GlobalEdgeIndex::Int64                  #global edge index
    PositiveOrientationFaceIndex::Int64     #global face index
    NegativeOrientationFaceIndex::Int64     #global face index
    PositiveOrientationLocalIndex::Int64    #local face index
    NegativeOrientationLocalIndex::Int64    #local face index
    IsBoundaryEdge::Bool                    #is this edge on the boundary
end

mutable struct SimpleMeshConnectivity
    
    edge_to_ele1D::Dict{Int64, Vector{Tuple{Int64, Int64}}}
    edge_to_ele2D::Dict{Int64, Vector{Tuple{Int64, Int64}}}
    edge_to_ele3D::Dict{Int64, Vector{Tuple{Int64, Int64}}}
    face_to_ele2D::Dict{Int64, Vector{Tuple{Int64, Int64}}}
    face_to_ele3D::Dict{Int64, Vector{Tuple{Int64, Int64}}}

    ele2D_to_ele2D::Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64}}
    ele3D_to_ele3D::Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64}}
    ele2D_to_ele1D::Dict{Tuple{Int64, Int64}, Int64}
    ele3D_to_ele2D::Dict{Tuple{Int64, Int64}, Int64}

    # Trace::Dict{Vector{Int64}, ElementEdgeProperties}
end

function setTracesForEdges(filename::String = nothing, order::Int64 = 1)

    if filename !== nothing
        gmsh.open(filename)
    end

    # access groups for physics
    gmsh.model.mesh.setOrder(order)

    # define meshing element
    elementType = gmsh.model.mesh.getElementType("triangle", order)

    #Check meshed elements properties
    elementProperties=gmsh.model.mesh.getElementProperties(elementType)

    #Get Edge Nodes in the mesh and create the edges
    edgeNodes = gmsh.model.mesh.getElementEdgeNodes(elementType)
    primaryEdgeNodes = gmsh.model.mesh.getElementEdgeNodes(elementType, -1, true)
    gmsh.model.mesh.createEdges()


    #Retrive the edge tags and Element Tags
    edgeTags,_ = gmsh.model.mesh.getEdges(primaryEdgeNodes)
    elementTags,_ = gmsh.model.mesh.getElementsByType(elementType)

    #Perform Resizing of Arrays to Edge wise ordering
    num_elements = length(elementTags)
    edgeNodesPerElement = length(edgeNodes) ÷ num_elements
    edgeNodes = reshape(edgeNodes, (edgeNodesPerElement ÷ 3, 3, num_elements))
    edgeTags = reshape(edgeTags, (3, num_elements))

    s = gmsh.model.addDiscreteEntity(1)


    #Identify and assign a unique tag to each new edge
    uniqueEdgeTags = Set()
    tagsForLines = []
    EdgeNodesForLines = []
    for i in eachindex(elementTags)
        for j in 1:3
            edge = edgeTags[j, i]
            if !in(edge, uniqueEdgeTags)
                push!(uniqueEdgeTags, edge)
                push!(tagsForLines, edge )
                append!(EdgeNodesForLines,edgeNodes[:,j,i])
            end
        end
    end

    #specficy the order and shape of the Elements being added
    elementType1D = gmsh.model.mesh.getElementType("Line", order)

    #Add the correspinding Elements to the mesh
    gmsh.model.mesh.addElementsByType(s, elementType1D, tagsForLines,EdgeNodesForLines)

    edges2Elements = Dict()
    Elements2Nbhs = Dict{Int64,Array}()
    EdgeProperties = Dict{Int64,ElementEdgeProperties}()

    # ElementWiseAverageCoordinates = Array{Float64}(undef, 3, num_elements)
    # for i in 1:1:num_elements
    #     coordinates = Array{Float64}(undef, 3, 3)
    #     for j in 1:2:6
    #         # println("primaryEdgeNodes[3*(i-1)+j] is ", primaryEdgeNodes[6*(i-1)+j])
            # coord, _= gmsh.model.mesh.getNode(primaryEdgeNodes[6*(i-1)+j])
    #         ind = floor(Int,(j-1)/2)+1
    #         coordinates[:,ind]=coord[:]
    #     end
    #     AverageCoordinates = sum(coordinates, dims=2)/3
    #     # println("Average Coordinate is ", AverageCoordinates)
    #     ElementWiseAverageCoordinates[:,i] = AverageCoordinates[:]
    # end


    for i in 1:num_elements
        for j in 1:3
            jacobian, _, global_edge_midpoint = gmsh.model.mesh.getJacobian(edgeTags[j,i], [0.0,0.0,0.0])
            # println("global_edge_midpoint is ", global_edge_midpoint)
            jacobian = reshape(jacobian, 3, 3)
            ZOrientation = jacobian[3,3]
            # println("Jacobian is ", jacobian)
            if !haskey(edges2Elements, edgeTags[j,i])
                edges2Elements[edgeTags[j,i] ] = [elementTags[i]]
                if(ZOrientation == 1)
                    EdgeProperties[edgeTags[j,i] ] = ElementEdgeProperties(edgeTags[j,i],elementTags[i],0,j,0,true)
                else
                    EdgeProperties[edgeTags[j,i] ] = ElementEdgeProperties(edgeTags[j,i],0,elementTags[i],0,j,true)
                end
            else
                push!(edges2Elements[edgeTags[j,i]],elementTags[i])
                if(ZOrientation == 1)
                    EdgeProperties[edgeTags[j,i] ].NegativeOrientationFaceIndex = elementTags[i]
                    EdgeProperties[edgeTags[j,i] ].NegativeOrientationLocalIndex = j
                else
                    EdgeProperties[edgeTags[j,i] ].PositiveOrientationFaceIndex = elementTags[i]
                    EdgeProperties[edgeTags[j,i] ].PositiveOrientationLocalIndex = j
                end
                EdgeProperties[edgeTags[j,i] ].IsBoundaryEdge = false
            end
        end
    end

    for i in 1:length(EdgeProperties)
        # if(haskey(Elements2Nbhs,EdgeProperties[i].PositiveOrientationFaceIndex))
            # append!(Elements2Nbhs[EdgeProperties[i].PositiveOrientationFaceIndex],EdgeProperties[i].NegativeOrientationFaceIndex)
        # else
        #     Elements2Nbhs[EdgeProperties[i].PositiveOrientationFaceIndex]=[EdgeProperties[i].NegativeOrientationFaceIndex]
        # end

        # if(haskey(Elements2Nbhs,EdgeProperties[i].NegativeOrientationFaceIndex))
        #     append!(Elements2Nbhs[EdgeProperties[i].NegativeOrientationFaceIndex],EdgeProperties[i].PositiveOrientationFaceIndex)
        # else
        #     Elements2Nbhs[EdgeProperties[i].NegativeOrientationFaceIndex]=[EdgeProperties[i].PositiveOrientationFaceIndex]
        # end
        if(haskey(Elements2Nbhs,EdgeProperties[i].PositiveOrientationFaceIndex))
            Elements2Nbhs[EdgeProperties[i].PositiveOrientationFaceIndex][EdgeProperties[i].PositiveOrientationLocalIndex]=EdgeProperties[i].NegativeOrientationFaceIndex
        elseif (EdgeProperties[i].PositiveOrientationFaceIndex != 0)
            Elements2Nbhs[EdgeProperties[i].PositiveOrientationFaceIndex]=zeros(Int64,3)
            Elements2Nbhs[EdgeProperties[i].PositiveOrientationFaceIndex][EdgeProperties[i].PositiveOrientationLocalIndex]=EdgeProperties[i].NegativeOrientationFaceIndex
        end
        if(haskey(Elements2Nbhs,EdgeProperties[i].NegativeOrientationFaceIndex))
            Elements2Nbhs[EdgeProperties[i].NegativeOrientationFaceIndex][EdgeProperties[i].NegativeOrientationLocalIndex]=EdgeProperties[i].PositiveOrientationFaceIndex
        elseif (EdgeProperties[i].NegativeOrientationFaceIndex != 0)
            Elements2Nbhs[EdgeProperties[i].NegativeOrientationFaceIndex]=zeros(Int64,3)
            Elements2Nbhs[EdgeProperties[i].NegativeOrientationFaceIndex][EdgeProperties[i].NegativeOrientationLocalIndex]=EdgeProperties[i].PositiveOrientationFaceIndex
        end
    end
    
    # println("Edge Properties are ", EdgeProperties)
    # gmsh.fltk.run()
    return EdgeProperties,edgeTags,elementTags,Elements2Nbhs
end

function loadMesh(filename::String = nothing)::SimpleMesh

    if filename !== nothing
        gmsh.open(filename)
    end

    # access groups for physics
    groups = gmsh.model.getPhysicalGroups()
    
    # get the elements of highest dimension only
    max_dim = maximum(getfield.(groups, 1))

    # get all elements
    properties, element_tags, element_node_tags = getAllElements(max_dim)

    num_elements = length(element_tags)
    num_nodes_per_element = length(element_node_tags) ÷ num_elements
    element_node_tags = reshape(element_node_tags, (num_nodes_per_element, num_elements))
    elements = Elements(element_tags, element_node_tags)

    # get all mesh nodes
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_tags = convert.(Int64, node_tags)
    node_coords = reshape(node_coords, (3, length(node_tags)))
    order = sortperm(node_tags) #reorder as gmsh doesn't always have tags ordered although it probably should
    nodes = Nodes(node_tags[order], node_coords[:, order])
    
    # setup physical group mappings
    physics_map = Dict{Tuple{Int32, Int32}, Vector{Int64}}()

    for (dim, tag) in groups
        element_tags_in_group, _ = getElements(dim, tag)
        physics_map[(dim,tag)] = element_tags_in_group
    end

    return SimpleMesh(properties, nodes, elements, physics_map)
end


function setupMeshConnectivity(mesh::SimpleMesh)::SimpleMeshConnectivity

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

    line_edge_tags, _ = gmsh.model.mesh.getEdges(line_edge_nodes)
    tri_edge_tags, _ = gmsh.model.mesh.getEdges(tri_edge_nodes)
    tet_edge_tags, _ = gmsh.model.mesh.getEdges(tet_edge_nodes)

    tri_face_tags, _ = gmsh.model.mesh.getFaces(3, tri_face_nodes)
    tet_face_tags, _ = gmsh.model.mesh.getFaces(3, tet_face_nodes)

    #create 
    line_tags, _ = gmsh.model.mesh.getElementsByType(line_element_type)
    tri_tags, _ = gmsh.model.mesh.getElementsByType(tri_element_type)
    tet_tags, _ = gmsh.model.mesh.getElementsByType(tet_element_type)

    #edge to ele dictionaries
    edge_to_ele1D = Dict{Int64, Vector{Tuple{Int64, Int64}}}()  #maps edge tag to pair (ele 1d tag, local edge index) (where local edge index is redundant for 1D)
    edge_to_ele2D = Dict{Int64, Vector{Tuple{Int64, Int64}}}()  #maps edge tag to pair (ele 2d tag, local edge index)
    edge_to_ele3D = Dict{Int64, Vector{Tuple{Int64, Int64}}}()  #maps edge tag to pair (ele 3d tag, local edge index)
    
    #face to ele dictionaries
    face_to_ele2D = Dict{Int64, Vector{Tuple{Int64, Int64}}}()  #maps face tag to pair (ele 2d tag, local face index) (where local face index is redundant for 2D)
    face_to_ele3D = Dict{Int64, Vector{Tuple{Int64, Int64}}}()  #maps face tag to pair (ele 3d tag, local face index)

    #ele to ele dictionaries
    ele2D_to_ele2D = Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64}}() #maps (ele 2d tag, local face index) to (ele 2d tag, local face index)
    ele3D_to_ele3D = Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64}}() #maps (ele 3d tag, local face index) to (ele 3d tag, local face index)
    ele2D_to_ele1D = Dict{Tuple{Int64, Int64}, Int64}()
    ele3D_to_ele2D = Dict{Tuple{Int64, Int64}, Int64}()

    #compute Edge to Ele 1D connectivity
    for iedge in eachindex(line_edge_tags)
        line_index = iedge
        local_edge_index = 1 #always for 1d elements
        if haskey(edge_to_ele1D,line_edge_tags[iedge]) == false
            edge_to_ele1D[line_edge_tags[iedge]] = [(line_tags[line_index], local_edge_index)] # connected on edge 1
        else
            push!(edge_to_ele1D[line_edge_tags[iedge]], (line_tags[line_index], local_edge_index)) #each line has one edge connected on edge 1 so this shouldn't happen ??
        end
    end

    #compute Edge to Ele 2D connectivity
    for iedge in eachindex(tri_edge_tags)
        tri_index = ((iedge - 1) ÷ 3) + 1
        local_edge_index = mod(iedge - 1, 3) + 1

        if haskey(edge_to_ele2D, tri_edge_tags[iedge]) == false
            edge_to_ele2D[tri_edge_tags[iedge]] = [(tri_tags[tri_index], local_edge_index)]
        else
            push!(edge_to_ele2D[tri_edge_tags[iedge]], (tri_tags[tri_index], local_edge_index)) #each tri has three edges
        end
    end

    #compute Edge to Ele 3D connectivity
    for iedge in eachindex(tet_edge_tags)
        tet_index = ((iedge - 1) ÷ 6) + 1
        local_edge_index = mod(iedge - 1, 6) + 1
        if haskey(edge_to_ele3D, tet_edge_tags[iedge]) == false
            edge_to_ele3D[tet_edge_tags[iedge]] = [(tet_tags[tet_index], local_edge_index)]
        else
            push!(edge_to_ele3D[tet_edge_tags[iedge]], (tet_tags[tet_index], local_edge_index)) #each tet has three edges)
        end
    end

    #compute Face to Ele 2D connectivity
    for iface in eachindex(tri_face_tags)
        tri_index = iface
        local_face_index = 1
        if haskey(face_to_ele2D, tri_face_tags[iface]) == false
            face_to_ele2D[tri_face_tags[iface]] = [(tri_tags[tri_index], local_face_index)]
        else
            push!(face_to_ele2D[tri_face_tags[iface]], (tri_tags[tri_index], local_face_index)) #each tri can have only one face
        end
    end

    #compute Face to Ele 3D connectivity
    for iface in eachindex(tet_face_tags)
        tet_index = ((iface - 1) ÷ 4) + 1
        local_face_index = mod(iface - 1, 4) + 1
        if haskey(face_to_ele3D, tet_face_tags[iface]) == false
            face_to_ele3D[tet_face_tags[iface]] = [(tet_tags[tet_index], local_face_index)]
        else
            push!(face_to_ele3D[tet_face_tags[iface]], (tet_tags[tet_index], local_face_index)) #each tet has four faces
        end
    end

    #compute Ele 2D to Ele 2D connectivity
    for iedge in eachindex(edge_to_ele2D)
        edge_element_connectivity = edge_to_ele2D[iedge]
        for iele in edge_element_connectivity, jele in edge_element_connectivity  #note these pull tag/local index pairs
            if iele[1] != jele[1]
                if haskey(ele2D_to_ele2D, iele) == false
                    ele2D_to_ele2D[iele] = jele #here we are adding (tag, local index) pairs
                else
                    #push!(EleToEle2D[(ielement,iface)], (jelement, jface))
                    @assert(0==1)
                end
            end
        end
    end

    #compute Ele 3D to Ele 3D connectivity
    for iface in eachindex(face_to_ele3D)
        face_element_connectivity = face_to_ele3D[iface]        
        for iele in face_element_connectivity, jele in face_element_connectivity #note these pull tag/local index pairs
            if iele[1] != jele[1]
                if haskey(ele3D_to_ele3D, iele) == false
                    ele3D_to_ele3D[iele] = jele #here we are adding (tag, local index) pairs
                else
                    #push!(EleToEle2D[(ielement,iface)], (jelement, jface))
                    @assert(0==1)
                end
            end
        end
    end

    #compute Ele 2D to Ele 1D connectivity
    for iedge in eachindex(edge_to_ele2D) #each index over a dictionary gives us the entries that exist!
        edge_ele2D_connectivity = edge_to_ele2D[iedge]
        if haskey(edge_to_ele1D, iedge) == false
            continue
        else
            elements_1d = edge_to_ele1D[iedge]
            @assert(length(elements_1d) == 1)
            element_1d = elements_1d[1]

            for iele in edge_ele2D_connectivity
                if haskey(ele2D_to_ele1D, iele) == false
                    ele2D_to_ele1D[iele] = element_1d[1] #add tag
                else
                    @assert(0==1)
                end
            end
        end
    end
    
    #compute Ele 3D to Ele 2D connectivity
    for iface in eachindex(face_to_ele3D) #each index over a dictionary gives us the entries that exist!
        face_ele3D_connectivity = face_to_ele3D[iface]
        if haskey(face_to_ele3D, iface) == false
            continue
        else
            elements_2d = face_to_ele2D[iface]
            @assert(length(elements_2d) == 1)
            element_2d = elements_2d[1]

            for iele in face_ele3D_connectivity
                if haskey(ele3D_to_ele2D, iele) == false
                    ele3D_to_ele2D[iele] = element_2d[1] #add tag
                else
                    @assert(0==1)
                end
            end
        end
    end

    connectivity = SimpleMeshConnectivity(edge_to_ele1D, edge_to_ele2D, edge_to_ele3D, face_to_ele2D, face_to_ele3D, ele2D_to_ele2D, ele3D_to_ele3D, ele2D_to_ele1D, ele3D_to_ele2D)

    return connectivity
end
