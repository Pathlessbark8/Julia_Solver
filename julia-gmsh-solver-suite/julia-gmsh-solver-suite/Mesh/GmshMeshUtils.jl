
struct ElementProperties
    type::Int32
    element_name::String
    dim::Int32
    order::Int32
    num_nodes::Int32
    local_node_coord
    num_primary_nodes::Int32 
end



function Properties(element_type)
    name, dim, order, nNodes, localCoord, nPrimary = gmsh.model.mesh.getElementProperties(element_type)
    return ElementProperties(element_type, name, dim, order, nNodes, localCoord, nPrimary)
end

function getElements(dim, physical_group)
    entities = gmsh.model.getEntitiesForPhysicalGroup(dim, physical_group)

    elementTags = Vector{UInt64}(undef, 0)
    elementNodeTags = Vector{UInt64}(undef, 0)

    for e in entities
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, e)

        # we expect only a single element type
        @assert size(elemTypes, 1) == 1

        append!(elementTags, elemTags[1])
        append!(elementNodeTags, elemNodeTags[1])
    end

    nElem = size(elementTags, 1)
    nNodes = size(elementNodeTags, 1)
    nodesPerElem = nNodes ÷ nElem

    elementNodeTags = reshape(elementNodeTags, (nodesPerElem, nElem))
    return elementTags, elementNodeTags
end

function createTraces(dim)

end

"""
Get all elements in mesh of dimension dim

Currently assumes only one element type.
"""
function getAllElements(dim)
    entities = gmsh.model.getEntities(dim)

    println("entities: ", entities)
    elementTags = Int64[]
    elementNodeTags = Int64[]
    types = Int64[]

    for (d, t) in entities
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(d, t)
        
        println("elemTypes: ", d,t)
        @assert size(elemTypes, 1) == 1
        append!(types, elemTypes[1])

        append!(elementTags, elemTags[1])
        append!(elementNodeTags, elemNodeTags[1])
    end

    @assert all(types[1] .== types[2:end])


    return Properties(types[1]), elementTags, elementNodeTags
end

function createEdges(element_type::Integer, numElements::Integer)
    # create edge numbering
    gmsh.model.mesh.createEdges()

    # get edges for each element
    edgeNodes = gmsh.model.mesh.getElementEdgeNodes(element_type, -1, true)
    edgeTags, _ = gmsh.model.mesh.getEdges(edgeNodes)
    edgeTags = convert.(Int64, edgeTags)

    edgesPerElement = length(edgeTags) ÷ numElements
    edgeTags = reshape(edgeTags, (edgesPerElement, numElements))

    # test edge tags
    test = sort(unique(edgeTags))
    if ~all(test .== test[begin]:test[end])
        uniqueEdges, uniqueEdgeNodes = gmsh.model.mesh.getAllEdges()
        uniqueEdgeNodes = reshape(uniqueEdgeNodes, (2, length(uniqueEdges)))

        uniqueEdges = convert.(Int64, uniqueEdges)
        missingEdges = setdiff(uniqueEdges, test)
        println("missing edges: ")
        display(sort(missingEdge))
        @assert length(missingEdges) == 0
    end
    return edgeTags
end

"""
Assumes triangular faces.
"""
function createFaces(element_type::Integer, numElements::Integer)
    # create face numbering
    gmsh.model.mesh.createFaces()

    # get faces for each element
    faceNodes = gmsh.model.mesh.getElementFaceNodes(element_type, 3, -1, true)
    faceTags, _ = gmsh.model.mesh.getFaces(3, faceNodes)
    faceTags = convert.(Int64, faceTags)

    facesPerElement = length(faceTags) ÷ numElements
    faceTags = reshape(faceTags, (facesPerElement, numElements))

    test = sort(unique(faceTags))
    if ~all(test .== test[begin]:test[end])
        uniqueFaces, _ = gmsh.model.mesh.getAllFaces(3)
        uniqueFaces = convert.(Int64, uniqueFaces)
        missingFaces = setdiff(uniqueFaces, test)
        println("missing faces:")
        display(sort(missingFaces))
        @assert length(missingFaces) == 0
    end
    return faceTags
end

function createPartitions(N)
    gmsh.option.setNumber("Mesh.PartitionCreateTopology", 1)
    gmsh.option.setNumber("Mesh.PartitionCreateGhostCells", 0)
    gmsh.option.setNumber("Mesh.PartitionCreatePhysicals", 1)
    gmsh.model.mesh.partition(N)
end

function getPartitions()
    N = gmsh.model.getNumberOfPartitions()
    partitions = [Vector{Tuple}(undef, 0) for n = 1:N]

    entities = gmsh.model.getEntities()

    for e in entities
        parts = gmsh.model.getPartitions(e[1], e[2])
        if length(parts) > 0
            for p in parts
                push!(partitions[p], e)
            end
        end
    end
    return partitions
end

