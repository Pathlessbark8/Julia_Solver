struct HCurlElement{T} <: AbstractFiniteElement{T}
    order::Int32
    quad_weights::Vector{Float64}
    quad_points::Vector{Float64}
    functions::Array{Float64,4}
    curl_functions::Array{Float64,4}
end

# type aliases
const HCurlLine = HCurlElement{Line}
const HCurlTriangle = HCurlElement{Triangle}
const HCurlTetrahedron = HCurlElement{Tetrahedron}

function HCurlElement(dimension, order)
    @assert dimension >= 1 && dimension <= 3

    gaussOrder = 2 * (order + 1)

    # assumption: geo order doesn't matter for local function space
    # assume simplex mesh
    if dimension == 1
        typeNumber = gmsh.model.mesh.getElementType("line", 1)
        elementType = Line
    elseif dimension == 2
        typeNumber = gmsh.model.mesh.getElementType("triangle", 1)
        elementType = Triangle
    else
        typeNumber = gmsh.model.mesh.getElementType("tetrahedron", 1)
        elementType = Tetrahedron
    end

    # quad points for element type and gauss order
    quad_points, quad_weights = gmsh.model.mesh.getIntegrationPoints(typeNumber, "Gauss$gaussOrder")
    numPoints = length(quad_weights)

    # the basis functions themselves
    numComponents, local_functions, numOrientations =
        gmsh.model.mesh.getBasisFunctions(typeNumber, quad_points, "HcurlLegendre$order")
    @assert numComponents == 3
    numFunctions = length(local_functions) ÷ numComponents ÷ numOrientations ÷ numPoints
    local_functions = reshape(local_functions, (3, numFunctions, numPoints, numOrientations))

    # the curl of the basis functions
    numComponents, curl_functions, _ =
        gmsh.model.mesh.getBasisFunctions(typeNumber, quad_points, "CurlHcurlLegendre$order")
    @assert numComponents == 3
    curl_functions = reshape(curl_functions, (3, numFunctions, numPoints, numOrientations))

    return HCurlElement{elementType}(order, quad_weights, quad_points, local_functions, curl_functions)
end

# number of hierarchical hcurl functions

# for a triangle
numEdgeFunctions(local_space::HCurlTriangle) = 3 * (local_space.order + 1)
numFaceFunctions(local_space::HCurlTriangle) = local_space.order == 0 ? 0 :
                                               3 * (local_space.order - 1) + (local_space.order - 1) * (local_space.order - 2)
numBubbleFunctions(local_space::HCurlTriangle) = 0

# for a tetrahedron
numEdgeFunctions(local_space::HCurlTetrahedron) = 6 * (local_space.order + 1)
numFaceFunctions(local_space::HCurlTetrahedron) = local_space.order == 0 ? 0 :
                                                  12 * (local_space.order - 1) + 4 * (local_space.order - 1) * (local_space.order - 2)
numBubbleFunctions(local_space::HCurlTetrahedron) = local_space.order == 0 ? 0 :
                                                    (local_space.order - 1) * (local_space.order - 2) * (local_space.order - 3) ÷ 2 +
                                                    2 * (local_space.order - 2) * (local_space.order - 1)

# number of geometrical degrees of freedom
numEdges(local_space::HCurlLine) = 1
numFaces(local_space::HCurlLine) = 0
numVolumes(local_space::HCurlLine) = 0

numEdges(local_space::HCurlTriangle) = 3
numFaces(local_space::HCurlTriangle) = 1
numVolumes(local_space::HCurlTriangle) = 0

numEdges(local_space::HCurlTetrahedron) = 6
numFaces(local_space::HCurlTetrahedron) = 4
numVolumes(local_space::HCurlTetrahedron) = 1


function DOFMap(local_space::HCurlElement{T}) where {T}
    # doesn't make sense for lines yet
    @assert T == Triangle || T == Tetrahedron

    nEdgeFunctions = numEdgeFunctions(local_space)
    funcsPerEdge = nEdgeFunctions ÷ numEdges(local_space)

    nFaceFunctions = numFaceFunctions(local_space)
    funcsPerFace = nFaceFunctions ÷ numFaces(local_space)

    # numVolumes may be zero, so we should enforce than 0/0 = 0
    nBubbleFunctions = numBubbleFunctions(local_space)
    funcsPerVol = nBubbleFunctions == 0 ? 0 : nBubbleFunctions ÷ numVolumes(local_space)

    dofMapping = Tuple{Int32,Int32}[]
    for i = 1:nEdgeFunctions
        push!(dofMapping, (1, 1 + (i - 1) ÷ funcsPerEdge))
    end
    for i = 1:nFaceFunctions
        push!(dofMapping, (2, 1 + (i - 1) ÷ funcsPerFace))
    end
    for i = 1:nBubbleFunctions
        push!(dofMapping, (3, 1 + (i - 1) ÷ funcsPerVol))
    end
    return dofMapping
end

