using LinearAlgebra

#=  From gmsh source code, basis functions are ordered by increasing dimensionality, ie.
    vertex functions, edge functions, face functions, bubble functions.
    The basis function reference used for gmsh is: 
    https://www-taylorfrancis-com.uml.idm.oclc.org/books/mono/10.1201/9780203488041/higher-order-finite-element-methods-pavel-solin-ivo-dolezel-karel-segeth
    
    A bubble function is a function that has zero components on edges, but non zero components elsewhere.
    Despite having zero components at element interfaces, 
    these functions are still required to form a complete basis of H curl functions wihtin each element.

    For Hcurl triangle
    There are:
        vertexFunctions:    0
        edgeFunctions:      3 * (basisOrder + 1) 
        faceFunctions:      0 for basisOrder=0
                            3*(basisOrder-1) + (basisOrder-1)*(basisOrder-2)
        bubbleFunctions:    0
    
    NOTE: For a triangle a face function is also a bubble function


    For Hcurl tetrahedron:
        vertexFunctions:    0
        edgeFunctions:      6 * (basisOrder + 1) 
        faceFunctions:      0 for basisOrder=0
                            12*(basisOrder-1) + 4*(basisOrder-1)*(basisOrder-2)
        bubbleFunctions:    0 for basisOrder=0
                            (basisOrder-1)*(basisOrder-2)*(basisOrder-3) / 2 + 2 * (basisOrder-1)*(basisOrder-2) 
    
    To assign global degrees of freedom to these functions:
        the global number of a functions associated with an element of dimension d is related to the tag of that element of dimension d
        e.g. There are basisOrder+1 edges assocated with each global edge of an element for Hcurl functions
=#

abstract type HCurlElement end

#=
We can remove geo_order and triangleType to make struct gmsh independent
=#
struct HCurlTriangle <: HCurlElement
    geo_order::Int 
    basis_order::Int
    triangleType::Int
    quad_weights::Vector{Float64}
    quad_points::Vector{Float64}
    local_functions::Array{SVector{3,Float64}}
    curl_local_functions::Array{SVector{3,Float64}}
end

function createLocalTriangle(geoOrder::Int, basisOrder::Int)
    # order for gauss quadrature in order to integare N dot N
    gaussOrder = 2*(basisOrder+1)
    # gmsh triangle type needed to get functions
    triangleType = gmsh.model.mesh.getElementType("triangle", geoOrder)

    # quad points for triangle type and gauss order
    quad_points, quad_weights = gmsh.model.mesh.getIntegrationPoints(triangleType, "Gauss$gaussOrder")
    numPoints = length(quad_weights)

    # the basis functions themselves
    numComponents, local_functions, numOrientations = 
        gmsh.model.mesh.getBasisFunctions(triangleType, quad_points, "HcurlLegendre$basisOrder")
    numFunctions = length(local_functions) ÷ numComponents ÷ numOrientations ÷ numPoints
    local_functions = reshape(reinterpret(SVector{3, Float64}, local_functions), (numFunctions, numPoints, numOrientations))

    # the curl of the basis functions
    _, curl_functions, _ = 
        gmsh.model.mesh.getBasisFunctions(triangleType, quad_points, "CurlHcurlLegendre$basisOrder")
    curl_functions = reshape(reinterpret(SVector{3, Float64}, curl_functions), (numFunctions, numPoints, numOrientations))

    # construct and return struct
    return HCurlTriangle(geoOrder, basisOrder, triangleType, quad_weights, quad_points, local_functions, curl_functions)
end

function assembleLocalMatrices(
    triTags::Vector,
    localTriangle::HCurlTriangle)

    basisName = "HcurlLegendre"

    nPoints = length(localTriangle.quad_weights)
    nFunctions = size(localTriangle.local_functions, 1)

    M = zeros(nFunctions, nFunctions, length(triTags))
    S = zeros(nFunctions, nFunctions, length(triTags))

    globBasisFunctions = zeros(SVector{3,Float64}, nFunctions, nPoints)
    globCurlBasisFunctions = zeros(SVector{3,Float64}, nFunctions, nPoints)

    for t in eachindex(triTags)
        tag = triTags[t]

        orient = gmsh.model.mesh.getBasisFunctionsOrientationForElement(tag, basisName)
        orient += 1

        jac, det, _ = gmsh.model.mesh.getJacobian(tag, localTriangle.quad_points)
        jac_view = reinterpret(SMatrix{3,3,Float64,9}, jac)

        for f = 1:nFunctions, i = 1:nPoints
            globBasisFunctions[f, i] = jac_view[i]' \ localTriangle.local_functions[f, i, orient]
            globCurlBasisFunctions[f, i] = localTriangle.curl_local_functions[f, i, orient] / det[i]
        end

        for j = 1:nFunctions, i = 1:nFunctions
            tempM = 0
            tempS = 0
            for q = 1:nPoints
                # brackets on inner product of functions required to get exact symmetry
                tempM += localTriangle.quad_weights[q] * (globBasisFunctions[i, q]' * globBasisFunctions[j, q]) * det[q]
                tempS += localTriangle.quad_weights[q] * (globCurlBasisFunctions[i, q]' * globCurlBasisFunctions[j, q]) * det[q]
            end
            M[i,j,t] = tempM
            S[i,j,t] = tempS
        end
    end

    return M, S
end

function nEdgeFunctions(element::HCurlTriangle)
    @assert element.basis_order >=0
    return 3*(element.basis_order+1)
end

function nFaceFunctions(element::HCurlTriangle)
    @assert element.basis_order >= 0
    if element.basis_order == 0
        return 0
    else
        return 3*(element.basis_order-1) + (element.basis_order-1)*(element.basis_order-2)
    end
end

function nBubbleFunctions(element::HCurlTriangle)
    return 0
end

function getDOFMapping(element::HCurlTriangle)
    numEdgeFunctions = nEdgeFunctions(element)
    funcsPerEdge = numEdgeFunctions ÷ 3 # nEdges = 3 for triangle
    
    numFaceFunctions = nFaceFunctions(element)
    funcsPerFace = numFaceFunctions ÷ 1 # nFaces = 1 for triangle

    dofMapping = []
    for i=1:numEdgeFunctions
        push!(dofMapping, (1, 1 + (i-1) ÷ funcsPerEdge))
    end
    for i=1:numFaceFunctions
        push!(dofMapping, (2, 1 + (i-1) ÷ funcsPerFace))
    end
    return dofMapping
end

struct HCurlTetrahedron  <: HCurlElement 
    order::Int
end

function nEdgeFunctions(element::HCurlTetrahedron)
    @assert order >=0
    return 3*(order+1)
end

function nFaceFunctions(element::HCurlTetrahedron)
    @assert element.order >= 0
    if element.order == 0
        return 0
    else
        return 12*(element.order-1) + 4*(element.order-1)*(element.order-2)
    end
end

function nBubbleFunctions(element::HCurlTetrahedron)
    @assert element.order >= 0
    if element.order == 0
        return 0
    else
        return (element.order-1)*(element.order-2)*(element.order-3) ÷ 2 + 2 * (element.order-1)*(element.order-2) 
    end
end