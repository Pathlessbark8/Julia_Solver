import Gmsh: gmsh
using CairoMakie # plotting
# using GLMakie
using LinearAlgebra
using SparseArrays
using StaticArrays
using Arpack
using KrylovKit
using Statistics
using WriteVTK
using Printf

include("gmsh_mesh_utils.jl")
include("gmsh_basis_utils.jl")
include("visualization.jl")
include("analytic.jl")

const ϵ0::Float64 = 8.854e-12
const μ0::Float64 = 4 * π * 1e-7
const c0::Float64 = 3.0e8

ω::Float64 = 300e6 / 2 / π
k2::Float64 = ω^2 * ϵ0 * μ0

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
                tempM += localTriangle.quad_weights[q] * globBasisFunctions[i, q]' * globBasisFunctions[j, q] * det[q]
                tempS += localTriangle.quad_weights[q] * globCurlBasisFunctions[i, q]' * globCurlBasisFunctions[j, q] * det[q]
            end
            if abs(tempM) > 1e-15
                M[i, j, t] = tempM
            end
            if abs(tempS) > 1e-15
                S[i, j, t] = tempS
            end
        end
    end

    return M, S
end

"""
Assemble global sparse matrix into default julia sparse format (CSC).

Inputs:

    triTags - vector of gmsh tags for every triangle in global domain (no hybrid mesh support)
    edgeTags - edgeTags returned from setupEdges()
    pecTags - vector of gmsh tag for edges on pec boundaries
    basisOrder - order of the basis function used (Alternate approach could be to pass in dofMapping)
    M - local element mass matrices for every element, in same order as triTags
    S - local element curl matrices for every element, in same order as triTags,
    
Returns:

    Kglob - assembled sparse matrix of S-k^2M
    Mglob - assembled sparse matrix of k^2M
    globalPecTags - potentially high order pec tags used to create a reduced matrix/apply pec boundaries
    globalFunctionIdxs - Table of global indexes of functions for each triangle, useful for plotting

    Note: my general approach to functions in this program is to ignore temporary allocations as
    I assume global memory usage will be dominated by the lu factorization.
"""
function assembleGlobalMatrices(
    triTags::Vector,
    edgeTags::Matrix,
    localTriangle::HCurlTriangle,
    M::Array,
    S::Array;
    pecTags::Vector=[],
    abcTags::Vector=[]
)
    # we need to know the max edge id in order to assign face dof to after edge dof
    maxEdgeId = maximum(edgeTags)

    # get the associated local degrees of freedom for each basis function
    # used to map local edges to global edges
    dofMapping = getDOFMapping(localTriangle)
    localEdgeIdxs = [t for (d, t) in dofMapping if d == 1]
    nEdgeFuncs = length(localEdgeIdxs)
    funcsPerEdge = nEdgeFuncs ÷ 3 # for a triangle

    localFaceIdxs = [t for (d, t) in dofMapping if d == 2]
    nFaceFuncs = length(localFaceIdxs)
    funcsPerFace = nFaceFuncs ÷ 1 # 1 for a triangle

    # create global DOF table
    globalFunctionIdxs = Matrix{Int32}(undef, nEdgeFuncs + nFaceFuncs, length(triTags))
    edgeOffset = funcsPerEdge * maxEdgeId
    for t in eachindex(triTags)
        # for each edge function
        for ifunc = 1:nEdgeFuncs
            funcOrder = 1 + (ifunc - 1) % funcsPerEdge
            zeroOrderEdgeIndex = edgeTags[localEdgeIdxs[ifunc], t]
            highOrderEdgeIndex = funcsPerEdge * (zeroOrderEdgeIndex - 1) + funcOrder
            globalFunctionIdxs[ifunc, t] = highOrderEdgeIndex
        end
        # for each face function
        for ifunc = 1:nFaceFuncs
            funcOrder = 1 + (ifunc - 1) % funcsPerFace
            highOrderFaceIndex = edgeOffset + funcsPerFace * (t - 1) + funcOrder
            globalFunctionIdxs[nEdgeFuncs+ifunc, t] = highOrderFaceIndex
        end
    end

    # pec global edge tags
    # pec boundaries not associated with face basis functions for triangles
    globalPecTags = Int32[]
    for t in eachindex(pecTags)
        # for each edge function
        for ifunc = 1:funcsPerEdge
            highOrderEdgeIndex = funcsPerEdge * (pecTags[t] - 1) + ifunc
            push!(globalPecTags, highOrderEdgeIndex)
        end
    end

    # create COO lists
    I = Int64[]
    J = Int64[]
    VK = ComplexF64[]
    VM = ComplexF64[]
    # we could add k2[t] here for inhomogenous permittivity
    for t in eachindex(triTags)
        for jfunc = 1:nEdgeFuncs+nFaceFuncs
            for ifunc = 1:nEdgeFuncs+nFaceFuncs
                vk = S[ifunc, jfunc, t] - k2 * M[ifunc, jfunc, t]
                vm = k2 * M[ifunc, jfunc, t]

                push!(I, globalFunctionIdxs[ifunc, t])
                push!(J, globalFunctionIdxs[jfunc, t])
                push!(VK, vk)
                push!(VM, vm)
            end
        end
    end

    Kglob = sparse(I, J, VK)
    Mglob = sparse(I, J, VM)

    # contribution from ABC bc's
    Kabc = nothing
    if length(abcTags) > 0
        println("Integrating abc's here")

        # edge to element
        edgeToEleMap = Dict{Int32,Int32}()

        # find the associated element for each abc tag
        # assume only at most one element is associated
        sort!(abcTags)
        for itri in eachindex(triTags)
            for t in edgeTags[:, itri]
                e = searchsortedfirst(abcTags, t)
                if e <= length(abcTags)
                    edgeToEleMap[t] = itri
                end
            end
        end

        # for each abc line, we need to integrate
        # ∫jk(n x ϕ)⋅(n x E) dl
        basisOrder = localTriangle.basis_order
        gaussOrder = 2 * (basisOrder + 1)

        lineType = gmsh.model.mesh.getElementType("line", localTriangle.geo_order)
        coord, weights = gmsh.model.mesh.getIntegrationPoints(lineType, "Gauss$gaussOrder")
        numPoints = length(weights)

        # get functions on edge
        numComponents, local_functions, numOrientations = gmsh.model.mesh.getBasisFunctions(lineType, coord, "HcurlLegendre$basisOrder")
        numFunctions = length(local_functions) ÷ numComponents ÷ numOrientations ÷ numPoints
        local_functions = reshape(reinterpret(SVector{3,Float64}, local_functions), (numFunctions, numPoints, numOrientations))

        # curl local functions
        _, curl_functions, _ = gmsh.model.mesh.getBasisFunctions(lineType, coord, "CurlHcurlLegendre$basisOrder")
        curl_functions = reshape(reinterpret(SVector{3,Float64}, curl_functions), (numFunctions, numPoints, numOrientations))

        # convert local to global and integrate
        k = ω * sqrt(ϵ0 * μ0)
        Iabc = Int32[]
        Jabc = Int32[]
        Vabc = ComplexF64[]

        for iedge in eachindex(abcTags)
            edge = abcTags[iedge]
            orient = 1 + gmsh.model.mesh.getBasisFunctionsOrientationForElement(edge, "HcurlLegendre$basisOrder")

            jac, dets, globPoints = gmsh.model.mesh.getJacobian(edge, coord)
            jac = reinterpret(SMatrix{3,3,Float64,9}, jac)
            globPoints = reshape(globPoints, (3, numPoints))

            # nhat = J^-T * yhat (line element is xhat plus direction)
            nhats = [jac[i]' \ SA_F64[0, 1, 0] for i = 1:numPoints]
            nhats ./= norm.(nhats) # normalize (not sure if needed)

            globBasisFunctions = [cross(nhats[i], jac[i]' \ local_functions[f, i, orient]) for f = 1:numFunctions, i = 1:numPoints]
            globCurlBasisFunctions = [dot(nhats[i], curl_functions[f, i, orient] / dets[i]) for f = 1:numFunctions, i = 1:numPoints]

            Kabc = zeros(ComplexF64, (numFunctions, numFunctions))
            # tempM = first order part
            # tempS = second order abc part
            for j = 1:numFunctions, i = 1:numFunctions
                tempM::Float64 = 0.0
                # tempS::ComplexF64 = 0.0
                for q = 1:numPoints
                    tempM += weights[q] * globBasisFunctions[i, q]' * globBasisFunctions[j, q] * dets[q]

                    # R = norm(globPoints[:,q])
                    # # curvature of surface, not necessarily radial distance
                    # κ = 1/R
                    # tempS += 1 / (2*(1im*k + κ)) * weights[q] * globCurlBasisFunctions[i,q] * globCurlBasisFunctions[j,q] * dets[q]
                end
                Kabc[i, j] = 1im * k * tempM# + tempS
            end

            # assemble COO part
            globalTri = edgeToEleMap[edge]
            for j = 1:numFunctions
                globJ = globalFunctionIdxs[j, globalTri]
                for i = 1:numFunctions
                    globI = globalFunctionIdxs[i, globalTri]
                    push!(Iabc, globI)
                    push!(Jabc, globJ)
                    push!(Vabc, Kabc[i, j])
                end
            end
        end


        # append!(I, Iabc)
        # append!(J, Jabc)

        # append!(VK, Vabc)
        # append!(VM, -Vabc)

        Kabc = sparse(Iabc, Jabc, Vabc, size(Kglob, 1), size(Kglob, 2))
        #Kglob += Kabc

        abcDOF = unique(Iabc)
        # Mglob[abcDOF, :] .= 0
        # Mglob[:, abcDOF] .= 0
        #Mglob -= Kabc
    end


    return Kglob, Mglob, globalPecTags, globalFunctionIdxs, Kabc
end

function createGeometry()
    square_id = gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    gmsh.model.addPhysicalGroup(2, [square_id], 3000, "freespace")
    square_boundary = gmsh.model.getBoundary([(2, square_id)], true, false, false)
    gmsh.model.addPhysicalGroup(1, getfield.(square_boundary, 2), 2000, "pec")
    # gmsh.fltk.run()
end

function createCircle()
    circle_id = gmsh.model.occ.addCircle(0, 0, 0, 1)
    curve_loop_id = gmsh.model.occ.addCurveLoop([circle_id])
    surface_id = gmsh.model.occ.addPlaneSurface([curve_loop_id])

    p1 = gmsh.model.occ.addPoint(-0.25, -0.25, 0)
    p2 = gmsh.model.occ.addPoint(-0.25, 0.25, 0)
    p3 = gmsh.model.occ.addPoint(0.25, 0.25, 0)
    p4 = gmsh.model.occ.addPoint(0.25, -0.25, 0)

    l1 = gmsh.model.occ.addLine(p1, p2)
    l2 = gmsh.model.occ.addLine(p2, p3)
    l3 = gmsh.model.occ.addLine(p3, p4)

    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    gmsh.model.mesh.embed(1, [l1, l2, l3], 2, surface_id)

    gmsh.model.addPhysicalGroup(1, [l1, l2, l3], 2000, "pec")

    gmsh.model.addPhysicalGroup(2, [surface_id], 3000, "freespace")
    circle_boundary = gmsh.model.getBoundary((2, circle_id), true, false, false)
    circle_boundary = getfield.(circle_boundary, 2)
    gmsh.model.addPhysicalGroup(1, circle_boundary, 2001, "abc")

    # gmsh.fltk.run()
end

function createMesh(geoOrder::Integer, lc::Float64)
    # mesh size
    points = gmsh.model.getEntities(0)
    gmsh.model.mesh.setSize(points, lc)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(geoOrder)
    # gmsh.write("PECCircle.msh")
    # gmsh.fltk.run()
end

function setupEdges(triangleType::Integer, nTris::Integer)
    # create edge numbering
    gmsh.model.mesh.createEdges()

    # get edges for each element
    edgeNodes = gmsh.model.mesh.getElementEdgeNodes(triangleType, -1, true)
    edgeTags, _ = gmsh.model.mesh.getEdges(edgeNodes)
    edgeTags = convert.(Int32, edgeTags)

    edgesPerElement = length(edgeTags) ÷ nTris
    edgeTags = reshape(edgeTags, (edgesPerElement, nTris))

    return edgeTags
end


function main()
    geoOrder::Int = 2
    basisOrder::Int = 2
    basisName::String = "HcurlLegendre$basisOrder"
    curlBasisName::String = "CurlHcurlLegendre$basisOrder"

    gmsh.initialize()

    wavelength = c0 / (2 * π * ω)
    lc = wavelength / 10
    println("lambda by 10 ", lc)

    # createGeometry()
    createCircle()
    createMesh(geoOrder, lc)

    triTags, triNodes = getElements(2, 3000)
    pecTags, _ = getElements(1, 2000)
    abcTags, _ = getElements(1, 2001)
    # convert abc to pec for testing
    # append!(pecTags, abcTags)
    # abcTags = []

    nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
    nodeCoords = reshape(nodeCoords, (3, length(nodeTags)))

    nTris = length(triTags)

    localTriangle = createLocalTriangle(geoOrder, basisOrder)
    edgeTags = setupEdges(localTriangle.triangleType, nTris)

    # assemble per element matrices
    println("Number of triangles: $(length(triTags))")
    @time M, S = assembleLocalMatrices(triTags, localTriangle)

    @time K, M, globalPecTags, globalFunctionIdxs, Kabc = assembleGlobalMatrices(triTags, edgeTags, localTriangle, M, S, abcTags=abcTags, pecTags=pecTags)

    K += Kabc
    # M -= Kabc
    
    # integral, tag = addPointSource(0.0, 0.0, 0, basisOrder)
    # idx = searchsortedfirst(triTags, tag)

    # b = zeros(size(K,1))
    # b[globalFunctionIdxs[:,idx]] = integral
    # @time u = K \ b
    # plotHighOrder(3, localTriangle.triangleType, basisName, triTags, triNodes, nodeCoords, u, globalFunctionIdxs)

    """
    Computing eigenpairs of operator.
    We want the smallest pairs of the system Ax=λBx.
    the largest eigenvalues of the function y=A^-1Bx,
    correspond to the same eigenvectors of the original problem.
    """
    N = size(K, 1)
    # remove PEC dof before factorizing matrix
    nonPECDOF = setdiff(1:N, globalPecTags)
    Kpec = K[nonPECDOF, nonPECDOF]
    Mpec = M[nonPECDOF, nonPECDOF]
    Kabcpec = Kabc[nonPECDOF, nonPECDOF]

    # we have to use K^-1 * M * x as geneigsolve expects M to be hermitian positive definite (not possible with ABC or loss)
    # use factorize as it checks properties of matrix
    f = factorize(Kpec)
    matvec(x) = f \ (Mpec * x)

    # create big matvec
    bigN = 2 * size(Kpec, 1)
    k0::Float64 = sqrt(k2)
    
    C = (-2*Mpec + Kabcpec) ./ k0
    Mpec ./= k0^2
    
    function bigmatvec(x)
        N = size(x, 1)
        halfN = N ÷ 2

        vn = @view x[1:halfN]
        zn = @view x[1+halfN:end]

        # solve first block equation for zn'
        # M zn' = M vn
        result = zeros(ComplexF64, size(x))
        result[1+halfN:end] .= vn

        # solve second block equation for vn'
        # K vn' + C zn' = M zn
        # form rhs = M zn - C zn'
        # note x is [vn; zn]
        # copy rhs to first block of result
        result[1:halfN] = Mpec * zn .- C * vn
        # overwrite result using factorization of K
        ldiv!(f, @view result[1:halfN])
        return result
    end

    println("k0 ", k0)

    # key to efficiency of this algorithm is the EigSorter to filter out lambdas with real ≈ -1.0
    # these modes are degenerate and degrade performance significantly -> high multiplicity of eigenvalue=-1.0
    # λ, ϕpec, info = eigsolve(matvec, size(Kpec, 1), 25, EigSorter(f -> (real(1/f) <= -0.99)*abs(f)), ishermitian=false, krylovdim=60)
    λ, bigϕpec, info = eigsolve(bigmatvec, bigN, 25, EigSorter(f -> (real(1 / f) <= -0.99*k0*k0) * abs(f)), ishermitian=false, krylovdim=60)


    # put non pec degrees of freedom back into whole system degrees of freedom
    ϕ = zeros(ComplexF64, N, size(bigϕpec, 1))
    for i in eachindex(bigϕpec)
        ϕ[nonPECDOF, i] = bigϕpec[i][1:bigN÷2]
    end
    λ .= 1 ./ λ
    display(λ)

    # fig = Figure()
    # ax1 = Axis(fig[1,1], ylabel="Re(λ)", xlabel="n")
    # ax2 = Axis(fig[1,2], ylabel="Im(λ)", xlabel="n")
    # plot!(ax1, real.(λ))
    # plot!(ax2, imag.(λ))
    # display(fig)

    # for i=1:10
    #     plotHighOrder(2, localTriangle.triangleType, basisName, triTags, triNodes, nodeCoords, ϕ[:,i], globalFunctionIdxs)
    # end

    saveVTK(localTriangle, triTags, ϕ[:, :], globalFunctionIdxs)

    gmsh.finalize()
end

try
    gmsh.finalize()
catch

end

main()