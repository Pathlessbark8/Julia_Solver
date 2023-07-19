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
using ForwardDiff
using FFTW

include("gmsh_mesh_utils.jl")
include("gmsh_basis_utils.jl")
include("visualization.jl")
include("analytic.jl")

ϵ0::Float64 = 8.854e-12
μ0::Float64 = 4 * π * 1e-7
c0::Float64 = 3.0e8
σ0::Float64 = 1e-4


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
function assembleGlobalMatricesTimeDomain(
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
    VT = Float64[]
    VR = Float64[]
    VS = Float64[]

    for t in eachindex(triTags)
        for jfunc = 1:nEdgeFuncs+nFaceFuncs
            for ifunc = 1:nEdgeFuncs+nFaceFuncs
                vt = ϵ0 * M[ifunc, jfunc, t]
                vr = σ0 * M[ifunc, jfunc, t]
                vs = 1/μ0 * S[ifunc, jfunc, t]
                push!(I, globalFunctionIdxs[ifunc, t])
                push!(J, globalFunctionIdxs[jfunc, t])

                push!(VT, vt)
                push!(VR, vr)
                push!(VS, vs)
            end
        end
    end

    T = sparse(I, J, VT)
    R = sparse(I, J, VR)
    S = sparse(I, J, VS)

    # contribution from ABC bc's
    # ABC contribute to R
    # PEC dof are removed from whole system after
    # admittance
    Y = sqrt(ϵ0 / μ0)
    if length(abcTags) > 0
        sortedABCs = sort(abcTags)
        edgeToEleMap = Dict{Int32,Int32}()
        for itri in eachindex(triTags)
            for t in edgeTags[:, itri]
                e = searchsortedfirst(sortedABCs, t)
                if e <= length(sortedABCs)
                    edgeToEleMap[t] = itri
                end
            end
        end 

        # setup local basis functions on edge
        basisOrder = localTriangle.basis_order
        gaussOrder = 2 * (basisOrder + 1)

        lineType = gmsh.model.mesh.getElementType("line", localTriangle.geo_order)
        coord, weights = gmsh.model.mesh.getIntegrationPoints(lineType, "Gauss$gaussOrder")
        numPoints = length(weights)

        # get functions on edge
        numComponents, local_functions, numOrientations = gmsh.model.mesh.getBasisFunctions(lineType, coord, "HcurlLegendre$basisOrder")
        numFunctions = length(local_functions) ÷ numComponents ÷ numOrientations ÷ numPoints
        local_functions = reshape(reinterpret(SVector{3,Float64}, local_functions), (numFunctions, numPoints, numOrientations))
        

        # for each abc edge
        Iabc = Int32[]
        Jabc = Int32[]
        Vabc = Float64[]
        for iedge in eachindex(sortedABCs)
            edgeTag = sortedABCs[iedge]
            edgeOrientation = 1 + gmsh.model.mesh.getBasisFunctionsOrientationForElement(edgeTag, "HcurlLegendre$basisOrder")
            
            # get jacobian along edge
            jac, dets, globPoints = gmsh.model.mesh.getJacobian(edgeTag, coord)
            jac = reinterpret(SMatrix{3,3,Float64,9}, jac)
            globPoints = reshape(globPoints, (3, numPoints))

            # nhat = J^-T * yhat (line element is xhat plus direction)
            nhats = [jac[i]' \ SA_F64[0, 1, 0] for i = 1:numPoints]
            nhats ./= norm.(nhats) # normalize (not sure if needed)

            globBasisFunctions = [cross(nhats[i], jac[i]' \ local_functions[f, i, edgeOrientation]) for f = 1:numFunctions, i = 1:numPoints]
            
            Rabc = zeros(numFunctions, numFunctions)
            for j = 1:numFunctions, i = 1:numFunctions
                tempM::Float64 = 0.0
                for q = 1:numPoints
                    tempM += weights[q] * (globBasisFunctions[i, q]' * globBasisFunctions[j, q]) * dets[q]
                end
                Rabc[i, j] = Y * tempM
            end

            # assemble COO part
            globalTri = edgeToEleMap[iedge]
            for j = 1:numFunctions
                globJ = globalFunctionIdxs[j, globalTri]
                for i = 1:numFunctions
                    globI = globalFunctionIdxs[i, globalTri]
                    push!(Iabc, globI)
                    push!(Jabc, globJ)
                    push!(Vabc, Rabc[i, j])
                end
            end
        end

        # assemble global abc part
        Rabc = sparse(Iabc, Jabc, Vabc, size(R,1), size(R,2))
        R += Rabc
    end

    return T, R, S, globalFunctionIdxs, globalPecTags
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

function gaussDeriv(t, σ)
    return 1/σ * t .* exp.(-(t./σ).^2 ./ 2)
end

function gaussSin(t, σ, ωc)
    return exp(-(t/σ)^2 / 2) * sin(ωc * t)
end

function main()
    geoOrder::Int = 2
    basisOrder::Int = 2
    basisName::String = "HcurlLegendre$basisOrder"
    curlBasisName::String = "CurlHcurlLegendre$basisOrder"

    gmsh.initialize()

    f = 100e6
    wavelength = c0 / (2 * π * f)
    lc = wavelength / 10
    println("lambda by 10 ", lc)

    
    # createGeometry()
    createCircle()
    createMesh(geoOrder, lc)
    # gmsh.fltk.run()
    # return

    triTags, triNodes = getElements(2, 3000)
    pecTags, _ = getElements(1, 2000)
    abcTags, _ = getElements(1, 2001)
    # convert abc to pec for testing
    # append!(pecTags, abcTags)
    # abcTags = []
    # pecTags = []

    nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
    nodeCoords = reshape(nodeCoords, (3, length(nodeTags)))

    nTris = length(triTags)

    localTriangle = createLocalTriangle(geoOrder, basisOrder)
    edgeTags = setupEdges(localTriangle.triangleType, nTris)

    # assemble per element matrices
    println("Number of triangles: $(length(triTags))")
    @time M, S = assembleLocalMatrices(triTags, localTriangle)

    @time T, R, S, globalDOFMap, globalPecTags = assembleGlobalMatricesTimeDomain(triTags, edgeTags, localTriangle, M, S, abcTags=abcTags, pecTags=pecTags)

    # Β ≥ 0.25 → unconditionally stable
    Beta = 0.25
    fs = 1e10
    Δt = 0.5 / fs
    Δx = 1e-2
    Δt = Δx / c0
    println("Δt=$Δt 1/Δt=$(1/Δt)")

    A = T ./ (Δt^2) + R ./ (2 * Δt) + Beta .* S
    B = T ./ (2 / Δt^2) - (1 - 2 * Beta) .* S
    C = -(T ./ (Δt^2) - R ./ (2 * Δt) + Beta .* S)

    # remove PEC DOF
    N = size(A, 1)
    A[globalPecTags, :] .= 0
    A[:, globalPecTags] .= 0
    A[[CartesianIndex(i,i) for i in globalPecTags]] .= 1

    B[globalPecTags, :] .= 0
    B[:, globalPecTags] .= 0

    C[globalPecTags, :] .= 0
    C[:, globalPecTags] .= 0
    # A * E^{n+1} = B * E^{n} + C * E^{n-1} + Β * f^{n+1} + (1-2*Β) * f^{n} + Β * f^{n-1}

    # source terms
    
    # J(t) = gaussDeriv(t, σ)
    ωc = 2*pi*300e6 
    σ = 1 / (2*pi * 50e6)
    J(t) = gaussSin(t, σ, ωc)
    Jt(t) = ForwardDiff.derivative(J, t)
    # Jt(t) = J(t)

    tstart = -4*σ
    t = tstart:Δt:-tstart
    
    display(length(t))

    signal = J.(t)
    F = fftshift(fft(signal))
    freqs = fftshift(fftfreq(length(t), fs))
    p = lines(t, signal, label="f(t)")
    display(p)
    p = lines(freqs, abs.(F), label="F(ω)")
    display(p)
    # return

    f = factorize(A)
    display(f)

    Ecurr = zeros(N)
    Eprev = zeros(N)
    Enext = zeros(N)
    rhs = zeros(N)

    # create source in center
    localR, tag = addPointSource(-0.1, -0.1, 0, localTriangle.basis_order, curl=false)
    idx = searchsortedfirst(triTags, tag)
    # return
    sourceIdx = globalDOFMap[:, idx]
    SourceIntegral = sparsevec(sourceIdx, localR, N)
    
    sol = [Float64[]]
    solTimes = []
    
    for (i,tnext) in enumerate(t)
        scalarSource = Beta * Jt(tnext) + (1-2*Beta)*Jt(tnext - Δt) + Beta*Jt(tnext - 2*Δt)
        rhs .= B * Ecurr + C * Eprev - SourceIntegral*scalarSource

        Enext = f \ rhs
        
        if i % 5 == 1
            push!(sol, Ecurr)
            push!(solTimes, tnext - Δt)
        end

        println("i=$(i) Enorm = $(norm(Ecurr))")

        # update old
        Eprev = Ecurr
        Ecurr = Enext
    end

    popfirst!(sol)
    
    println("saving solution to paraview")
    savePVD(localTriangle, triTags, globalDOFMap, sol, solTimes)
    gmsh.finalize()
end

try
    gmsh.finalize()
catch

end

main()