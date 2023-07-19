struct Triple{T}
    a::T
    b::T
    c::T
end
Base.:+(x::Triple, y::Triple) = Triple(x.a+y.a, x.b+y.b, x.c+y.c)
Base.:/(x::Triple, y::Number) = Triple(x.a/y, x.b/y, x.c/y)

function tesselateTriangle(levels)
    nodes = Triple{Float64}[
        Triple(0.0,0.0,0.0),
        Triple(1.0,0.0,0.0),
        Triple(0.0,1.0,0.0)
    ]
    tris = Triple{Int64}[Triple(1,2,3)]

    npop = 1
    for l=1:levels
        nadded = 0
        for iold = 1:npop
            tri = popfirst!(tris)

            # create 3 new nodes
            n4 = (nodes[tri.a] + nodes[tri.b]) / 2
            n5 = (nodes[tri.b] + nodes[tri.c]) / 2
            n6 = (nodes[tri.a] + nodes[tri.c]) / 2
            
            push!(nodes, n4)
            i4 = length(nodes)

            push!(nodes, n5)
            i5 = length(nodes)

            push!(nodes, n6)
            i6 = length(nodes)

            # add four new triangles
            push!(tris, Triple(tri.a, i4, i6))
            push!(tris, Triple(i4, tri.b, i5))
            push!(tris, Triple(i4, i5, i6))
            push!(tris, Triple(i6, i5, tri.c))
            nadded += 4
        end
        npop = nadded
    end
    return tris, nodes
end

function tesselateMesh(triTags, localNodeIdxs, tessPoints)
    globTessNodes = Float64[] # will reshape at end
    globTessTris = Int64[]
    
    for t in eachindex(triTags)
        # just evaluate the shape function
        _, _, coord = gmsh.model.mesh.getJacobian(triTags[t], tessPoints)
        offset = length(globTessNodes) ÷ 3
        globNodeIdxs = localNodeIdxs .+ offset

        append!(globTessNodes, coord)
        append!(globTessTris, globNodeIdxs)
    end
    globTessNodes = reshape(globTessNodes, (3, length(globTessNodes) ÷ 3))
    globTessTris = reshape(globTessTris, (3, length(globTessTris) ÷ 3))
    
    return globTessTris, globTessNodes
end

function saveVTK(
    localTriangle::HCurlTriangle,
    triTags::Vector,
    ϕ::VecOrMat,
    globalFunctionIdxs::Array)
    basisName = "HcurlLegendre$(localTriangle.basis_order)"
    
    localTessTris, localTessNodes = tesselateTriangle(2)
    centroids = [
        (localTessNodes[t.a] + localTessNodes[t.b] + localTessNodes[t.c])/3 
        for t in localTessTris
    ]
    centroids = reinterpret(Float64, centroids)
    tessPoints = reinterpret(Float64, localTessNodes)
    localNodeIdxs = reinterpret(Int64, localTessTris)

    nBasisPoints = length(localTessTris)
    
    globTessTris, globTessNodes = tesselateMesh(triTags, localNodeIdxs, tessPoints)
    
    # get local basis functions
    _, bfuncs, numOrient = gmsh.model.mesh.getBasisFunctions(
        localTriangle.triangleType, centroids, basisName)
    nFunctions = length(bfuncs) ÷ numOrient ÷ 3 ÷ nBasisPoints
    bfuncs = reshape(reinterpret(SVector{3,Float64}, bfuncs), (nFunctions, nBasisPoints, numOrient))
    
    # map basis function to global triangles in every triangle at center
    nFields = max(size(ϕ, 2), 1)
    fields = [ComplexF64[] for i=1:nFields]
    for t in eachindex(triTags)
        jac, _, _ = gmsh.model.mesh.getJacobian(triTags[t], tessPoints)
        jac = reinterpret(SMatrix{3,3,Float64,9}, jac)
        
        orient = gmsh.model.mesh.getBasisFunctionsOrientationForElement(triTags[t], basisName)
        orient += 1

        globFunctions = [jac[p]' \ bfuncs[f, p, orient] for f=1:nFunctions, p=1:nBasisPoints]
        for p=1:nBasisPoints
            # evalute the field at this point
            for i=1:nFields
                field = zeros(ComplexF64, 3)
            
                for f=1:nFunctions
                    field += globFunctions[f,p] .* ϕ[globalFunctionIdxs[f, t],i]
                end

                append!(fields[i], field)
            end
        end
    end
    
    vectorFields = [reshape(fields[i], (3, length(fields[i]) ÷ 3)) for i=1:nFields]
    
    cells = [
        MeshCell(VTKCellTypes.VTK_TRIANGLE, globTessTris[1:3, i]) 
        for i in axes(globTessTris, 2)
    ]

    vtk_grid("test.vtu", globTessNodes, cells) do vtk
        for count = 1:nFields
            
            vtk[@sprintf("E%02d_real", count)] = real.(vectorFields[count])
            vtk[@sprintf("E%02d_imag", count)] = imag.(vectorFields[count])
        end
    end
end

# time series
function savePVD(localTriangle::HCurlTriangle,
    triTags::Vector,
    globalFunctionIdxs::Array,
    solutions::Vector{Vector{Float64}},
    times)

    basisName = "HcurlLegendre$(localTriangle.basis_order)"
    
    # create local tesselation points
    localTessTris, localTessNodes = tesselateTriangle(0)
    centroids = [
        (localTessNodes[t.a] + localTessNodes[t.b] + localTessNodes[t.c])/3 
        for t in localTessTris
    ]
    centroids = reinterpret(Float64, centroids)
    tessPoints = reinterpret(Float64, localTessNodes)
    localNodeIdxs = reinterpret(Int64, localTessTris)

    nBasisPoints = length(localTessTris)
    
    # create global tesselation points
    globTessTris, globTessNodes = tesselateMesh(triTags, localNodeIdxs, tessPoints)
    cells = [
        MeshCell(VTKCellTypes.VTK_TRIANGLE, globTessTris[1:3, i]) 
        for i in axes(globTessTris, 2)
    ]

    # get local basis functions
    _, bfuncs, numOrient = gmsh.model.mesh.getBasisFunctions(
        localTriangle.triangleType, centroids, basisName)
    nFunctions = length(bfuncs) ÷ numOrient ÷ 3 ÷ nBasisPoints
    bfuncs = reshape(reinterpret(SVector{3,Float64}, bfuncs), (nFunctions, nBasisPoints, numOrient))
    
    # map basis function to global triangles in every triangle at center
    # for time series, assemble a global matrix
    # T x P X F?
    println("Creating interpolations")
    evaluations = zeros(SVector{3, Float64}, nFunctions, nBasisPoints, length(triTags))
    for t in eachindex(triTags)
        jac, _, _ = gmsh.model.mesh.getJacobian(triTags[t], tessPoints)
        jac = reinterpret(SMatrix{3,3,Float64,9}, jac)
        
        orient = gmsh.model.mesh.getBasisFunctionsOrientationForElement(triTags[t], basisName)
        orient += 1

        globFunctions = [jac[p]' \ bfuncs[f, p, orient] for f=1:nFunctions, p=1:nBasisPoints]
        for p=1:nBasisPoints
            # evalute a field at this point
            for f=1:nFunctions
                evaluations[f, p, t] = globFunctions[f,p]
            end
        end
    end

    println("evaluating each timestep")
    
    mkpath("results")
    paraview_collection("results/full_simulation") do pvd
        for (n, time) in enumerate(times)
            fields = zeros(SVector{3,Float64}, nBasisPoints, length(triTags))
            println("solution n: $(n)")
            for t in eachindex(triTags)
                localSolution = solutions[n][globalFunctionIdxs[:, t]]
                for p=1:nBasisPoints
                    fields[p, t] = evaluations[:, p, t]' * localSolution
                end
            end
            fields = reshape(fields, nBasisPoints*length(triTags))
            vtk_grid("results/timestep_$n", globTessNodes, cells) do vtk
                vtk["E"] = fields
                pvd[time] = vtk
            end
        end
    end
end


function plotHighOrder(
    nSample::Integer,
    triangleType::Integer,
    basisName::String,
    triTags::Vector,
    triNodes::Array,
    nodeCoords::Array,
    ϕ::Vector,
    globalFunctionIdxs::Array
    )
    # create local points to sample local triangle
    step = 1.0 / nSample
    x = step/3:step:1.0
    y = step/3:step:1.0
    points = [(xi, yj, 0.0) for xi in x, yj in y if xi+yj<=1.0]
    nPoints = length(points)
    centerPoints = reinterpret(Float64, points)

    x = Float64[]
    y = Float64[]
    u = ComplexF64[]
    v = ComplexF64[]
    
    # get local basis functions
    _, bfuncs, numOrient = gmsh.model.mesh.getBasisFunctions(triangleType, centerPoints, basisName)
    nFunctions = length(bfuncs) ÷ numOrient ÷ 3 ÷ nPoints
    bfuncs = reshape(reinterpret(SVector{3,Float64}, bfuncs), (nFunctions, nPoints, numOrient))
    
    # map basis function to global triangles in every triangle at centerPoints
    for t in eachindex(triTags)
        jac, _, coord = gmsh.model.mesh.getJacobian(triTags[t], centerPoints)
        jac = reinterpret(SMatrix{3,3,Float64,9}, jac)
        coord = reshape(coord, (3, nPoints))

        orient = gmsh.model.mesh.getBasisFunctionsOrientationForElement(triTags[t], basisName)
        orient += 1

        globFunctions = [jac[p]' \ bfuncs[f, p, orient] for f =1:nFunctions, p=1:nPoints]
        for p=1:nPoints
            field = @SVector ComplexF64[0.0, 0.0, 0.0]
            for i=1:nFunctions
                field += globFunctions[i,p] * ϕ[globalFunctionIdxs[i, t]]
            end
            push!(u, field[1])
            push!(v, field[2])
            push!(x, coord[1,p])
            push!(y, coord[2,p])
        end
    end
    # add triangles
    # first three nodes are primary nodes
    faces = reshape(triNodes[1:3,:], (3, length(triTags)))'
    verts = nodeCoords[1:2,:]'

    len = hypot.(u,v)
    avg = median(len)
    u[len .> 10*avg] .= 0
    v[len .> 10*avg] .= 0
    len = hypot.(u,v)
    
    u .= real.(u) ./ len
    v .= real.(v) ./ len
    
    f = Figure()
    ax = Axis(f[1,1], aspect=1)
    poly!(ax,verts, faces, strokewidth=1, transparency=true, shading=true, color=RGBAf(0.0,0.0,0.0,0.0))
    arrows!(ax,x, y, u, v, arrowcolor=len, linecolor=len, arrowsize=10, lengthscale=1e-1)

    display(f)
end
