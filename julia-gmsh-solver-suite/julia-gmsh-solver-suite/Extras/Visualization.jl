module Visualization

using ..FunctionSpaces: SimpleFiniteElementSpace, Triangle, Tetrahedron, HCurlTriangle, HCurlElement, HCurlTetrahedron
using Gmsh: gmsh
using StaticArrays
using SparseArrays
using WriteVTK
using HDF5
using Printf

const VecF = SVector{3,Float64}
const MatF = SMatrix{3,3,Float64,9}
const VecI{T} = SVector{T,Int64}

include("MakieViz.jl")

export plotHighOrder, MultiPlot


function interpolate(
    space::SimpleFiniteElementSpace{F},
    eval_points) where {F<:HCurlElement}
    x = Float64[]
    y = Float64[]
    z = Float64[]

    # get basic element information
    properties = space.mesh.properties
    elType = properties.type
    basis_name = "HcurlLegendre$(space.local_space.order)"

    # evaluate local basis functions
    num_points = length(eval_points) ÷ 3
    num_components, basis_functions, num_orientations = gmsh.model.mesh.getBasisFunctions(elType, eval_points, basis_name)
    @assert num_components == 3
    num_functions = length(basis_functions) ÷ num_components ÷ num_orientations ÷ num_points
    basis_functions = reshape(reinterpret(VecF, basis_functions), (num_functions, num_points, num_orientations))

    # global function idxs
    globalFunctionIdxs = space.dof_map

    rows = Int64[]
    cols = Int64[]
    vals = VecF[]

    currentRow = 1
    # map basis function to global triangles in every triangle at centerPoints
    for (t, tag) in enumerate(space.mesh.elements.tags)
        jac, _, coord = gmsh.model.mesh.getJacobian(tag, eval_points)
        jac = reinterpret(MatF, jac)
        coord = reshape(coord, (3, num_points))

        orient = gmsh.model.mesh.getBasisFunctionsOrientationForElement(tag, "HcurlLegendre")
        orient += 1

        globFunctions = [jac[p]' \ basis_functions[f, p, orient] for f = 1:num_functions, p = 1:num_points]
        for p = 1:num_points
            for f = 1:num_functions
                field = globFunctions[f, p]

                push!(rows, currentRow)
                push!(cols, globalFunctionIdxs[f, t])
                push!(vals, field)
            end

            push!(x, coord[1, p])
            push!(y, coord[2, p])
            push!(z, coord[3, p])
            currentRow += 1
        end
    end

    Interpolator = sparse(rows, cols, vals)

    if F == HCurlTriangle
        return Interpolator, x, y
    elseif F == HCurlTetrahedron
        return Interpolator, x, y, z
    end
end

function tesselateTriangle(levels::Integer)
    nodes = VecF[
        SA_F64[0, 0, 0],
        SA_F64[1, 0, 0],
        SA_F64[0, 1, 0]
    ]
    tris = VecI{3}[SA[1, 2, 3]]

    npop = 1
    for l = 1:levels
        nadded = 0
        for iold = 1:npop
            tri = popfirst!(tris)

            # create 3 new nodes
            n4 = (nodes[tri[1]] + nodes[tri[2]]) / 2
            n5 = (nodes[tri[2]] + nodes[tri[3]]) / 2
            n6 = (nodes[tri[1]] + nodes[tri[3]]) / 2

            push!(nodes, n4)
            i4 = length(nodes)

            push!(nodes, n5)
            i5 = length(nodes)

            push!(nodes, n6)
            i6 = length(nodes)

            # add four new triangles
            push!(tris, SA[tri[1], i4, i6])
            push!(tris, SA[i4, tri[2], i5])
            push!(tris, SA[i4, i5, i6])
            push!(tris, SA[i6, i5, tri[3]])
            nadded += 4
        end
        npop = nadded
    end

    return tris, nodes
end

function tesselateTetrahedron(levels::Integer)
    nodes = VecF[
        SA_F64[0, 0, 0],
        SA_F64[1, 0, 0],
        SA_F64[0, 1, 0],
        SA_F64[0, 0, 1]
    ]

    tets = VecI{4}[SA[1, 2, 3, 4]]

    npop = 1
    for l = 1:levels
        nadded = 0
        for iold = 1:npop
            tet = popfirst!(tets)

            i1, i2, i3, i4 = tet[1], tet[2], tet[3], tet[4]

            # new nodes
            n5 = (nodes[i1] + nodes[i2]) / 2
            n6 = (nodes[i2] + nodes[i3]) / 2
            n7 = (nodes[i1] + nodes[i3]) / 2
            n8 = (nodes[i1] + nodes[i4]) / 2
            n9 = (nodes[i2] + nodes[i4]) / 2
            n10 = (nodes[i3] + nodes[i4]) / 2

            push!(nodes, n5)
            i5 = length(nodes)
            push!(nodes, n6)
            i6 = length(nodes)
            push!(nodes, n7)
            i7 = length(nodes)
            push!(nodes, n8)
            i8 = length(nodes)
            push!(nodes, n9)
            i9 = length(nodes)
            push!(nodes, n10)
            i10 = length(nodes)

            # add new tetrahedrons
            # very bottom layer
            push!(tets, SA[i1, i5, i7, i8])
            push!(tets, SA[i5, i2, i6, i9])
            push!(tets, SA[i7, i6, i3, i10])
            # middle of tetra
            push!(tets, SA[i7, i8, i9, i10])
            push!(tets, SA[i6, i7, i9, i10])
            push!(tets, SA[i5, i7, i8, i9])
            push!(tets, SA[i5, i6, i7, i9])

            # top
            push!(tets, SA[i8, i9, i10, i4])


            nadded += 8
        end
        npop = nadded
    end

    return tets, nodes
end

function tesselateMesh(eleTags, localNodeIdxs, tessPoints)
    globTessNodes = Float64[] # will reshape at end
    globTessEles = Int64[]

    nodes_per_element = size(localNodeIdxs, 1)

    for t in eachindex(eleTags)
        # just evaluate the shape function
        _, _, coord = gmsh.model.mesh.getJacobian(eleTags[t], tessPoints)
        offset = length(globTessNodes) ÷ 3
        globNodeIdxs = localNodeIdxs .+ offset

        append!(globTessNodes, coord)
        append!(globTessEles, globNodeIdxs)
    end
    globTessNodes = reshape(globTessNodes, (3, length(globTessNodes) ÷ 3))
    globTessEles = reshape(globTessEles, (nodes_per_element, length(globTessEles) ÷ nodes_per_element))

    return globTessEles, globTessNodes
end

function saveVTK(
    file_name,
    array_name,
    space::SimpleFiniteElementSpace{HCurlTriangle},
    solution_coeffs::Vector{Vector{Tv}};
    levels=0,
    times::Vector=Float64[]) where {Tv<:Union{Float64,ComplexF64}}

    # grab stuff from FunctionSpace
    triTags = space.mesh.elements.tags

    ## tesselate the mesh
    tessTris, tessNodes = tesselateTriangle(levels)
    # create centroids (as VecF) in order to evaluate basis functions
    tessCentroids = [
        (tessNodes[tri[1]] + tessNodes[tri[2]] + tessNodes[tri[3]]) / 3
        for tri in tessTris
    ]

    tessTris = reshape(reinterpret(Int64, tessTris), (3, length(tessTris)))
    tessNodes = reinterpret(Float64, tessNodes)
    tessCentroids = reinterpret(Float64, tessCentroids)

    globalTris, globalNodes = tesselateMesh(triTags, tessTris, tessNodes)

    Interpolator, _, _ = interpolate(space, tessCentroids)

    # h5open("$file_name.hdf", "w") do file
    file = h5open("$file_name.hdf", "w")
    VTKHDF = create_group(file, "VTKHDF")
    HDF5.attributes(VTKHDF)["Version"] = [2, 0]

    type = "UnstructuredGrid"
    dspace = HDF5.dataspace(type)
    dtype = HDF5.datatype(type)
    HDF5.h5t_set_cset(dtype, HDF5.H5T_CSET_ASCII)
    attr = create_attribute(VTKHDF, "Type", dtype, dspace)
    write_attribute(attr, dtype, type)
    
    Nparts = max(1, length(times))

    VTKHDF["NumberOfConnectivityIds"] = [length(globalTris) for i=1:Nparts]
    VTKHDF["NumberOfPoints"] = [size(globalNodes, 2)  for i=1:Nparts]
    VTKHDF["NumberOfCells"] = [size(globalTris, 2)  for i=1:Nparts]
    VTKHDF["Points"] = globalNodes
    # triangle is type 5
    VTKHDF["Types"] = UInt8[5 for _ in axes(globalTris, 2)]
    VTKHDF["Connectivity"] = globalTris[:] .- 1
    VTKHDF["Offsets"] = [3 * i for i in 0:size(globalTris, 2)]

    CellData = create_group(VTKHDF, "CellData")

    if length(times) > 0
        # setup steps group
        NSteps = length(times)
        cell_data_offset = size(globalTris, 2)

        Steps = create_group(VTKHDF, "Steps")
        HDF5.attributes(Steps)["NSteps"] = NSteps
        Steps["Values"] = times

        # offsets that are just all zeros
        Steps["PartOffsets"] = zeros(Int64, NSteps)
        Steps["PointOffsets"] = zeros(Int64, NSteps)
        Steps["CellOffsets"] = zeros(Int64, 1, NSteps)
        Steps["ConnectivityIdOffsets"] = zeros(Int64, 1, NSteps)

        # offsets into cell arrays
        CellDataOffsets = create_group(Steps, "CellDataOffsets")
        # number of cells is also offset
        CellDataOffsets[array_name] = [cell_data_offset*(i-1) for i in 1:NSteps]

        TimeSeries = create_dataset(CellData, array_name, Tv, (3, NSteps * cell_data_offset), chunk=(3, cell_data_offset))

        for i in eachindex(solution_coeffs)
            interpolation = Interpolator * solution_coeffs[i]
            interpolation = reshape(reinterpret(Tv, interpolation), (3, length(interpolation)))

            chunk_start = 1 + cell_data_offset * (i - 1)
            chunk_end = cell_data_offset * i
            TimeSeries[:, chunk_start:chunk_end] = interpolation
        end
    else
        for i in eachindex(solution_coeffs)
            interpolation = Interpolator * solution_coeffs[i]
            interpolation = reshape(reinterpret(Tv, interpolation), (3, length(interpolation)))

            name = @sprintf "%s_%03d" array_name i
            CellData[name] = real.(interpolation)

        end
    end

    display(file)
    # end
    close(file)
end

function saveVTK(
    file_name,
    array_name,
    space::SimpleFiniteElementSpace{HCurlTetrahedron},
    solution_coeffs::Vector{Vector{Tv}};
    levels=0,
    times::Vector=Float64[]) where {Tv<:Union{Float64,ComplexF64}}

    tetTags = space.mesh.elements.tags

    tessTets, tessNodes = tesselateTetrahedron(levels)

    tessCentroids = [
        (tessNodes[tet[1]] + tessNodes[tet[2]] + tessNodes[tet[3]] + tessNodes[tet[4]]) / 4
        for tet in tessTets
    ]

    tets = reshape(reinterpret(Int64, tessTets), (4, length(tessTets)))
    tessNodes = reinterpret(Float64, tessNodes)
    tessCentroids = reinterpret(Float64, tessCentroids)

    globalTess, globalNodes = tesselateMesh(tetTags, tets, tessNodes)

    Interpolator, _, _, _ = interpolate(space, tessCentroids)

    h5open("$file_name.hdf", "w") do file
        VTKHDF = create_group(file, "VTKHDF")
        HDF5.attributes(VTKHDF)["Version"] = [2, 0]

        type = "UnstructuredGrid"
        dspace = HDF5.dataspace(type)
        dtype = HDF5.datatype(type)
        HDF5.h5t_set_cset(dtype, HDF5.H5T_CSET_ASCII)
        attr = create_attribute(VTKHDF, "Type", dtype, dspace)
        write_attribute(attr, dtype, type)

        VTKHDF["NumberOfConnectivityIds"] = [length(globalTess)]
        VTKHDF["NumberOfPoints"] = [size(globalNodes, 2)]
        VTKHDF["NumberOfCells"] = [size(globalTess, 2)]
        VTKHDF["Points"] = globalNodes
        # tetrahedron is type 10
        VTKHDF["Types"] = UInt8[10 for _ in axes(globalTess, 2)]
        VTKHDF["Connectivity"] = globalTess[:] .- 1
        VTKHDF["Offsets"] = [4 * i for i in 0:size(globalTess, 2)]

        CellData = create_group(VTKHDF, "CellData")

        for i in eachindex(solution_coeffs)
            interpolation = Interpolator * solution_coeffs[i]
            interpolation = reshape(reinterpret(Tv, interpolation), (3, length(interpolation)))

            name = @sprintf "%s_%03d" array_name i
            CellData[name] = real.(interpolation)
        end
        display(file)
    end
end

function savePVD(
    name,
    space::SimpleFiniteElementSpace{HCurlTriangle},
    solution_coeffs::Vector{Vector{Tv}},
    times::Vector{Float64};
    levels=0) where {Tv<:Union{Float64,ComplexF64}}

    # grab stuff from FunctionSpace
    triTags = space.mesh.elements.tags

    ## tesselate the mesh
    tessTris, tessNodes = tesselateTriangle(levels)
    # create centroids (as VecF) in order to evaluate basis functions
    tessCentroids = [
        (tessNodes[tri[1]] + tessNodes[tri[2]] + tessNodes[tri[3]]) / 3
        for tri in tessTris
    ]

    tessTris = reshape(reinterpret(Int64, tessTris), (3, length(tessTris)))
    tessNodes = reinterpret(Float64, tessNodes)
    tessCentroids = reinterpret(Float64, tessCentroids)

    globalTris, globalNodes = tesselateMesh(triTags, tessTris, tessNodes)

    Interpolator, _, _ = interpolate(space, tessCentroids)

    cells = [
        MeshCell(VTKCellTypes.VTK_TRIANGLE, globalTris[1:3, i])
        for i in axes(globalTris, 2)
    ]

    paraview_collection(name) do pvd
        for (it, t) in enumerate(times)
            field = Interpolator * solution_coeffs[it]
            vtk_grid("timestep_$it", globalNodes, cells) do vtk

                if Tv == ComplexF64
                    field = real.(field)
                end
                vtk["E"] = field
                pvd[t] = vtk
            end
        end
    end

end

export saveVTK, savePVD


end