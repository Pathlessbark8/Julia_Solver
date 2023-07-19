using StaticArrays
using SparseArrays
using LinearAlgebra

"""
There is currently no integrator for Dirichlet or Neumann boundary conditions.
Zero Neumann is the default boundary condition for the H curl basis functions.
Dirichlet BC's are likely to be encountered for PEC, which we just set to zero, or for
scattered fields, for which we don't need to integrate as we'd already have the coefficients.
"""


"""
Assemble a GlobalFunctionSpace for a HCurl triangle and a SimpleMesh mesh.

This function assigns local functions in each element a global degree of freedom to 
be used for the finite element method.
"""
function HCurlSpace(mesh::SimpleMesh, local_space::HCurlTriangle)
    nElements = numElements(mesh.elements)

    # we need to create gmsh edges now, which requres knowing geometric order of mesh
    element_type = mesh.properties.type
    edgeTags = createEdges(element_type, nElements)

    # create global dof mapping
    maxEdgeId = maximum(edgeTags)
    nEdgeFunctions = numEdgeFunctions(local_space)
    nFaceFunctions = numFaceFunctions(local_space)

    # how each local function is mapped to local geometric dof
    localDOFMap = DOFMap(local_space)
    nFunctions = length(localDOFMap)
    @assert nFunctions == nEdgeFunctions + nFaceFunctions

    localEdgeIdxs = [t for (d, t) in localDOFMap if d == 1]
    funcsPerEdge = nEdgeFunctions ÷ numEdges(local_space)
    funcsPerFace = nFaceFunctions ÷ numFaces(local_space)

    # actually fill the table
    edgeOffset = funcsPerEdge * maxEdgeId
    global_dof_map = Array{Int64}(undef, (nFunctions, nElements))
    for t in eachindex(mesh.elements.tags)
        # for each edge function
        for ifunc = 1:nEdgeFunctions
            funcOrder = 1 + (ifunc - 1) % funcsPerEdge
            zeroOrderEdgeIndex = edgeTags[localEdgeIdxs[ifunc], t]
            highOrderEdgeIndex = funcsPerEdge * (zeroOrderEdgeIndex - 1) + funcOrder
            global_dof_map[ifunc, t] = highOrderEdgeIndex
        end
        # for each face function
        for ifunc = 1:nFaceFunctions
            funcOrder = 1 + (ifunc - 1) % funcsPerFace
            highOrderFaceIndex = edgeOffset + funcsPerFace * (t - 1) + funcOrder
            global_dof_map[nEdgeFunctions+ifunc, t] = highOrderFaceIndex
        end
    end

    # create global dof list for physical groups of dimension lower than 2
    physical_group_mapping = Dict{Tuple{Int32,Int32},Array{Int32,2}}()

    # could be abstracted, but this is currently only for edge dof
    for ((dim, physicsTag), elementTags) in mesh.physics_map
        # only for boundary conditions
        if dim == 1
            list = Int32[]
            for tag in elementTags
                zeroOrderEdgeIndex = tag
                for ifunc = 1:funcsPerEdge
                    funcOrder = 1 + (ifunc - 1) % funcsPerEdge
                    highOrderEdgeIndex = funcsPerEdge * (zeroOrderEdgeIndex - 1) + funcOrder

                    # use local side as an offset
                    # localSide = edgeSideMaping[tag]
                    # highOrderEdgeIndex += (localSide - 1) * funcsPerEdge

                    push!(list, highOrderEdgeIndex)
                end
            end
            physical_group_mapping[(dim, physicsTag)] = reshape(list, (funcsPerEdge, length(elementTags)))
        end
    end
    SimpleFiniteElementSpace{HCurlTriangle}(mesh, local_space, global_dof_map, physical_group_mapping)
end

"""
Assemble a GlobalFunctionSpace for a HCurlTetrahedron and a SimpleMesh mesh.

This function assigns local functions in each element a global degree of freedom to 
be used for the finite element method.
"""
function HCurlSpace(mesh::SimpleMesh, local_space::HCurlTetrahedron)
    nElements = numElements(mesh.elements)

    # we need to create gmsh edges now, which requres knowing geometric order of mesh
    element_type = mesh.properties.type
    edgeTags = createEdges(element_type, nElements)

    # we also need to create faces
    faceTags = createFaces(element_type, nElements)

    # create global dof mapping
    max_edge_id = maximum(edgeTags)
    max_face_id = maximum(faceTags)
    num_edge_funcs = numEdgeFunctions(local_space)
    num_face_funcs = numFaceFunctions(local_space)
    num_bubble_funcs = numBubbleFunctions(local_space)

    # how each local function is mapped to local geometric dof
    local_dof_map = DOFMap(local_space)
    num_functions = length(local_dof_map)
    @assert num_functions == num_edge_funcs + num_face_funcs + num_bubble_funcs

    local_edge_ids = [t for (d, t) in local_dof_map if d == 1]
    local_face_ids = [t for (d, t) in local_dof_map if d == 2]

    funcs_per_edge = num_edge_funcs ÷ numEdges(local_space)
    funcs_per_face = num_face_funcs ÷ numFaces(local_space)
    funcs_per_vol = num_bubble_funcs ÷ numVolumes(local_space)

    # actually fill the table
    edge_offset = funcs_per_edge * max_edge_id
    face_offset = funcs_per_face * max_face_id + edge_offset
    global_dof_map = Array{Int64}(undef, (num_functions, nElements))
    for t in eachindex(mesh.elements.tags)
        # for each edge function
        for ifunc = 1:num_edge_funcs
            funcOrder = 1 + (ifunc - 1) % funcs_per_edge
            mesh_edge_tag = edgeTags[local_edge_ids[ifunc], t]
            highOrderEdgeIndex = funcs_per_edge * (mesh_edge_tag - 1) + funcOrder
            global_dof_map[ifunc, t] = highOrderEdgeIndex
        end
        # for each face function
        for ifunc = 1:num_face_funcs
            funcOrder = 1 + (ifunc - 1) % funcs_per_face
            mesh_face_tag = faceTags[local_face_ids[ifunc], t]
            highOrderFaceIndex = funcs_per_face * (mesh_face_tag - 1) + funcOrder
            global_dof_map[num_edge_funcs+ifunc, t] = edge_offset + highOrderFaceIndex
        end
        # for each bubble function
        for ifunc = 1:num_bubble_funcs
            funcOrder = 1 + (ifunc - 1) % funcs_per_vol
            highOrderVolIndex = funcs_per_vol * (t - 1) + funcOrder
            global_dof_map[num_edge_funcs+num_face_funcs+ifunc, t] = face_offset + highOrderVolIndex
        end
    end

    # create global dof list for physical groups of dimension lower than 3
    physical_group_mapping = Dict{Tuple{Int32,Int32},Array{Int32,2}}()

    # loop over face elements
    mesh_order = mesh.properties.order
    for ((dim, physicsTag), elementTags) in mesh.physics_map
        # only for boundary conditions
        if dim == 2
            list = Int32[]
            # get the edges associated with triangular faces
            ents = gmsh.model.getEntitiesForPhysicalGroup(dim, physicsTag)
            face_type = gmsh.model.mesh.getElementType("Triangle", mesh_order)
            edgeNodeTags = UInt64[]
            for ent_tag in ents
                _tags = gmsh.model.mesh.getElementEdgeNodes(face_type, ent_tag, true)
                append!(edgeNodeTags, _tags)
            end
            face_edge_tags, _ = gmsh.model.mesh.getEdges(edgeNodeTags)
            face_edge_tags = reshape(face_edge_tags, (3, length(face_edge_tags) ÷ 3))

            for (iface, mesh_face_tag) in enumerate(elementTags)
                #for each edge function
                for ifunc = 1:3*funcs_per_edge # 3 edges per face
                    funcOrder = 1 + (ifunc - 1) % funcs_per_edge
                    iedge = 1 + (ifunc - 1) ÷ funcs_per_edge

                    # the actual edge tag of the gmsh edge 
                    mesh_edge_tag = face_edge_tags[iedge, iface]
                    # the high order edge tag
                    highOrderEdgeIndex = funcs_per_edge * (mesh_edge_tag - 1) + funcOrder

                    push!(list, highOrderEdgeIndex)
                end
                # for each face function
                for ifunc = 1:funcs_per_face
                    funcOrder = 1 + (ifunc - 1) % funcs_per_face
                    highOrderFaceIndex = edge_offset + funcs_per_face * (mesh_face_tag - 1) + funcOrder
                    push!(list, highOrderFaceIndex)
                end
            end
            physical_group_mapping[(dim, physicsTag)] = reshape(list, (3 * funcs_per_edge + funcs_per_face, length(elementTags)))
        end
    end
    SimpleFiniteElementSpace{HCurlTetrahedron}(mesh, local_space, global_dof_map, physical_group_mapping)
end

"""
Assemble local matrices.
"""
function integrate(space::SimpleFiniteElementSpace{HCurlElement{T}}) where {T<:AbstractElement}
    local_functions = space.local_space.functions
    curl_functions = space.local_space.curl_functions
    (_, num_functions, num_points, num_orientations) = length(local_functions)

    # convert to arrays of SVector
    local_functions = reinterpret(SVector{3,Float64}, local_functions)
    local_functions = reshape(local_functions, (num_functions, num_points, num_orientations))

    curl_functions = reinterpret(SVector{3,Float64}, curl_functions)
    curl_functions = reshape(curl_functions, (num_functions, num_points, num_orientations))

    # get quad points and weights
    quad_weights = space.local_space.quad_weights
    quad_points = space.local_space.quad_points

    # storage for M and S
    M = zeros(num_functions, num_functions, length(space.mesh.elements.tags))
    S = zeros(num_functions, num_functions, length(space.mesh.elements.tags))

    # temporary arrays
    basis_functions = zeros(SVector{3,Float64}, num_points, num_functions)
    curl_basis_functions = zeros(SVector{3,Float64}, num_points, num_functions)

    for (t, tag) in enumerate(space.mesh.elements.tags)
        # get orientation for this element
        orient = 1 + gmsh.model.mesh.getBasisFunctionsOrientationForElement(tag, "HcurlLegendre")

        # get jacobians for this element, convert to SMatrix
        jacobians, determinants, _ = gmsh.model.mesh.getJacobian(tag, quad_points)
        jacobians = reinterpret(SMatrix{3,3,Float64,9}, jacobians)

        for i = 1:num_points, f = 1:num_functions
            basis_functions[i, f] = jacobians[i]' \ local_functions[f, i, orient]
            curl_basis_functions[i, f] = (jacobians[i] * curl_functions[f, i, orient]) / determinants[i]
        end

        for j = 1:num_functions, i = 1:num_functions
            tempM = 0
            tempS = 0
            for q = 1:num_points
                # brackets on inner product of functions required to get exact symmetry
                tempM += quad_weights[q] * (basis_functions[q, i]' * basis_functions[q, j]) * determinants[q]
                tempS += quad_weights[q] * (curl_basis_functions[q, i]' * curl_basis_functions[q, j]) * determinants[q]
            end
            M[i, j, t] = tempM
            S[i, j, t] = tempS
        end
    end
    return M, S
end

"""
Assemble a global mass matrix:

M(V)_ij = ∫ coefficients[V] ϕ_i ⋅ ϕ_j dV
"""
function massMatrix(
    space::SimpleFiniteElementSpace{F} where {F<:HCurlElement},
    ; coefficients::Vector{Tv}=Float64[],
) where {Tv<:Number}

    local_functions = space.local_space.functions
    (_, num_functions, num_points, num_orientations) = size(local_functions)

    # convert to arrays of SVector
    local_functions = reinterpret(SVector{3,Float64}, local_functions)
    local_functions = reshape(local_functions, (num_functions, num_points, num_orientations))

    # get quad points and weights
    quad_weights = space.local_space.quad_weights
    quad_points = space.local_space.quad_points

    # temporary array to reduce allocations in loop
    basis_functions = zeros(SVector{3,Float64}, num_functions, num_points)

    # coo lists
    M_coo = Tv[]
    I_coo = Int64[]
    J_coo = Int64[]

    for (t, tag) in enumerate(space.mesh.elements.tags)
        # get orientation for this element
        orient = 1 + gmsh.model.mesh.getBasisFunctionsOrientationForElement(tag, "HcurlLegendre")

        # get jacobians for this element, convert to SMatrix
        jacobians, determinants, _ = gmsh.model.mesh.getJacobian(tag, quad_points)
        jacobians = reinterpret(SMatrix{3,3,Float64,9}, jacobians)

        for f = 1:num_functions, p = 1:num_points
            basis_functions[f, p] = jacobians[p]' \ local_functions[f, p, orient]
        end

        for j = 1:num_functions, i = 1:num_functions
            tempM::Tv = 0
            for q = 1:num_points
                # brackets on inner product of functions required to get exact symmetry
                tempM += quad_weights[q] * (basis_functions[i, q]' * basis_functions[j, q]) * determinants[q]
            end

            if length(coefficients) > 0
                tempM = tempM * coefficients[t]
            end

            push!(M_coo, tempM)
            push!(I_coo, space.dof_map[i, t])
            push!(J_coo, space.dof_map[j, t])
        end
    end
    N = maximum(space.dof_map)
    M = sparse(I_coo, J_coo, M_coo, N, N)
    return M
end

"""
Assemble a global curl matrix:

S(V)_ij = ∫ coefficients[V] curl(ϕ_i) ⋅ curl(ϕ_j) dV
"""
function curlMatrix(
    space::SimpleFiniteElementSpace{F} where {F<:HCurlElement},
    ; coefficients::Vector{Tv}=Float64[]
) where {Tv<:Number}

    curl_functions = space.local_space.curl_functions
    (_, num_functions, num_points, num_orientations) = size(curl_functions)

    # convert to arrays of SVector
    curl_functions = reinterpret(SVector{3,Float64}, curl_functions)
    curl_functions = reshape(curl_functions, (num_functions, num_points, num_orientations))

    # get quad points and weights
    quad_weights = space.local_space.quad_weights
    quad_points = space.local_space.quad_points

    # temporary array to reduce allocations in loop
    curl_basis_functions = zeros(SVector{3,Float64}, num_functions, num_points)

    # coo lists
    S_coo = eltype(coefficients)[]
    I_coo = Int64[]
    J_coo = Int64[]

    for (t, tag) in enumerate(space.mesh.elements.tags)
        # get orientation for this element
        orient = 1 + gmsh.model.mesh.getBasisFunctionsOrientationForElement(tag, "HcurlLegendre")

        # get jacobians for this element, convert to SMatrix
        jacobians, determinants, _ = gmsh.model.mesh.getJacobian(tag, quad_points)
        jacobians = reinterpret(SMatrix{3,3,Float64,9}, jacobians)

        for f = 1:num_functions, p = 1:num_points
            curl_basis_functions[f, p] = (jacobians[p] * curl_functions[f, p, orient]) / determinants[p]
        end

        for j = 1:num_functions, i = 1:num_functions
            tempS::Tv = 0
            for q = 1:num_points
                # brackets on inner product of functions required to get exact symmetry
                tempS += quad_weights[q] * (curl_basis_functions[i, q]' * curl_basis_functions[j, q]) * determinants[q]
            end

            if length(coefficients) > 0
                tempS = tempS * coefficients[t]
            end

            push!(S_coo, tempS)
            push!(I_coo, space.dof_map[i, t])
            push!(J_coo, space.dof_map[j, t])
        end
    end

    N = maximum(space.dof_map)
    S = sparse(I_coo, J_coo, S_coo, N, N)
    return S
end


"""
Integrate boundary physics so they can be applied to a global system of equations.
Potentially creates a matrix to be added to global matrix, and potentially a right hand side vector.

For Frequency domain problems: Y=jk
For time domain problems:
    Y=√(ϵ/μ) if using μ for curl part, or
    Y=√(ϵ*μ0) if using μr for curl part

Note: we can handle second order ABC in this function,
if we allow either a vector of coefficients Β, 
or allow the function to scale the coefficients using the radius.


The Robin condition for the Electric field is:

Frequency Domain:

̂n × (μr^-1 ∇ × E) + jk0 / ηr ̂n × (̂n × E) = K_N

Time Domain:

̂n × (μ^-1 ∇ × E(t)) + Y ̂n × (̂n × d/dt E(t)) = K_N(t)

Where ̂n, A, and P are vectors.

The difference between the Frequency and Time-domain methods is that Frequency domain absorbs μ0 and a jω term into the admittance.
"""
function integrateRobinBoundary(
    space::SimpleFiniteElementSpace{F},
    physics_tag::Integer,
    Y::Union{Float64,ComplexF64},
) where {F<:HCurlElement}
    dim = space.mesh.properties.dim
    bc_dim = dim - 1

    # reference map for nhat
    if F == HCurlTriangle
        nhat_reference = SA_F64[0, 1, 0]
    elseif F == HCurlTetrahedron
        nhat_reference = SA_F64[0, 0, 1]
    else
        throw("Function space: $F doesn't support robin boundary condition.")
    end

    # local boundary element
    boundary_element = HCurlElement(bc_dim, space.local_space.order)
    _, num_functions, num_points, num_orientations = size(boundary_element.functions)

    local_functions = reinterpret(SVector{3,Float64}, boundary_element.functions)
    local_functions = reshape(local_functions, (num_functions, num_points, num_orientations))

    # get the actual elements and degrees of freedom for boundary functions
    key = (bc_dim, physics_tag)
    boundary_tags = space.mesh.physics_map[key]
    global_dof = space.physical_group_mapping[key]

    # arrays for storing intermediary values
    nhat = Array{SVector{3,Float64}}(undef, num_points)
    basis_functions = Array{SVector{3,Float64}}(undef, num_points, num_functions)

    K_coo = typeof(Y)[]
    I_coo = []
    J_coo = []
    for (t, tag) in enumerate(boundary_tags)
        # orientation
        orient = 1 + gmsh.model.mesh.getBasisFunctionsOrientationForElement(tag, "HcurlLegendre")

        jacobians, determinants, _ = gmsh.model.mesh.getJacobian(tag, boundary_element.quad_points)
        jacobians = reinterpret(SMatrix{3,3,Float64,9}, jacobians)

        # map reference function to global element
        for p = 1:num_points
            nhat[p] = jacobians[p]' \ nhat_reference
        end
        for p = 1:num_points, f = 1:num_functions
            basis_functions[p, f] = cross(nhat[p], jacobians[p]' \ local_functions[f, p, orient])
        end

        # integrate ABC
        for j = 1:num_functions, i = 1:num_functions
            temp = 0.0
            for p = 1:num_points
                temp += boundary_element.quad_weights[p] * (basis_functions[p, i]' * basis_functions[p, j]) * determinants[p]
            end
            temp = temp * Y

            push!(K_coo, temp)
            push!(I_coo, global_dof[i, t])
            push!(J_coo, global_dof[j, t])
        end
    end

    N = maximum(space.dof_map)
    Krobin = sparse(I_coo, J_coo, K_coo, N, N)
    return Krobin
end

"""
Outputs need to be projected in a direction.
"""
function integratePointSource(
    space::SimpleFiniteElementSpace{F} where {F<:HCurlElement},
    x::Float64, y::Float64, z::Float64;
    curl=false
)

    properties = space.mesh.properties

    if curl == false
        basisName = "HcurlLegendre$(space.local_space.order)"
    else
        basisName = "CurlHcurlLegendre$(space.local_space.order)"
    end

    tag, type, _, u, v, w = gmsh.model.mesh.getElementByCoordinates(x, y, z)

    # if types don't match, then the integral won't work
    if type == properties.type
        # get the local functions first
        _, localFuncs, numOrientations = gmsh.model.mesh.getBasisFunctions(type, [u, v, w], basisName)
        numFunctions = length(localFuncs) ÷ numOrientations ÷ 3
        localFuncs = reshape(reinterpret(SVector{3,Float64}, localFuncs), (numFunctions, numOrientations))

        # now get global functions
        jac, det, _ = gmsh.model.mesh.getJacobian(tag, [u, v, w])
        jac = reinterpret(SMatrix{3,3,Float64,9}, jac)

        orient = 1 + gmsh.model.mesh.getBasisFunctionsOrientationForElement(tag, basisName)

        # map local to global
        if curl == false
            integral = [jac[1]' \ localFuncs[f, orient] for f = 1:numFunctions]
        else
            integral = [(jac[1] * localFuncs[f, orient]) / det[1] for f = 1:numFunctions]
        end

        # now get dof associated with integral
        idx = searchsortedfirst(space.mesh.elements.tags, tag)
        @assert space.mesh.elements.tags[idx] == tag
        dof_idxs = space.dof_map[:, idx]

        return integral, dof_idxs
    else
        printstyled("Could not find element\n", color=:red)
    end
    return SVector{3,Float64}[], Int64[]
end