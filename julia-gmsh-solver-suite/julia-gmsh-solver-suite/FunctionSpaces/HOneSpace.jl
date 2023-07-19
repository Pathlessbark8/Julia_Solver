using StaticArrays
using SparseArrays
using LinearAlgebra
using NodesAndModes
using StructArrays

"""
There is currently no integrator for Dirichlet or Neumann boundary conditions.
Zero Neumann is the default boundary condition for the H curl basis functions.
Dirichlet BC's are likely to be encountered for PEC, which we just set to zero, or for
scattered fields, for which we don't need to integrate as we'd already have the coefficients.
"""


"""
Assemble a GlobalFunctionSpace for a HOne triangle and a SimpleMesh mesh.

This function assigns local functions in each element a global degree of freedom to 
be used for the finite element method.
"""
function HOneSpace(mesh::SimpleMesh, local_space::HOneElement{T}) where {T}
        
    num_elements = numElements(mesh.elements)
    solution_order = local_space.order
    
    #our HOne space currently assumes a Lagrange basis 
    properties = mesh.properties
    geometric_order = properties.order
    element_type = properties.type

    #local nodes associated with nodal basis
    if T == Triangle
        r, s = nodes(Tri(), solution_order) #call to NodesAndModes
        t = 0 .* r
        element_type = gmsh.model.mesh.getElementType("triangle",geometric_order)
    elseif T == Tetrahedron
        r, s, t = nodes(Tet(), solution_order) #call to NodesAndModes
        element_type = gmsh.model.mesh.getElementType("tetrahedron",geometric_order)
    else
        @assert(0==1)
    end

    #map from [-1 1] to expected gmsh interval [0 1]
    r = (r .+ 1)./2.0                     
    s = (s .+ 1)./2.0
    t = (s .+ 1)./2.0
    n_local_points = length(r)

    #create local points array
    local_points = Array{Float64, 1}(undef, 3*n_local_points)
    local_points[1:3:end] = r
    local_points[2:3:end] = s
    local_points[3:3:end] = t
    local_points = reshape(local_points, 3, n_local_points) #storing as 3 x n_local_points

    # generate global points in each element so we can determine unique dofs
    n_comp_basis::Int64, basis_functions, n_orient::Int64 = gmsh.model.mesh.getBasisFunctions(element_type, local_points[:], "Lagrange$solution_order")
    n_basis = (((length(basis_functions) ÷ n_comp_basis) ÷ n_local_points) ÷ n_orient)
    basis_functions = reshape(basis_functions, n_comp_basis, n_basis, n_local_points, n_orient)

    # sort local coordinates so that basis are ordered properly (could this go in local dof map?)
    mask = abs.(basis_functions[1,:,:,1])
    mask = mask .> 1e-12
    nonzeros = findall(x -> x > 0, mask) #this sort will go by columns
    @assert(length(nonzeros) == n_local_points)
    permutearray = [nonzeros[i][1] for i in eachindex(nonzeros)]

    r = r[permutearray]
    s = s[permutearray]
    t = t[permutearray]

    local_points = Array{Float64, 1}(undef, 3*n_local_points)
    local_points[1:3:end] = r
    local_points[2:3:end] = s
    local_points[3:3:end] = t
    local_points = reshape(local_points, 3, n_local_points) #storing as 3 x n_local_points

    #global points at local points for each element - this comes from gmsh one element at a time
    global_points = Array{Float64,2}(undef,3,0)
    _, _, global_points = gmsh.model.mesh.getJacobians(element_type, local_points[:])
    n_global_points = length(global_points) ÷ 3
    global_points = reshape(global_points, 3, n_global_points) #storing as 3 x n_global_points
   
    #create a list of unique nodes from global points as follows:
    # round the points for comparison
    rounded_global_points = round.(global_points, digits=10)
    rounded_global_points = StructArray((rounded_global_points[1,:], rounded_global_points[2,:], rounded_global_points[3,:]))

    # sort the points
    sort_indexes = sortperm(rounded_global_points)
    sorted_rounded_global_points = rounded_global_points[sort_indexes]

    # find unique points (hence the rounding) and create a idx2unique map - just march through points until
    # some change is found, then you have a new point
    node_idx2unique = Array{Int64, 1}(undef, n_global_points)
    next_unique_index = 1
    for i in eachindex(sorted_rounded_global_points)
        if (i == 1)
            node_idx2unique[sort_indexes[i]] = next_unique_index
        else
            if (sorted_rounded_global_points[i] != sorted_rounded_global_points[i-1])
                next_unique_index+=1
            end
            node_idx2unique[sort_indexes[i]] = next_unique_index
        end
    end
    n_unique_nodes = next_unique_index #the number of degrees of freedom in this nodal basis
    
    global_dof_map = Array{Int64}(undef, (n_basis, num_elements))

    for iele in range(1,num_elements,step=1)
        global_ids = (1:n_basis) .+ (iele-1)*n_basis
        unique_ids = [node_idx2unique[global_id] for global_id in global_ids]
        global_dof_map[:,iele] = unique_ids
        #ele_idx2global[:,iele] = global_ids
    end

    # # create global dof list for physical groups of dimension lower than 2
    # this is only because SimpleFunctionSpace requires it - to be discussed
    globalPhysicsMap = Dict{Tuple{Int32,Int32},Array{Int32,2}}()

    SimpleFiniteElementSpace{HOneTriangle}(mesh, local_space, global_dof_map, globalPhysicsMap)  #last argument is globalPhysicsMap
end

