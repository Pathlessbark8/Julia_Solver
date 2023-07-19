using LinearAlgebra
using StaticArrays
using SparseArrays

#=------------------------------------
Element Matrices Struct
-------------------------------------=#

struct EMTMElementVolumeMatrices <: AbstractElementVolumeMatrices
    B::Array{Float64,3}                     #nbasis x nquad_points x nelements
    M::Array{Float64,3}                     #nbasis x nbasis x nelements
    Sx::Array{Float64,3}                    #nbasis x nbasis x nelements
    Sy::Array{Float64,3}                    #nbasis x nbasis x nelements
end

struct EMTMBoundaryFluxMatrices <: AbstractBoundaryFluxMatrices
    F_minus::Array{Complex{Float64},3}
    ele_indexes::Vector{Int64}
    #Fplus::Array{Float64,3}
    #Frhs::Array{Float64,3}
end


function elementVolumeMatrices(problem_type::EMTMProblemType, space::SimpleFiniteElementSpace{HOneTriangle})
    
    local_functions = space.local_space.functions[1,:,:,:] #only one component for scalar functions
    gradient_functions = space.local_space.gradient_functions
    (_, num_functions, num_points, num_orientations) = size(gradient_functions)

    # convert to arrays of SVector
    gradient_functions = reinterpret(SVector{3,Float64}, gradient_functions)
    gradient_functions = reshape(gradient_functions, (num_functions, num_points, num_orientations))

    # get quad points and weights
    quad_weights = space.local_space.quad_weights
    quad_points = space.local_space.quad_points
    num_quad_points = length(quad_weights)
    #println("Num points = ", num_points, " num quad points = ", num_quad_points)
    @assert(num_points == num_quad_points)

    # storage for M and S
    B = zeros(num_functions, num_quad_points, length(space.mesh.elements.tags))
    M = zeros(num_functions, num_functions, length(space.mesh.elements.tags))
    Sx = zeros(num_functions, num_functions, length(space.mesh.elements.tags))
    Sy = zeros(num_functions, num_functions, length(space.mesh.elements.tags))

    #dBx, dBy, dBz #these would be useful

    # temporary arrays
    basis_functions = zeros(Float64, num_points, num_functions)
    gradient_basis_functions = zeros(SVector{3,Float64}, num_points, num_functions)

    for (t, tag) in enumerate(space.mesh.elements.tags)
        # get orientation for this element
        orient = 1 + gmsh.model.mesh.getBasisFunctionsOrientationForElement(tag, "Lagrange")

        # get jacobians for this element, convert to SMatrix
        jacobians, determinants, _ = gmsh.model.mesh.getJacobian(tag, quad_points)
        jacobians = reinterpret(SMatrix{3,3,Float64,9}, jacobians)

        for i = 1:num_points, f = 1:num_functions
            basis_functions[i, f] = local_functions[f, i, orient] # not fantastic for efficiency but done for now to be consistent with HCurl formulations. orient should be 1 and there is nothing to do for the value of the basis function
            gradient_basis_functions[i, f] = jacobians[i]' \ gradient_functions[f, i, orient]
        end

        for q = 1:num_quad_points, i = 1:num_functions
            B[i, q, t] = basis_functions[q,i]*determinants[q]*quad_weights[q]
        end

        for j = 1:num_functions, i = 1:num_functions
            tempM = 0
            tempSx = 0
            tempSy = 0

            for q = 1:num_points
                # brackets on inner product of functions required to get exact symmetry
                tempM += quad_weights[q] * (basis_functions[q, i]' * basis_functions[q, j]) * determinants[q]
                tempSx += quad_weights[q] * (gradient_basis_functions[q, i][1] * gradient_basis_functions[q, j][1]) * determinants[q]
                tempSy += quad_weights[q] * (gradient_basis_functions[q, i][2] * gradient_basis_functions[q, j][2]) * determinants[q]
            
            end

            M[i, j, t] = tempM
            Sx[i, j, t] = tempSx
            Sy[i, j, t] = tempSy
        
        end
    end

    element_matrices = EMTMElementVolumeMatrices(B, M, Sx, Sy)
    return element_matrices
end

function elementFluxMatrices(problem_type::EMTMProblemType, space::SimpleFiniteElementSpace{HOneTriangle}, boundary_conditions::Vector{EMElementBoundaryCondition}, constitutives::ComplexEMElementConstitutives, frequency::Float64)

    # Note: this function is currently written assuming element type triangle0 - it may or may not work otherwise
    
    # get basic setup information
    element_type = space.mesh.properties.type
    solution_order = space.local_space.order
    integration_order = solution_order*solution_order + 1
    num_comp, num_basis, num_quad_points, num_orient = size(space.local_space.functions)

    # # get edges 
    # edge_nodes = convert.(Int64,gmsh.model.mesh.getElementEdgeNodes(element_type))   
    # element_edge_tags, _ = gmsh.model.mesh.getEdges(edge_nodes)   # 3 * number of elements for triangles

    # # create gmsh elements for edges (this could be unified in one)
    # new_mesh_entity = gmsh.model.addDiscreteEntity(1) #add a 1D mesh entity for the 1D edge elements
    # line_element_tag_offset = gmsh.model.mesh.getMaxElementTag() # new element tags (which are unique in the mesh) must be offset from existing tags
    
    # unique_line_tags = Set()    # keep track of added lines
    # tags_for_lines = []         # list of new line tags
    # nodes_for_lines = []        # list of nodes for lines
    
    # # loop over element edges to create a unique set of edges to add as lines
    # for iedge in eachindex(element_edge_tags)  
    #     if !(in(element_edge_tags[iedge] , unique_line_tags))   # if we haven't added this line yet
    #         push!(unique_line_tags, element_edge_tags[iedge])   # add the line tag to the set 
    #         push!(tags_for_lines, element_edge_tags[iedge] + line_element_tag_offset) # add the line tag to the list, offset appropriately
    #         append!(nodes_for_lines, edge_nodes[2*(iedge-1) .+ (1:2)]) # add the node indexes in the mesh that make up the line    
    #     end
    # end

    line_type = gmsh.model.mesh.getElementType("line",1)

    # get a local integration rule on line elements for evaluating integrals
    line_quad_rule = gmsh.model.mesh.getIntegrationPoints(line_type,"Gauss$integration_order")
    line_n_quad_points = length(line_quad_rule[1]) ÷ 3
    line_quad_points = line_quad_rule[1]
    line_quad_coordinate = 0.5 .* (line_quad_points[1:3:end] .+ 1) #only the first coordinate matters
    line_quad_weights = line_quad_rule[2]

    W = diagm(line_quad_weights)
    nxdBx = Array{Complex{Float64}, 2}(undef, num_basis, line_n_quad_points) #derivative of basis in x
    nydBy = Array{Complex{Float64}, 2}(undef, num_basis, line_n_quad_points) #derivative of basis in y
    dBx = Array{Complex{Float64}, 2}(undef, num_basis, line_n_quad_points) #derivative of basis in x
    dBy = Array{Complex{Float64}, 2}(undef, num_basis, line_n_quad_points) #derivative of basis in y
    
    nxnxE = Array{Complex{Float64}, 2}(undef, num_basis, line_n_quad_points)
    F = Array{Complex{Float64}, 3}(undef,  num_basis, num_basis, length(boundary_conditions)) #nbasis x nbasis x nelements
    F2 = Array{Complex{Float64}, 3}(undef,  num_basis, num_basis, length(boundary_conditions)) #nbasis x nbasis x nelements
    ele_indexes = Vector{Int64}(undef, length(boundary_conditions))
    # process boundary conditions

    nxvec = Vector{Float64}(undef,0)
    nyvec = Vector{Float64}(undef,0)
    xvec = Vector{Float64}(undef,0)
    yvec = Vector{Float64}(undef,0)

    for ibd in eachindex(boundary_conditions)
        
        boundary = boundary_conditions[ibd]

        # get element properties (tag, index, face) from boundary -- TODO: check if this will properly support PEC/PMC sheets (internal boundaries)
        element_tag = boundary.ele_tag
        element_index = boundary.ele_index
        element_face_index = boundary.ele_face_index
        boundary_type = boundary.boundary_type
        boundary_physics_tag = boundary.boundary_physics_tag 
        boundary_element_tag = boundary.boundary_element_tag

        # map line quadrature points to points inside the element -- this mapping is element dependent
        # and a triangle is assumed here
        if element_face_index == 1
            element_edge_quad_points = hcat(line_quad_coordinate, 0*line_quad_coordinate, 0*line_quad_coordinate)' # edge 1 = (u, 0)
        elseif element_face_index == 2
            element_edge_quad_points = hcat(line_quad_coordinate, 1 .- line_quad_coordinate, 0*line_quad_coordinate)' # edge 2 = (u, 1-u)
        elseif element_face_index == 3
            element_edge_quad_points = hcat(0*line_quad_coordinate, line_quad_coordinate, 0* line_quad_coordinate)' #edge 3 = (0,v)
            @assert(false,"Unexpected face index")
        else
            @assert(false,"Unexpected face index")
        end
        
        # get the 1d line jacobian
        line_jacobian, line_determinant, line_global_edge_points = gmsh.model.mesh.getJacobian(boundary_element_tag, line_quad_points[:])
        line_jacobian = reshape(line_jacobian, 3, 3, line_n_quad_points)

        # get the 2d element jacobian
        ele_jacobian, _, ele_global_edge_points = gmsh.model.mesh.getJacobian(element_tag, element_edge_quad_points[:])
        ele_jacobian = reshape(ele_jacobian, 3, 3, line_n_quad_points)

        @assert(norm(ele_global_edge_points - line_global_edge_points) < 1e-12)

        # get the outward normal to the boundary
        surface_normal = ele_jacobian[:,3,:]
        edge_tangent = line_jacobian[:,1,:] # this is +/- 1*element edge vector but we will adjust

        # calculate edge normal at each point along the edge
        nhat = [cross(edge_tangent[:,ipoint], surface_normal[:,ipoint]) for ipoint in 1:line_n_quad_points]
        nhat = [nhat[ipoint]/norm(nhat[ipoint]) for ipoint in 1:line_n_quad_points]

        for ipoint in eachindex(nhat)

            push!(nxvec, nhat[ipoint][1])
            push!(nyvec, nhat[ipoint][2])
            push!(xvec, ele_global_edge_points[3*(ipoint-1) + 1])
            push!(yvec, ele_global_edge_points[3*(ipoint-1) + 2])

        end

        # orient edge normal outwards
        _, _, ele_centroid = gmsh.model.mesh.getJacobian(element_tag, [1/3, 1/3, 0])
        ele_centroid_to_edge = ele_global_edge_points[1:3] - ele_centroid #just take the first point arbitrarily
        mesh_centroid_to_edge = ele_global_edge_points[1:3]/norm(ele_global_edge_points[1:3]) #not quite but close
        for ipoint in 1:line_n_quad_points
            if dot(ele_centroid_to_edge, nhat[ipoint]) < 0 # change projection if negative dot product
                nhat[ipoint] = nhat[ipoint] .* -1
            end
            # proj = dot(mesh_centroid_to_edge, nhat[ipoint])
            # println("Projection = ", proj)
        end


        # basis functions along the boundary
        n_comp_basis::Int64, basis_on_surface, n_orient::Int64 = gmsh.model.mesh.getBasisFunctions(element_type, element_edge_quad_points[:], "Lagrange$solution_order")
        n_basis_on_surface = (((length(basis_on_surface) ÷ n_comp_basis) ÷ line_n_quad_points) ÷ n_orient)
        basis_on_surface = reshape(basis_on_surface, n_comp_basis, n_basis_on_surface, line_n_quad_points, n_orient)
        B = basis_on_surface[1,:,:,1]

        # gradient of basis functions along the surface
        n_comp_grad_basis::Int64, grad_basis_on_surface, n_orient_grad::Int64 = gmsh.model.mesh.getBasisFunctions(element_type, element_edge_quad_points[:], string("GradLagrange$solution_order")) 
        grad_basis_on_surface = reshape(grad_basis_on_surface, n_comp_grad_basis, n_basis_on_surface, line_n_quad_points, n_orient_grad)
        

        #------------------------------------------------
        # Pull out constitutives at points along the edge
        #------------------------------------------------

        α = 1 # upwind parameter
        ω = 2*π*frequency
        ε_minus = [constitutives.element_ε_r[element_index]*ε0 for ipoint in 1:line_n_quad_points]
        μ_minus = [constitutives.element_μ_r[element_index]*μ0 for ipoint in 1:line_n_quad_points]
        Z_minus = [sqrt(μ_minus[ipoint]/ε_minus[ipoint]) for ipoint in 1:line_n_quad_points]
        Y_minus = [1/Z_minus[ipoint] for ipoint in 1:line_n_quad_points]
        
        if boundary_type == ABCEMBoundaryCondition()
            α_abc = 1
            α = α_abc
            c1 = -1im * ω .* (μ_minus .* α_abc) ./ (2 .* Z_minus)
            c2 = Z_minus ./ (2 .* Z_minus)

            Zm = Z_minus[1]
            Zp = Z_minus[1]
        elseif boundary_type == PECEMBoundaryCondition()
            α_pec = 1
            α = α_pec
            c1 = -1im * ω .* (μ_minus .* α_pec) ./ Z_minus
            c2 = 1 .* Z_minus ./ Z_minus
            Zm = Z_minus[1]
            Zp = 0
        elseif boundary_type == PMCEMBoundaryCondition()
            c1 = 0 .* Z_minus
            c2 = 0 .* Z_minus
        # elseif boundary_type == IMP::BoundaryConditionType
        #     Z_plus = boundary.Z
        #     Y_plus = 1/boundary.Z
        else
            @assert(false,"Unknown Boundary Condition Type")
        end

        #------------------------------------------------
        # For TM wave equation we need to integrate n dot gradE = -n x curlE = jωμ n x H
        # = -jωμ (-α/(Z^- + Z^+) n x n x E - Z^-/(Z^- + Z^+) n x H) for one-sided conditions
        # with n x H = -1/(jωμ) n x curlE = +1/(jωμ) n dot gradE
        #------------------------------------------------

        # for TM problems we need nx*dBx and ny*dBy
        for ipoint in 1:line_n_quad_points
            #jacobian^-T (at a point) * gradient of local basis = gradient in global coordinates (with a transpose or something)?
            #T = ele_jacobian[:,:,ipoint]'\grad_basis_on_surface[:,:,ipoint,1]
            T = ele_jacobian[:,:,ipoint]'\grad_basis_on_surface[:,:,ipoint,1]

            nxnxE[:,ipoint] =  c1[ipoint]*B[:,ipoint]
            nxdBx[:,ipoint] =  c2[ipoint]*nhat[ipoint][1]*T[1,:] #nx*dBx at each point
            nydBy[:,ipoint] =  c2[ipoint]*nhat[ipoint][2]*T[2,:] #ny*dBy at each point

            dBx[:,ipoint] = T[1,:]
            dBy[:,ipoint] = T[2,:]

            
        end
            

        J = diagm(line_determinant) #diagonal matrix of quad point jacobians for this element

        nx = nhat[1][1]
        ny = nhat[1][2]

        # calculate Fminus matrix for one-sided boundary conditions
        #F[:,:,ibd] = B*(J*W)*nxnxE' + B*(J*W)*(nxdBx' + nydBy')
        F[:,:,ibd] = B*(J*W)*transpose(nxnxE) + B*(J*W)*(transpose(nxdBx) + transpose(nydBy))
        F2[:,:,ibd] = -1im*ω*μ_minus[1]*α/(Zm + Zp)*(B*(J*W)*B') + Zm/(Zm + Zp)*(nx*B*(J*W)*dBx' + ny*B*(J*W)*dBy')

        @assert(norm(F[:,:,ibd] - F2[:,:,ibd]) < 1e-12)
        
        ele_indexes[ibd] = element_index

        # wonder if we should do these for "actual physics" and "background physics" at the same time        
            
        # might need some extra stuff for PEC/PMC for scattered field formulations
    
    end
       
    boundary_flux_matrices = EMTMBoundaryFluxMatrices(F,ele_indexes)

    return boundary_flux_matrices

end

#=------------------------------------
createFEMCOOList()
-- 2D TM Problems
-------------------------------------=#
function createFEMCOOList(problem_type::EMTMProblemType, space::SimpleFiniteElementSpace{HOneTriangle}, element_volume_matrices::EMTMElementVolumeMatrices, boundary_flux_matrices::EMTMBoundaryFluxMatrices, constitutives::ComplexEMElementConstitutives, frequency::Float64)
    
    K_COO_rows = Array{Int64,1}(undef,0)
    K_COO_cols = Array{Int64,1}(undef,0)    
    K_COO_vals = Array{ComplexF64,1}(undef,0)

    B_COO_rows = Array{Int64,1}(undef,0)
    B_COO_cols = Array{Int64,1}(undef,0)
    B_COO_vals = Array{ComplexF64,1}(undef,0)
    
    ω = (2*π*frequency)
    ω2 = ω*ω
    jωμ0 = 1.0im*ω*μ0
    
    n_elements = size(element_volume_matrices.M,3)
    num_comp, num_basis, num_quad_points, num_orient = size(space.local_space.functions)
    
    local_quad_point_indexes = 1:num_quad_points
    for iele in range(1, n_elements)
        ele_unique_points = space.dof_map[:, iele]        
        ele_global_quad_points = (iele-1)*num_quad_points .+ local_quad_point_indexes  # assumes a constant number of quad points on each element for the testing integral
        num_points = length(ele_unique_points)
        
        K_row_idxs = repeat(ele_unique_points,1,num_points)
        K_col_idxs = K_row_idxs'
        K_row_idxs = reshape(K_row_idxs, num_points*num_points)        
        K_col_idxs = reshape(K_col_idxs, num_points*num_points)
        B_row_idxs = repeat(ele_unique_points, 1, num_quad_points)
        B_col_idxs = repeat(ele_global_quad_points, num_points, 1)
        B_row_idxs = reshape(B_row_idxs, num_points*num_quad_points)
        B_col_idxs = reshape(B_col_idxs, num_points*num_quad_points)
        M = element_volume_matrices.M[:,:,iele]
        Sx = element_volume_matrices.Sx[:,:,iele]
        Sy = element_volume_matrices.Sy[:,:,iele]
        B = element_volume_matrices.B[:,:,iele]
        k2 = constitutives.element_ε_r[iele]*constitutives.element_μ_r[iele]*ε0*μ0*ω2
        K = k2*M - Sx - Sy 
        B *= jωμ0*constitutives.element_μ_r[iele]
        K_vals = reshape(K, num_points*num_points)
        B_vals = reshape(B, num_points*num_quad_points)
        append!(B_COO_rows, B_row_idxs)
        append!(B_COO_cols, B_col_idxs)
        append!(K_COO_rows, K_row_idxs)
        append!(K_COO_cols, K_col_idxs)
        
        append!(K_COO_vals, K_vals)
        append!(B_COO_vals, B_vals)
    
    end
    for iele in eachindex(boundary_flux_matrices.ele_indexes)
        ele_index = boundary_flux_matrices.ele_indexes[iele]
        ele_unique_points = space.dof_map[:, ele_index]        
        num_points = length(ele_unique_points)
        
        F_row_idxs = repeat(ele_unique_points,1,num_points)
        F_col_idxs = F_row_idxs'
        F_row_idxs = reshape(F_row_idxs, num_points*num_points)        
        F_col_idxs = reshape(F_col_idxs, num_points*num_points)
        F = boundary_flux_matrices.F_minus[:,:,iele]
        
        F_vals = reshape(F, num_points*num_points)
        append!(K_COO_rows, F_row_idxs)
        append!(K_COO_cols, F_col_idxs)
        append!(K_COO_vals, F_vals)
    end
    K_coolist = COOList(K_COO_rows, K_COO_cols, K_COO_vals)
    B_coolist = COOList(B_COO_rows, B_COO_cols, B_COO_vals)
     
    return K_coolist, B_coolist
end

#=------------------------------------
buildFEMGlobalMatrices()
-- Build the Global Matrix from the COO Lists
-------------------------------------=#
function buildFEMGlobalMatrices(K_coolist::COOList, B_coolist::COOList)
    K = sparse(K_coolist.rows, K_coolist.cols, K_coolist.vals);
    B = sparse(B_coolist.rows, B_coolist.cols, B_coolist.vals);
    return K, B
end


