

#=------------------------------------
setupElementConstitutives()
-- assign constitutives based on tags and available element_physics_tags
-------------------------------------=#
function setupElementConstitutives(problem_type::AbstractProblemType, mesh::SimpleMesh, constitutives_table::Dict{Int64,Any}) 

    num_elements = length(mesh.elements.tags)

    ele_tag_to_idx = Dict{Int64, Int64}()

    next_index = 1
    for ele_tag in mesh.elements.tags
        ele_tag_to_idx[ele_tag] = next_index
        next_index = next_index + 1
    end

    element_complex_relative_permittivity = Vector{Complex{Float64}}(undef,num_elements)
    element_complex_relative_permeability = Vector{Complex{Float64}}(undef,num_elements)

    for key in keys(mesh.physics_map)
        dim = key[1]
        physics_tag = key[2]
        if dim == 2  #2d physics entry -- should be based on problem type
            constitutives = constitutives_table[physics_tag]
            element_tags = mesh.physics_map[key]
            for element_tag in element_tags
                element_index = ele_tag_to_idx[element_tag]
                element_complex_relative_permittivity[element_index] = constitutives.ε
                element_complex_relative_permeability[element_index] = constitutives.μ
            end
        end
    end

    element_constitutives = ComplexEMElementConstitutives(element_complex_relative_permittivity, element_complex_relative_permeability)
    return element_constitutives
end


function setupBoundaryConditions(problem_type::Union{EMTMProblemType, EMTEProblemType}, mesh::SimpleMesh, mesh_connectivity::SimpleMeshConnectivity, boundary_condition_table::Dict{Int64, AbstractEMBoundaryCondition})
    
    pec_conditions = Vector{EMElementBoundaryCondition}(undef, 0)
    pmc_conditions = Vector{EMElementBoundaryCondition}(undef, 0)
    abc_conditions = Vector{EMElementBoundaryCondition}(undef, 0)
    
    # need to be able to look up index of an element based on tag
    ele_tag_to_idx = Dict{Int64, Int64}()
    next_index = 1
    for ele_tag in mesh.elements.tags
        ele_tag_to_idx[ele_tag] = next_index
        next_index = next_index + 1
    end

    # need to be able to look up the physics tag of the boundary of an element
    boundary_ele_tag_to_phys = Dict{Int64, Int64}()
    for key in keys(mesh.physics_map)
        dim = key[1]
        physics_tag = key[2]
        if dim == 1
            for element_tag in mesh.physics_map[key]
                boundary_ele_tag_to_phys[element_tag] = physics_tag   # some check should be implemented to check that the element_tag hasn't been previously set             
            end
        end
    end


    E2D2D = mesh_connectivity.ele2D_to_ele2D
    
    # count missing element neighbours
    n_missing_neighbours = 0
    for ele_tag in(mesh.elements.tags)  
        for iface in(1:3)
            if !haskey(E2D2D,(ele_tag,iface))
                n_missing_neighbours+= 1
            end
        end
    end

    println("Number of missing neighbours = ", n_missing_neighbours)

    #actual boundary condition setup is done here

    E2D1D = mesh_connectivity.ele2D_to_ele1D

    for connection in eachindex(E2D1D)

        element_2d_tag = connection[1]
        element_2d_index = ele_tag_to_idx[element_2d_tag]
        face_index = connection[2]
        element_1d_tag = E2D1D[connection]
        element_1d_physics_tag = boundary_ele_tag_to_phys[element_1d_tag]
        
        if haskey(boundary_condition_table, element_1d_physics_tag)
            boundary_condition = boundary_condition_table[element_1d_physics_tag]
            element_boundary_condition = EMElementBoundaryCondition(element_2d_tag, element_2d_index, face_index, boundary_condition, element_1d_tag, element_1d_physics_tag)
            if (typeof(boundary_condition) == ABCEMBoundaryCondition)
                push!(abc_conditions, element_boundary_condition)
            elseif (typeof(boundary_condition) == PECEMBoundaryCondition)
                push!(pec_conditions, element_boundary_condition)
            elseif (typeof(boundary_condition) == PMCEMBoundaryCondition)
                push!(pmc_conditions, element_boundary_condition)
            end
        else
            println("Could not find a physics tag corresponding to ", element_1d_physics_tag)
            @assert(0==1)
        end

    end

    println("Found ", length(abc_conditions), " ABC boundary conditions.")
    println("Found ", length(pec_conditions), " PEC boundary conditions.")
    println("Found ", length(pmc_conditions), " PMC boundary conditions.")
    # println("Found ", length(imp_conditions), " IMP boundary conditions.")


    element_boundary_conditions = vcat(pec_conditions, pmc_conditions, abc_conditions)#, imp_conditions)
    return element_boundary_conditions
end
