
#=------------------------------------
evaluateSources()
-- Evaluate source terms for rhsB
-------------------------------------=#

function evaluateSources2DTM(mesh, nodal_basis, element_constitutives, frequency)
    
    n_dof = nodal_basis.n_unique_global_points
    n_sources = 1

    B = zeros(n_dof, n_sources)
    B[3,1] = 1.0
    return B
end

function evaluateSources2DTE(mesh, edge_basis, element_constitutives, frequency)
    
    n_basis_per_edge = edge_basis.n_basis ÷ length(edge_basis.element_edge_tags[:, 1]) #we assume the same number of basis functions per edge for now
    n_dof = length(edge_basis.unique_edge_tags)*n_basis_per_edge
    n_sources = 1
    
    B = zeros(n_dof, n_sources)
    B[300,1] = 1.0
    return B
end


function evaluateSources3D(mesh, edge_basis, element_constitutives, frequency)
    
    n_basis_per_edge = edge_basis.n_basis ÷ length(edge_basis.element_edge_tags[:, 1]) #we assume the same number of basis functions per edge for now
    n_dof = length(edge_basis.unique_edge_tags)*n_basis_per_edge
    n_sources = 1
    
    B = zeros(n_dof, n_sources)
    B[600,1] = 1.0
    return B
end
