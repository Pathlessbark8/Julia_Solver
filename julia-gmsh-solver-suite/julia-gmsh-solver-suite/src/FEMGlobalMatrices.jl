using Parameters
using SparseArrays

#=------------------------------------
COO Struct
-------------------------------------=#
@with_kw struct COOList
    rows::Array{Float64, 1}
    cols::Array{Float64, 1}
    vals::Array{Float64, 1}
end

#=------------------------------------
createFEMCOOLists2DTM()
-- Create coordinate lists for FEM global matrices - 2D TM Problems
-------------------------------------=#
function createFEMCOOLists2DTM(basis::NodalBasis, element_matrices::ElementMatrices2DTM, constitutives::ComplexConstitutives, frequency::Float64)
    
    COO_rows = Array{Float64,1}(undef,0)
    COO_cols = Array{Float64,1}(undef,0)
    
    K_COO_vals = Array{Float64,1}(undef,0)
    M_COO_vals = Array{Float64,1}(undef,0)

    ω = (2*π*frequency)
    ω2 = ω*ω

    for iele in range(1, basis.n_elements)
        ele_unique_points = basis.ele_idx2unique[:, iele]
        n_points = length(ele_unique_points);
        row_idxs = repeat(ele_unique_points,1,n_points)
        col_idxs = row_idxs';
        row_idxs = reshape(row_idxs, n_points*n_points)        
        col_idxs = reshape(col_idxs, n_points*n_points)

        M = element_matrices.M[:,:,iele]
        Sx = element_matrices.Sx[:,:,iele]
        Sy = element_matrices.Sy[:,:,iele]

        k2 = constitutives.element_ε_r[iele]*constitutives.element_μ_r[iele]*ε0*μ0*ω2

        K = k2*M - Sx - Sy 
        K_vals = reshape(K, n_points*n_points)
        M_vals = reshape(M, n_points*n_points)
        
        append!(COO_rows, row_idxs)
        append!(COO_cols, col_idxs)
        append!(K_COO_vals, K_vals)
        append!(M_COO_vals, M_vals)

    end

    K_coolist = COOList(COO_rows, COO_cols, K_COO_vals)
    M_coolist = COOList(COO_rows, COO_cols, M_COO_vals)
    
    return K_coolist, M_coolist

end

#=------------------------------------
createFEMCOOLists2DTE()
-- Create coordinate lists for FEM global matrices - 2D TE Problems
-------------------------------------=#
function createFEMCOOLists2DTE(basis::EdgeBasis, element_matrices::ElementMatrices2DTE, constitutives::ComplexConstitutives, frequency::Float64)
    
    COO_rows = Array{Float64,1}(undef,0)
    COO_cols = Array{Float64,1}(undef,0)
    
    K_COO_vals = Array{Float64,1}(undef,0)
    M_COO_vals = Array{Float64,1}(undef,0)

    ω = (2*π*frequency)
    ω2 = ω*ω

    n_basis_per_edge = basis.n_basis/length(basis.element_edge_tags[:, 1]) #we assume the same number of basis functions per edge for now

    for iele in range(1, basis.n_elements)
        ele_dof_indexes = Vector{Int64}(undef,0)
        
        for iedge in eachindex(basis.element_edge_tags[:, iele])
            
            edge_tag = basis.element_edge_tags[iedge, iele]
            dof_offset = (edge_tag - 1)*n_basis_per_edge
            dof_for_edge = Vector{Int64}(dof_offset .+ (1:n_basis_per_edge))
            ele_dof_indexes = cat(ele_dof_indexes, dof_for_edge, dims = 1)
        end
        @assert(length(ele_dof_indexes) == basis.n_basis)            
        row_idxs = repeat(ele_dof_indexes,1,basis.n_basis)
        col_idxs = row_idxs';
        row_idxs = reshape(row_idxs, basis.n_basis*basis.n_basis)        
        col_idxs = reshape(col_idxs, basis.n_basis*basis.n_basis)

        Mx = element_matrices.Mx[:,:,iele]
        My = element_matrices.My[:,:,iele]
        S = element_matrices.S[:,:,iele]

        k2 = constitutives.element_ε_r[iele]*constitutives.element_μ_r[iele]*ε0*μ0*ω2

        K = -k2*(Mx + My) + S
        K_vals = reshape(K, basis.n_basis*basis.n_basis)
        M_vals = reshape(Mx + My, basis.n_basis*basis.n_basis)
        
        append!(COO_rows, row_idxs)
        append!(COO_cols, col_idxs)
        append!(K_COO_vals, K_vals)
        append!(M_COO_vals, M_vals)

    end

    K_coolist = COOList(COO_rows, COO_cols, K_COO_vals)
    M_coolist = COOList(COO_rows, COO_cols, M_COO_vals)
    
    return K_coolist, M_coolist

end


#=------------------------------------
createFEMCOOLists3D()
-- Create coordinate lists for FEM global matrices - 3D Problems
-------------------------------------=#
function createFEMCOOLists3D(basis::EdgeBasis, element_matrices::ElementMatrices3D, constitutives::ComplexConstitutives, frequency::Float64)
    
    COO_rows = Array{Float64,1}(undef,0)
    COO_cols = Array{Float64,1}(undef,0)
    
    K_COO_vals = Array{Float64,1}(undef,0)
    M_COO_vals = Array{Float64,1}(undef,0)

    ω = (2*π*frequency)
    ω2 = ω*ω

    n_basis_per_edge = basis.n_basis/length(basis.element_edge_tags[:, 1]) #we assume the same number of basis functions per edge for now

    for iele in range(1, basis.n_elements)
        ele_dof_indexes = Vector{Int64}(undef,0)
        
        for iedge in eachindex(basis.element_edge_tags[:, iele])
            
            edge_tag = basis.element_edge_tags[iedge, iele]
            dof_offset = (edge_tag - 1)*n_basis_per_edge
            dof_for_edge = Vector{Int64}(dof_offset .+ (1:n_basis_per_edge))
            ele_dof_indexes = cat(ele_dof_indexes, dof_for_edge, dims = 1)
        end
        @assert(length(ele_dof_indexes) == basis.n_basis)            
        row_idxs = repeat(ele_dof_indexes,1,basis.n_basis)
        col_idxs = row_idxs';
        row_idxs = reshape(row_idxs, basis.n_basis*basis.n_basis)        
        col_idxs = reshape(col_idxs, basis.n_basis*basis.n_basis)

        Mx = element_matrices.Mx[:,:,iele]
        My = element_matrices.My[:,:,iele]
        Mz = element_matrices.Mz[:,:,iele]
        Sx = element_matrices.Sx[:,:,iele]
        Sy = element_matrices.Sy[:,:,iele]
        Sz = element_matrices.Sz[:,:,iele]

        k2 = constitutives.element_ε_r[iele]*constitutives.element_μ_r[iele]*ε0*μ0*ω2

        K = -k2*(Mx + My + Mz) + Sx + Sy + Sz
        K_vals = reshape(K, basis.n_basis*basis.n_basis)
        M_vals = reshape(Mx + My + Mz, basis.n_basis*basis.n_basis)
        
        append!(COO_rows, row_idxs)
        append!(COO_cols, col_idxs)
        append!(K_COO_vals, K_vals)
        append!(M_COO_vals, M_vals)

    end

    K_coolist = COOList(COO_rows, COO_cols, K_COO_vals)
    M_coolist = COOList(COO_rows, COO_cols, M_COO_vals)
    
    return K_coolist, M_coolist

end



#=------------------------------------
buildFEMGlobalMatrices2DTM()
-- Build the Global Matrix from the COO Lists
-------------------------------------=#
function buildFEMGlobalMatrices2DTM(K_coolist::COOList, M_coolist::COOList)
    K = sparse(K_coolist.rows, K_coolist.cols, K_coolist.vals);
    M = sparse(M_coolist.rows, M_coolist.cols, M_coolist.vals);
    return K, M
end


#=------------------------------------
buildFEMGlobalMatrices2DTE()
-- Build the Global Matrix from the COO Lists
-------------------------------------=#
function buildFEMGlobalMatrices2DTE(K_coolist::COOList, M_coolist::COOList)
    K = sparse(K_coolist.rows, K_coolist.cols, K_coolist.vals);
    M = sparse(M_coolist.rows, M_coolist.cols, M_coolist.vals);
    return K, M
end

#=------------------------------------
buildFEMGlobalMatrices3D()
-- Build the Global Matrix from the COO Lists
-------------------------------------=#
function buildFEMGlobalMatrices3D(K_coolist::COOList, M_coolist::COOList)
    K = sparse(K_coolist.rows, K_coolist.cols, K_coolist.vals);
    M = sparse(M_coolist.rows, M_coolist.cols, M_coolist.vals);
    return K, M
end
