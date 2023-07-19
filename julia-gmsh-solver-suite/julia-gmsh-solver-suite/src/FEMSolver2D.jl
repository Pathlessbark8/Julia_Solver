#=============================================================
2D TM Electromagnetics FEM Solver using the Gmsh API
==============================================================#

if isdefined(@__MODULE__, :EIL)
    
else
    include("EIL.jl")
end

using LinearAlgebra
using Plots
using Parameters
using NodesAndModes
using StructArrays
using .EIL

#=-------------------------------------------------------------------------------------
 Main Program
--------------------------------------------------------------------------------------=#

#=------------------------------------
 User Input (?)
------------------------------------=#
mesh_filename = "./PECCircle.msh"
problem_type = EIL.EMTMType()
basis_type = EIL.NodalBasisType()
solution_order = 1

frequency = 1e8

constitutives_table = Dict{Int64,Any}()
constitutives_table[3000] = SimpleComplexConstitutives("Freespace", 3000, EIL.simple_complex::EIL.ConstitutivesType, 1.0 + 0.0im, 1.0 + 0.0im)
constitutives_table[3001] = SimpleComplexConstitutives("Target", 3001, EIL.simple_complex::EIL.ConstitutivesType, 2.0 - 0.5im, 1.0 + 0.0im)
constitutives_table[3002] = SimpleComplexConstitutives("Other", 3002, EIL.simple_complex::EIL.ConstitutivesType, 2.0 + 0.0im, 1.0 + 0.0im)

boundary_condition_table = Dict{Int64, Any}()
boundary_condition_table[1000] = BoundaryCondition("PEC", 1000, EIL.PEC::EIL.BoundaryConditionType)
boundary_condition_table[1001] = BoundaryCondition("ABC", 1001, EIL.ABC::EIL.BoundaryConditionType)

#=------------------------------------
 Mesh Setup
------------------------------------=#

setupGmsh()
mesh = readMesh(mesh_filename)
mesh_connectivity = setupMeshConnectivity(mesh)

#=------------------------------------
 Create the basis functions
------------------------------------=#

local_basis = setupLocalBasis(mesh, solution_order, problem_type)
dof_map::EIL.AbstractDOFMap = setupDOFMap(local_basis.element_type, local_basis.solution_order, basis_type)

#=------------------------------------
 Construct local matrices (no constitutives)
------------------------------------=#

element_matrices = elementMassAndStiffnessMatrices(local_basis, problem_type)


# # element_matrices_2DTE::ElementMatrices2DTE = elementMassAndStiffnessMatrices2DTE(mesh_2D, edge_basis_2D) #for 2D TE
# # element_matrices_3D::ElementMatrices3D = elementMassAndStiffnessMatrices3D(mesh_3D, edge_basis_3D) #for 3D

# #=------------------------------------
#  Setup element constitutives
# ------------------------------------=#

# element_constitutives_2D = setupElementConstitutives(mesh_2D, nodal_basis_2D, constitutives_table) #these could be provided on nodes instead
# # element_constitutives_3D = setupElementConstitutives3D(mesh_3D, edge_basis_3D, constitutives_table) #these could be provided on nodes instead

# #=------------------------------------
#  Create matrix entries
# ------------------------------------=#

# K_coolist_2DTM::COOList, M_coolist_2DTM::COOList = createFEMCOOLists2DTM(nodal_basis_2D, element_matrices_2DTM, element_constitutives_2D, frequency) #for 2D TM
# # K_coolist_2DTE::COOList, M_coolist_2DTE::COOList = createFEMCOOLists2DTE(edge_basis_2D, element_matrices_2DTE, element_constitutives_2D, frequency) #for 2D TE
# # K_coolist_3D::COOList, M_coolist_3D::COOList = createFEMCOOLists3D(edge_basis_3D, element_matrices_3D, element_constitutives_3D, frequency) #for 3D

# #=------------------------------------
#  Build global matrix
# ------------------------------------=#

# K_2DTM, M_2DTM = buildFEMGlobalMatrices2DTM(K_coolist_2DTM, M_coolist_2DTM) #for 2D TM
# # K_2DTE, M_2DTE = buildFEMGlobalMatrices2DTE(K_coolist_2DTE, M_coolist_2DTE) #for 2D TE
# # K_3D, M_3D = buildFEMGlobalMatrices3D(K_coolist_3D, M_coolist_3D) #for 2D TE

# #=------------------------------------
#  Construct RHS
# ------------------------------------=#

# B_2DTM = evaluateSources2DTM(mesh_2D, nodal_basis_2D, element_constitutives_2D, frequency)
# # B_2DTE = evaluateSources2DTE(mesh_2D, edge_basis_2D, element_constitutives_2D, frequency)
# # B_3D = evaluateSources3D(mesh_3D, edge_basis_3D, element_constitutives_3D, frequency)

# #=------------------------------------
#  Fill and Factor Operator
# ------------------------------------=#

# lu_solver_2DTM = lu(K_2DTM)
# # lu_solver_2DTE = lu(K_2DTE)
# # lu_solver_3D = lu(K_3D)

# # s = Plots.spy(lu_solver_3D.L)
# # display(s)
# # s = Plots.spy(lu_solver_3D.U)
# # display(s)

# #=------------------------------------
#  Solve System
# ------------------------------------=#
# ω = 2*π*frequency
# jωμ = 1.0im*ω*μ0

# solution_2DTM = lu_solver_2DTM \ ((jωμ)*(M_2DTM*B_2DTM))
# # solution_2DTM = solution_2DTM[:,1]

# # solution_2DTE = lu_solver_2DTE \ ((-jωμ)*(M_2DTE*B_2DTE))
# # solution_2DTE = solution_2DTE[:,1] #edge coefficient solution

# # solution_3D = lu_solver_3D \ ((-jωμ)*(M_3D*B_3D))

# #=-----------------------------------#
#  Reorder 2D TM solution
# ------------------------------------=#
# # assuming a first-order solution, the nodes in the nodal basis
# # align with gmsh nodes, but we need to reorder them
# gmsh_solution_2DTM = 0*solution_2DTM
# gmsh_node_processed = zeros(Int64, length(solution_2DTM))
# for iele = 1:mesh_2D.elements[3].count
#     for inode = 1:mesh_2D.elements[3].nodes_per_element 
#         global_node_index = (iele-1)*3 + inode
#         gmsh_node_index = nodal_basis_2D.node_idx2gmsh[global_node_index]
#         gmsh_node = mesh_2D.nodes.coords[:, gmsh_node_index]
        
#         unique_node_index = nodal_basis_2D.node_idx2unique[global_node_index]
#         unique_node = nodal_basis_2D.global_points[:, global_node_index]
        
#         @assert(norm(gmsh_node - unique_node) < 1e-12)

#         if (gmsh_node_processed[gmsh_node_index] == 0)
#             gmsh_node_processed[gmsh_node_index] = 1
#             gmsh_solution_2DTM[gmsh_node_index] = solution_2DTM[unique_node_index]
#         end
#     end
# end

# # #=-----------------------------------#
# #  Evaluate 2D TE Solution at centroid of each element
# # ------------------------------------=#

# # triangle_element_type = mesh_2D.element_types[3]

# # solution_Ex_2DTE = Vector{Float64}(undef,0)
# # solution_Ey_2DTE = Vector{Float64}(undef,0)
# # n_basis_per_edge = edge_basis_2D.n_basis ÷ length(edge_basis_2D.element_edge_tags[:, 1]) #we assume the same number of basis functions per edge for now

# # local_centroid = [1.0/3.0, 1.0/3.0, 0]

# # global_centroids = Array{Float64, 3}(undef, 3, 1, 0)
# # jacobians_at_centroids = Array{Float64, 4}(undef, 3, 3, 1, 0)
# # determinants_at_centroids = Array{Float64, 2}(undef, 1, 0)    
# # barycentres = Array{Float64,2}(undef,3, 0)
# # #get basis functions at centroids of each element
# # for i in eachindex(mesh_2D.physical_groups)

# #     group_dim = mesh_2D.physical_groups[i][1]
# #     group_tag = mesh_2D.physical_groups[i][2]
# #     group_entities = gmsh.model.getEntitiesForPhysicalGroup(group_dim, group_tag)

# #     #loop over entries in physical group
# #     for entity in group_entities
        
# #         entity_element_types = gmsh.model.mesh.getElementTypes(group_dim, entity)
        
# #         for itype in eachindex(entity_element_types)

# #             if entity_element_types[itype] != mesh_2D.element_types[3]
# #                 continue #only interested in these elements for this basis
# #             end
            
# #             entity_barycentres = gmsh.model.mesh.getBarycenters(triangle_element_type,entity, false, true)
# #             entity_jacobians_at_centroids, entity_determinants_at_centroids, entity_global_centroids = gmsh.model.mesh.getJacobians(entity_element_types[itype], local_centroid, entity)
# #             n_entity_elements = (length(entity_global_centroids) ÷ 3) ÷ 1
            
# #             entity_barycentres = reshape(entity_barycentres,3, n_entity_elements)
# #             entity_global_centroids = reshape(entity_global_centroids, 3, 1, n_entity_elements)
# #             entity_determinants_at_centroids = reshape(entity_determinants_at_centroids, 1, n_entity_elements) #scalar for each point and each element
# #             entity_jacobians_at_centroids = reshape(entity_jacobians_at_centroids, 3, 3, 1, n_entity_elements) #3x3 matrix for each point and each element

# #             global barycentres = cat(barycentres, entity_barycentres, dims=2)
# #             global global_centroids = cat(global_centroids, entity_global_centroids, dims=3)
# #             global jacobians_at_centroids = cat(jacobians_at_centroids, entity_jacobians_at_centroids, dims=4)
# #             global determinants_at_centroids = cat(determinants_at_centroids, entity_determinants_at_centroids, dims=2)

# #         end #each type in entity
# #     end # group entities
# # end #physical groups

# # n_comp_basis::Int64, basis_functions_at_centroids, n_orient::Int64 = gmsh.model.mesh.getBasisFunctions(mesh_2D.element_types[3], local_centroid[:], "HcurlLegendre0")
# # n_basis = (((length(basis_functions_at_centroids) ÷ n_comp_basis) ÷ 1) ÷ n_orient)
# # basis_functions_at_centroids = reshape(basis_functions_at_centroids, n_comp_basis, n_basis, 1, n_orient)


# # # now use basis functions to reconstruct solution at the centroid
# # reconstruction_points = Vector{Float64}(undef,0)
# # solution_at_points_2DTE = Vector{Complex{Float64}}(undef,0) #vector solution at array of points

# # for iele in range(1,edge_basis_2D.n_elements)
    
# #     basis_orientation = edge_basis_2D.element_basis_orientations[iele]
# #     basis_orientation += 1 #gmsh to julia offset
    
# #     for ipoint in 1:1 #only one centroid

# #         solution_at_point = zeros(Float64,3)
        
# #         reconstruction_point = global_centroids[:,ipoint, iele] 
                     
# #         global reconstruction_points = cat(reconstruction_points, reconstruction_point, dims=1)
# #         #solution_at_point = barycentres[:,iele]
# #         for iedge in eachindex(edge_basis_2D.element_edge_tags[:, iele])
            
# #             edge_tag = edge_basis_2D.element_edge_tags[iedge, iele]
            
# #             for ibasis_per_edge in 1:n_basis_per_edge
# #                 basis_index = (iedge-1)*n_basis_per_edge + ibasis_per_edge
# #                 @assert basis_index == iedge
# #                 # if (basis_index != mod(edge_tag, 3) + 1)
# #                 #     println("edge tag = ", edge_tag)
# #                 #     println("basis_index = ", basis_index)
# #                 #     @assert basis_index == mod(edge_tag-1, 3) + 1
# #                 # end
# #                 jacobian = jacobians_at_centroids[:, :, ipoint, iele]
# #                 basis_function = (jacobian')\basis_functions_at_centroids[:,basis_index,ipoint,basis_orientation]  

# #                 coefficient_index = (edge_tag - 1)*n_basis_per_edge + ibasis_per_edge;

# #                 @assert coefficient_index == edge_tag

# #                 solution_at_point += solution_2DTE[coefficient_index]*basis_function;
# #                 #solution_at_point = global_centroids[:,1,iele]
                

# #                 #println("adding basis contribution")
# #             end

# #         end
# #         if abs(imag(solution_at_point[1])) > 10 
# #             solution_at_point[1] = 0.0
# #         end
# #         if abs(imag(solution_at_point[2])) > 10 
# #             solution_at_point[2] = 0.0
# #         end
# #         if abs(imag(solution_at_point[3])) > 10 
# #             solution_at_point[3] = 0.0
# #         end


# #         global solution_at_points_2DTE = cat(solution_at_points_2DTE, solution_at_point, dims=1)
# #         #println("total solution at point done")
# #     end
# # end
# # solution_at_points_2DTE = reshape(solution_at_points_2DTE, 3, length(solution_at_points_2DTE[:]) ÷ 3)

# #=-----------------------------------#
#  Evaluate 3D Solution at centroid of each element
# ------------------------------------=#

# # tetrahedron_element_type = mesh_3D.element_types[4]

# # n_basis_per_edge = edge_basis_3D.n_basis ÷ length(edge_basis_3D.element_edge_tags[:, 1]) #we assume the same number of basis functions per edge for now

# # local_centroid = [1.0/3.0, 1.0/3.0, 1.0/3.0]

# # global_centroids = Array{Float64, 3}(undef, 3, 1, 0)
# # jacobians_at_centroids = Array{Float64, 4}(undef, 3, 3, 1, 0)
# # determinants_at_centroids = Array{Float64, 2}(undef, 1, 0)    
# # barycentres = Array{Float64,2}(undef,3, 0)
# # #get basis functions at centroids of each element
# # for i in eachindex(mesh_3D.physical_groups)

# #     group_dim = mesh_3D.physical_groups[i][1]
# #     group_tag = mesh_3D.physical_groups[i][2]
# #     group_entities = gmsh.model.getEntitiesForPhysicalGroup(group_dim, group_tag)

# #     #loop over entries in physical group
# #     for entity in group_entities
        
# #         entity_element_types = gmsh.model.mesh.getElementTypes(group_dim, entity)
        
# #         for itype in eachindex(entity_element_types)

# #             if entity_element_types[itype] != tetrahedron_element_type
# #                 continue #only interested in these elements for this basis
# #             end
            
# #             entity_barycentres = gmsh.model.mesh.getBarycenters(tetrahedron_element_type,entity, false, true)
# #             entity_jacobians_at_centroids, entity_determinants_at_centroids, entity_global_centroids = gmsh.model.mesh.getJacobians(entity_element_types[itype], local_centroid, entity)
# #             n_entity_elements = (length(entity_global_centroids) ÷ 3) ÷ 1
            
# #             entity_barycentres = reshape(entity_barycentres,3, n_entity_elements)
# #             entity_global_centroids = reshape(entity_global_centroids, 3, 1, n_entity_elements)
# #             entity_determinants_at_centroids = reshape(entity_determinants_at_centroids, 1, n_entity_elements) #scalar for each point and each element
# #             entity_jacobians_at_centroids = reshape(entity_jacobians_at_centroids, 3, 3, 1, n_entity_elements) #3x3 matrix for each point and each element

# #             global barycentres = cat(barycentres, entity_barycentres, dims=2)
# #             global global_centroids = cat(global_centroids, entity_global_centroids, dims=3)
# #             global jacobians_at_centroids = cat(jacobians_at_centroids, entity_jacobians_at_centroids, dims=4)
# #             global determinants_at_centroids = cat(determinants_at_centroids, entity_determinants_at_centroids, dims=2)

# #         end #each type in entity
# #     end # group entities
# # end #physical groups

# # n_comp_basis::Int64, basis_functions_at_centroids, n_orient::Int64 = gmsh.model.mesh.getBasisFunctions(tetrahedron_element_type, local_centroid[:], "HcurlLegendre0")
# # n_basis = (((length(basis_functions_at_centroids) ÷ n_comp_basis) ÷ 1) ÷ n_orient)
# # basis_functions_at_centroids = reshape(basis_functions_at_centroids, n_comp_basis, n_basis, 1, n_orient)


# # # now use basis functions to reconstruct solution at the centroid
# # reconstruction_points = Vector{Float64}(undef,0)
# # solution_at_points_3D = Vector{Complex{Float64}}(undef,0) #vector solution at array of points

# # for iele in range(1,edge_basis_3D.n_elements)
    
# #     basis_orientation = edge_basis_3D.element_basis_orientations[iele]
# #     basis_orientation += 1 #gmsh to julia offset
    
# #     for ipoint in 1:1 #only one centroid

# #         solution_at_point = zeros(Float64,3)
        
# #         reconstruction_point = global_centroids[:,ipoint, iele] 
                     
# #         global reconstruction_points = cat(reconstruction_points, reconstruction_point, dims=1)
# #         #solution_at_point = barycentres[:,iele]
# #         for iedge in eachindex(edge_basis_3D.element_edge_tags[:, iele])
            
# #             edge_tag = edge_basis_3D.element_edge_tags[iedge, iele]
            
# #             for ibasis_per_edge in 1:n_basis_per_edge
# #                 basis_index = (iedge-1)*n_basis_per_edge + ibasis_per_edge
# #                 @assert basis_index == iedge
# #                 # if (basis_index != mod(edge_tag, 3) + 1)
# #                 #     println("edge tag = ", edge_tag)
# #                 #     println("basis_index = ", basis_index)
# #                 #     @assert basis_index == mod(edge_tag-1, 3) + 1
# #                 # end
# #                 jacobian = jacobians_at_centroids[:, :, ipoint, iele]
# #                 basis_function = (jacobian')\basis_functions_at_centroids[:,basis_index,ipoint,basis_orientation]  

# #                 coefficient_index = (edge_tag - 1)*n_basis_per_edge + ibasis_per_edge;

# #                 @assert coefficient_index == edge_tag

# #                 solution_at_point += solution_3D[coefficient_index]*basis_function;
# #                 #solution_at_point = global_centroids[:,1,iele]
                

# #                 #println("adding basis contribution")
# #             end

# #         end
# #         if abs(imag(solution_at_point[1])) > 10 
# #             solution_at_point[1] = 0.0
# #         end
# #         if abs(imag(solution_at_point[2])) > 10 
# #             solution_at_point[2] = 0.0
# #         end
# #         if abs(imag(solution_at_point[3])) > 10 
# #             solution_at_point[3] = 0.0
# #         end


# #         global solution_at_points_3D = cat(solution_at_points_3D, solution_at_point, dims=1)
# #         #println("total solution at point done")
# #     end
# # end
# # solution_at_points_3D = reshape(solution_at_points_3D, 3, length(solution_at_points_3D[:]) ÷ 3)





# # # # # glyph_solution_2DTE_u = imag(glyph_solution_2DTE[1:3:end])
# # # # # glyph_solution_2DTE_v = imag(glyph_solution_2DTE[2:3:end])
# # # # # solution_plot = quiver(points_solution_2DTE[1:3:end], points_solution_2DTE[2:3:end], quiver2d=(0.001*glyph_solution_2DTE_u,0.001*glyph_solution_2DTE_v), label="2D TE Solution")
# # # # # display(solution_plot)

# # # # #let's try to plot some basis functions to see if they are correct

# # # # test_point_order = 6 

# # # # test_quad_rule = gmsh.model.mesh.getIntegrationPoints(triangle_element_type,"Gauss$test_point_order")
# # # # n_test_points = length(test_quad_rule[1]) ÷ 3
# # # # test_points = reshape(test_quad_rule[1], 3, n_test_points) #store as 3 x n_quad_points
# # # # test_weights = test_points[2] #not used

# # # # n_comp_basis::Int64, basis_functions_at_test_points, n_orient::Int64 = gmsh.model.mesh.getBasisFunctions(triangle_element_type, test_points[:], "HcurlLegendre0")
# # # # n_basis = (((length(basis_functions_at_test_points) ÷ n_comp_basis) ÷ n_test_points) ÷ n_orient)
# # # # basis_functions_at_test_points = reshape(basis_functions_at_test_points, n_comp_basis, n_basis, n_test_points, n_orient)



# # # # #get global test points
# # # # global_test_points = Array{Float64, 3}(undef, 3, n_test_points, 0)
# # # # jacobians_at_test_points = Array{Float64, 4}(undef, 3, 3, n_test_points, 0)
# # # # determinants_at_test_points = Array{Float64, 2}(undef, n_test_points, 0)  

# # # # for i in eachindex(mesh.physical_groups)

# # # #     group_dim = mesh.physical_groups[i][1]
# # # #     group_tag = mesh.physical_groups[i][2]
# # # #     group_entities = gmsh.model.getEntitiesForPhysicalGroup(group_dim, group_tag)

# # # #     #loop over entries in physical group
# # # #     for entity in group_entities
        
# # # #         entity_element_types = gmsh.model.mesh.getElementTypes(group_dim, entity)
        
# # # #         for itype in eachindex(entity_element_types)

# # # #             #TODO: expand this to a larger comparison as needed/desired
# # # #             if entity_element_types[itype] != triangle_element_type
# # # #                 continue #only interested in these elements for this basis
# # # #             end

# # # #             entity_jacobians_at_test_points, entity_determinants_at_test_points, entity_global_test_points = gmsh.model.mesh.getJacobians(triangle_element_type, test_points[:], entity)
# # # #             n_entity_elements = (length(entity_global_test_points) ÷ 3) ÷ n_test_points
# # # #             entity_global_test_points = reshape(entity_global_test_points, 3, n_test_points, n_entity_elements)
# # # #             entity_determinants_at_test_points = reshape(entity_determinants_at_test_points, n_test_points, n_entity_elements) #scalar for each point and each element
# # # #             entity_jacobians_at_test_points = reshape(entity_jacobians_at_test_points, 3, 3, n_test_points, n_entity_elements) #3x3 matrix for each point and each element

# # # #             global global_test_points = cat(global_test_points, entity_global_test_points, dims=3)
# # # #             global jacobians_at_test_points = cat(jacobians_at_test_points, entity_jacobians_at_test_points, dims=4)
# # # #             global determinants_at_test_points = cat(determinants_at_test_points, entity_determinants_at_test_points, dims=2)

# # # #         end #each type in entity
# # # #     end # group entities
# # # # end #physical groups

# # # # ele_index = 144
# # # # basis_index = 3
# # # # basis_function_at_all_points = Vector{Float64}(undef,0)
# # # # local_basis_functions = 
# # # # for iele in range(ele_index,ele_index)
    
# # # #     basis_orientation = edge_basis.element_basis_orientations[iele]
# # # #     basis_orientation += 1 #gmsh to julia offset
    
# # # #     for ipoint in 1:n_test_points

# # # #         basis_at_point = zeros(Float64,3)
# # # #         jacobian = jacobians_at_test_points[:, :, ipoint, iele]
# # # #         basis_function_at_point = jacobian'\basis_functions_at_test_points[:,basis_index,ipoint,basis_orientation]  
# # # #         global basis_function_at_all_points = cat(basis_function_at_all_points, basis_function_at_point, dims=1)        
# # # #     end
# # # # end
# # # # glyph_basis_u = basis_function_at_all_points[1:3:end]
# # # # glyph_basis_v = basis_function_at_all_points[2:3:end]
# # # # basis_plot = quiver(global_test_points[1,:,ele_index], global_test_points[2,:,ele_index], quiver2d=(0.0001*glyph_basis_u,0.0001*glyph_basis_v), label="Basis")
# # # # display(basis_plot)


# #=------------------------------------
#  Extract Solution and plot it in Gmsh 2D TM
# ------------------------------------=#

# t = gmsh.view.add("Test")
# gmsh.view.addHomogeneousModelData(t, 0, gmsh.model.getCurrent(), "NodeData", 1:length(gmsh_solution_2DTM), imag(gmsh_solution_2DTM))
# gmsh.view.write(t, "./solution.msh")
# gmsh.open("./solution.msh")
# gmsh.option.setNumber("Mesh.SurfaceEdges", 0)
# gmsh.fltk.run()
# gmsh.finalize()


# # # #=------------------------------------
# # #  Extract Solution and plot it in Gmsh 2D TM
# # # ------------------------------------=#

# # # #my unique nodes are not in the same order as gmsh's unique nodes !! we need to create a mapping to get that straightened out.
# # # TEx = gmsh.view.add("TEx")
# # # TEy = gmsh.view.add("TEy")
# # # TMz = gmsh.view.add("TMz")
# # # gmsh.view.addHomogeneousModelData(TEx, 0, gmsh.model.getCurrent(), "ElementData", mesh_2D.elements[3].tags, imag(solution_at_points_2DTE[1,:]))
# # # gmsh.view.addHomogeneousModelData(TEy, 0, gmsh.model.getCurrent(), "ElementData", mesh_2D.elements[3].tags, imag(solution_at_points_2DTE[2,:]))
# # # gmsh.view.write(TEy, "./solution.msh")
# # # gmsh.view.write(TEy, "./solution.msh")
# # # gmsh.open("./solution.msh")
# # # gmsh.option.setNumber("Mesh.SurfaceEdges", 0)
# # # gmsh.fltk.run()
# # # gmsh.finalize()




# # #my unique nodes are not in the same order as gmsh's unique nodes !! we need to create a mapping to get that straightened out.
# # Ex = gmsh.view.add("Ex")
# # Ey = gmsh.view.add("Ey")
# # Ez = gmsh.view.add("Ez")
# # gmsh.view.addHomogeneousModelData(Ex, 0, gmsh.model.getCurrent(), "ElementData", mesh_3D.elements[4].tags, imag(solution_at_points_3D[1,:]))
# # gmsh.view.addHomogeneousModelData(Ey, 0, gmsh.model.getCurrent(), "ElementData", mesh_3D.elements[4].tags, imag(solution_at_points_3D[2,:]))
# # gmsh.view.addHomogeneousModelData(Ez, 0, gmsh.model.getCurrent(), "ElementData", mesh_3D.elements[4].tags, imag(solution_at_points_3D[3,:]))
# # gmsh.view.write(Ex, "./solution.msh")
# # gmsh.view.write(Ey, "./solution.msh")
# # gmsh.view.write(Ez, "./solution.msh")
# # gmsh.open("./solution.msh")
# # gmsh.option.setNumber("Mesh.SurfaceEdges", 0)
# # gmsh.fltk.run()
# # gmsh.finalize()
