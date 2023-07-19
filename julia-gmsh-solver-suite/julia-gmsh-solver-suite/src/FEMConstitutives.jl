#=------------------------------------
Constants
-------------------------------------=#

global const ε0=8.854e-12
global const μ0=4*π*1e-7
global const c0=1/sqrt(ε0*μ0)
global const η0=sqrt(μ0/ε0)


#=------------------------------------
Constitutives Types Enum
-------------------------------------=#
@enum ConstitutivesType begin
    simple_real=1
    simple_complex=2
    debye_real=3
    debye_complex=4
end

#=------------------------------------
Simple Real Consititutives 
-------------------------------------=#
struct SimpleRealConstitutives
    name::String
    tag::Int64
    type::ConstitutivesType
    ε::Float64
    μ::Float64
    σ::Float64
end

#=------------------------------------
Simple Complex Consititutives 
-------------------------------------=#
struct SimpleComplexConstitutives
    name::String
    tag::Int64
    type::ConstitutivesType
    ε::Complex{Float64}
    μ::Complex{Float64}
end

#=------------------------------------
Constitutives Struct
-------------------------------------=#
struct ComplexConstitutives
    element_ε_r::Vector{Complex{Float64}}
    element_μ_r::Vector{Complex{Float64}}
end

# #=------------------------------------
# setupElementConstitutives()
# -- assign constitutives based on tags and available physics_tags
# -------------------------------------=#
# function setupElementConstitutives(mesh::Mesh, basis::NodalBasis, constitutives_table::Dict{Int64,Any})

#     n_elements = mesh.elements[3].count #assuming triangles here
#     element_complex_relative_permittivity = Vector{Complex{Float64}}(undef,n_elements)
#     element_complex_relative_permeability = Vector{Complex{Float64}}(undef,n_elements)

#     for iele in 1:n_elements
#         ele_physics_tag = mesh.elements[3].physics_tags[iele]
#         constitutives = constitutives_table[ele_physics_tag] #will throw a Key Error if key not found

#         if constitutives.type != simple_complex::ConstitutivesType
#             @error("Only simple complex constitutives are currently supported")
#         end
#         element_complex_relative_permittivity[iele] = constitutives.ε
#         element_complex_relative_permeability[iele] = constitutives.μ

#     end
#     element_constitutives = ComplexConstitutives(element_complex_relative_permittivity, element_complex_relative_permeability)
#     return element_constitutives
# end


# function setupElementConstitutives3D(mesh::Mesh, basis::EdgeBasis, constitutives_table::Dict{Int64,Any})

#     n_elements = mesh.elements[4].count #assuming triangles here
#     element_complex_relative_permittivity = Vector{Complex{Float64}}(undef,n_elements)
#     element_complex_relative_permeability = Vector{Complex{Float64}}(undef,n_elements)

#     for iele in 1:n_elements
#         ele_physics_tag = mesh.elements[4].physics_tags[iele]
#         constitutives = constitutives_table[ele_physics_tag] #will throw a Key Error if key not found

#         if constitutives.type != simple_complex::ConstitutivesType
#             @error("Only simple complex constitutives are currently supported")
#         end
#         element_complex_relative_permittivity[iele] = constitutives.ε
#         element_complex_relative_permeability[iele] = constitutives.μ

#     end
#     element_constitutives = ComplexConstitutives(element_complex_relative_permittivity, element_complex_relative_permeability)
#     return element_constitutives
# end

