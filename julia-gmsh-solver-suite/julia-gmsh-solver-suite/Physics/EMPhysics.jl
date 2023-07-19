
#=-------------------------------------------------------------------
 Fundamental Constants
-------------------------------------------------------------------=#
const c0 = 299_792_458 # precise constant
const μ0 = 4*π*1e-7 # precise constant
const ε0 = 1 / (μ0 * c0^2)
const Z0 = sqrt(μ0/ε0)
const Y0 = 1/Z0

#=-------------------------------------------------------------------
 Problem Types
-------------------------------------------------------------------=#

struct EMTMProblemType <: AbstractProblemType end
struct EMTEProblemType <: AbstractProblemType end
struct EM3DProblemType <: AbstractProblemType end

#=-------------------------------------------------------------------
 Field Types
-------------------------------------------------------------------=#

struct ElectricField <: AbstractFieldType end
struct MagneticField <: AbstractFieldType end

#=-------------------------------------------------------------------
 Constitutives Enums and Structs
-------------------------------------------------------------------=#

abstract type AbstractEMConstitutives <: AbstractConstitutives end
abstract type AbstractEMElementConstitutives <: AbstractElementConstitutives end
abstract type AbstractEMNodalConstitutives <: AbstractNodalConstitutives end

#=------------------------------------
Constitutives Types Enum
-------------------------------------=#
@enum EMConstitutivesEnum begin
    simple_real=1
    simple_complex=2
    debye_real=3
    debye_complex=4
end
export EMConstitutivesEnum

#=------------------------------------
Simple Real Consititutives 
-------------------------------------=#
struct SimpleRealEMConstitutives <: AbstractEMConstitutives
    name::String
    tag::Int64
    ε::Float64
    μ::Float64
    σ::Float64
end

#=------------------------------------
Simple Complex Consititutives 
-------------------------------------=#
struct SimpleComplexEMConstitutives <: AbstractEMConstitutives
    name::String
    tag::Int64
    ε::Complex{Float64}
    μ::Complex{Float64}
end

#=------------------------------------
Constitutives Struct
-------------------------------------=#
struct ComplexEMElementConstitutives <: AbstractEMElementConstitutives
    element_ε_r::Vector{Complex{Float64}}
    element_μ_r::Vector{Complex{Float64}}
end

#=-------------------------------------------------------------------
 Boundary Condition Enums and Structs
-------------------------------------------------------------------=#

abstract type AbstractEMBoundaryCondition <: AbstractBoundaryCondition end
abstract type AbstractEMElementBoundaryCondition end #not sure if this should inherit

#=------------------------------------
Boundary Condition Type Structs
-------------------------------------=#

struct PECEMBoundaryCondition <: AbstractEMBoundaryCondition end
struct PMCEMBoundaryCondition <: AbstractEMBoundaryCondition end
struct ABCEMBoundaryCondition <: AbstractEMBoundaryCondition end
struct IMPEMBoundaryCondition{T<:Union{Float64, ComplexF64}} <: AbstractEMBoundaryCondition
    impedance::T    
end
EMABCBoundaryCondition(T) = EMIMPBoundaryCondition(T)

struct EMElementBoundaryCondition <: AbstractEMElementBoundaryCondition

    ele_tag::Int64
    ele_index::Int64
    ele_face_index::Int64
    
    boundary_type::AbstractEMBoundaryCondition
    boundary_physics_tag::Int64 
    boundary_element_tag::Int64 # lower dimension
    
end


# # electromagnetic specifc boundary condition aliases
# ABC() = Robin() 
# PEC(E::Electric) = Dirichlet()
# PEC(H::Magnetic) = Neumann()

# PMC(E::Electric) = Neumann()
# PMC(H::Magnetic) = Dirichlet()
# # Higher than first order ABC doesn't necessarily map to Robin boundary conditions
# struct ABC2 <: AbstractPhysics end
# export PEC, PMC, ABC, ABC2



