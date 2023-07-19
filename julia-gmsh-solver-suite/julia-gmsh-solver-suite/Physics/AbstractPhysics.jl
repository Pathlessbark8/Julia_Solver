#=-------------------------------------------------------------------
 Problem Type
-------------------------------------------------------------------=#

abstract type AbstractProblemType end

#=-------------------------------------------------------------------
 Field Types
-------------------------------------------------------------------=#

abstract type AbstractFieldType end

#=-------------------------------------------------------------------
 Constitutive Types
-------------------------------------------------------------------=#

abstract type AbstractConstitutives end
abstract type AbstractElementConstitutives <: AbstractConstitutives end
abstract type AbstractNodalConstitutives <: AbstractConstitutives end


#=-------------------------------------------------------------------
 Boundary Condition Types
-------------------------------------------------------------------=#

abstract type AbstractBoundaryCondition end
struct DirichletBoundaryCondition <: AbstractBoundaryCondition end  # boundary condition of the first kind
struct NeumannBoundaryCondition <: AbstractBoundaryCondition end    # boundary condition of the second kind
struct RobinBoundaryCondition <: AbstractBoundaryCondition end       # boundary condition of the third kind

