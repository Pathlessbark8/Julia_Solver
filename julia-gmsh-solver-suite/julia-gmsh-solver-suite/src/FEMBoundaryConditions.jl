using Parameters

#=------------------------------------
Boundary Condition Types Enum
-------------------------------------=#
@enum BoundaryConditionType begin
    PEC=1
    ABC=2
end

@with_kw struct BoundaryCondition
    name::String
    tag::Int64
    type::BoundaryConditionType
end
