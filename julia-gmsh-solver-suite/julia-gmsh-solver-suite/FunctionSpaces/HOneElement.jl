struct HOneElement{T} <: AbstractFiniteElement{T}
    order::Int32
    quad_weights::Vector{Float64}
    quad_points::Vector{Float64}
    functions::Array{Float64,4}
    gradient_functions::Array{Float64,4}
end

# type aliases
const HOneLine = HOneElement{Line}
const HOneTriangle = HOneElement{Triangle}
const HOneTetrahedron = HOneElement{Tetrahedron}

function HOneElement(dimension, order)
    @assert dimension >= 1 && dimension <= 3

    quadrature_order = 2 * (order + 1)

    # assumption: geo order doesn't matter for local function space
    # assume simplex mesh
    if dimension == 1
        type_number = gmsh.model.mesh.getElementType("line", 1)
        element_type = Line
    elseif dimension == 2
        type_number = gmsh.model.mesh.getElementType("triangle", 1)
        element_type = Triangle
    else
        type_number = gmsh.model.mesh.getElementType("tetrahedron", 1)
        element_type = Tetrahedron
    end

    # quad points for element type and gauss order
    quad_points, quad_weights = gmsh.model.mesh.getIntegrationPoints(type_number, "Gauss$quadrature_order")
    num_points = length(quad_weights)

    # the basis functions themselves
    num_components, local_functions, num_orientations =
        gmsh.model.mesh.getBasisFunctions(type_number, quad_points, "Lagrange$order")
    @assert num_components == 1
    num_functions = length(local_functions) ÷ num_components ÷ num_orientations ÷ num_points
    println("num_points = $num_points, order = $order, num_functions = $num_functions")
    local_functions = reshape(local_functions, (num_components, num_functions, num_points, num_orientations))
    println("size of local_functions is " * string(size(local_functions)))
    # the curl of the basis functions
    num_components, gradient_functions, _ =
        gmsh.model.mesh.getBasisFunctions(type_number, quad_points, "GradLagrange$order")
    @assert num_components == 3
    gradient_functions = reshape(gradient_functions, (num_components, num_functions, num_points, num_orientations))

    return HOneElement{element_type}(order, quad_weights, quad_points, local_functions, gradient_functions)
end

