using GLMakie

function plotHighOrder(
    space::SimpleFiniteElementSpace{HCurlTriangle},
    solution_coeffs::Vector{Tv}) where {Tv<:Union{Float64,ComplexF64}}
    
    triNodes = space.mesh.elements.elementNodes
    nodeCoords = space.mesh.nodes.coords
    faces = reshape(triNodes[1:3, :], (3, size(triNodes,2)))'
    verts = nodeCoords[1:2, :]'

    Interpolator, x, y = interpolate(space, space.local_space.quad_points)
    projection = Interpolator * solution_coeffs

    u = [real(p[1]) for p in projection]
    v = [real(p[2]) for p in projection]

    color = hypot.(u,v)
    
    f = Figure()
    ax = Axis(f[1, 1], aspect=1)
    poly!(ax, verts, faces, strokewidth=1, transparency=true, shading=true, color=RGBAf(0.0, 0.0, 0.0, 0.0))
    arrows!(ax, x, y, u, v, linecolor=color, arrowcolor=color, lengthscale=1e-1)
    display(f)
end

function MultiPlot(
    space::SimpleFiniteElementSpace{HCurlTriangle},
    solution_coeffs::Vector{Vector{Tv}}) where {Tv<:Union{Float64,ComplexF64}}

    N = length(solution_coeffs)
    # generic stuff for all plots
    # add triangles
    # first three nodes are primary nodes
    triNodes = space.mesh.elements.elementNodes
    nodeCoords = space.mesh.nodes.coords
    faces = reshape(triNodes[1:3, :], (3, size(triNodes,2)))'
    verts = nodeCoords[1:2, :]'

    x, y, u, v = interpolate(space)
    u = real.(u)
    v = real.(v)

    points = [Point2f(i,j) for (i,j) in zip(x,y)]
    directions = [Point2f(i,j) for (i,j) in zip(u,v)]
    colors = hypot.(u, v)

    f = Figure()
    ax = Axis(f[1, 1], aspect=1)
    
    slider = Slider(f[2,1], range=1:N, startvalue = 1, snap=true)
    
    directions = lift(slider.value) do idx
        _,_,u,v = interpolate(space, solution_coeffs[idx])
        u = real.(u)
        v = real.(v)
        [Point2f(i,j) for (i,j) in zip(u,v)]
    end

    poly!(ax, verts, faces, strokewidth=1, transparency=true, shading=true, color=RGBAf(0.0, 0.0, 0.0, 0.0))
    arrows!(ax, points, directions, arrowcolor = colors, lengthscale=1e-1)
    display(f)
end