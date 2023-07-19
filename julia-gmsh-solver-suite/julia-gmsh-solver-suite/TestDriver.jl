using Revise
includet("Mesh/Mesh.jl")
using .Meshes
Revise.track(Meshes, "Mesh/GmshMeshUtils.jl")
Revise.track(Meshes, "Mesh/SimpleMesh.jl")

includet("Physics/Physics.jl")
using .Physics

includet("FunctionSpaces/FunctionSpaces.jl")
using .FunctionSpaces
Revise.track(FunctionSpaces, "FunctionSpaces/HCurlElement.jl")
Revise.track(FunctionSpaces, "FunctionSpaces/HCurlSpace.jl")

includet("Extras/Visualization.jl")
using .Visualization
Revise.track(Visualization, "Extras/MakieViz.jl")

using Gmsh: gmsh
using SparseArrays
using StaticArrays
using LinearAlgebra

include("Extras/Analytic.jl")

using GLMakie

# eigensolvers
using KrylovKit
function createWeirdChamber(;target_radius = 0.4, tumour_radius = 0.1)

    N = 11
    fibs = [1.0,1.0]
    for i in 1:N-1
        push!(fibs, fibs[end]+fibs[end-1])
    end
    
    phi = 2*pi*fibs ./ fibs[end]
    R0 = 2
    points = []
    for (i, angle) in enumerate(phi[begin:end-1])
        p = gmsh.model.occ.addPoint(R0*cos(angle), R0*sin(angle), 0)
        push!(points, p)
    end

    push!(points, points[begin])

    spline = gmsh.model.occ.addBSpline(points)
    # # create a cylindrical boundary
    # circle_id = gmsh.model.occ.addCircle(0, 0, 0, 1)
    curve_loop_id = gmsh.model.occ.addCurveLoop([spline])
    surface_id = gmsh.model.occ.addPlaneSurface([curve_loop_id])

    # create fat
    fat_id = gmsh.model.occ.addCircle(0, 0, 0, target_radius)
    fat_loop_id = gmsh.model.occ.addCurveLoop([fat_id])
    fat_surface_id = gmsh.model.occ.addPlaneSurface([fat_loop_id])

    # create tumour
    tumour_id = gmsh.model.occ.addCircle(-0.2, -0.1, 0, tumour_radius)
    tumour_loop_id = gmsh.model.occ.addCurveLoop([tumour_id])
    tumour_surface_id = gmsh.model.occ.addPlaneSurface([tumour_loop_id])

    gmsh.model.occ.synchronize()

    # cut out circles from each other

    outTags, outMap = gmsh.model.occ.fragment([(2, surface_id)], [(2, fat_surface_id), (2, tumour_surface_id)])
    gmsh.model.occ.synchronize()

    gmsh.model.addPhysicalGroup(2, [4], 2000, "freespace")
    gmsh.model.addPhysicalGroup(2, [5], 2001, "fat")
    gmsh.model.addPhysicalGroup(2, [3], 2002, "tumour")

    gmsh.model.addPhysicalGroup(1, [curve_loop_id], 1000, "pec")

    
end
function createTETest(;target_radius = 0.6, tumour_radius = 0.1)
    # create a cylindrical boundary
    circle_id = gmsh.model.occ.addCircle(0, 0, 0, 1)
    curve_loop_id = gmsh.model.occ.addCurveLoop([circle_id])
    surface_id = gmsh.model.occ.addPlaneSurface([curve_loop_id])

    # create fat
    # create spiky interface
    lines = []
    phi = LinRange(0, 2*pi,12)[begin:end-1]
    dphi = phi[2]-phi[1]
    h = 0.3
    r1 = target_radius + h/2
    r2 = target_radius - h/2
    for i=1:11
        p1 = gmsh.model.occ.addPoint(r1*cos(phi[i]), r1*sin(phi[i]), 0)
        p2 = gmsh.model.occ.addPoint(r2*cos(phi[i] + dphi/2), r2*sin(phi[i] + dphi/2), 0)
        p3 = gmsh.model.occ.addPoint(r1*cos(phi[i]+dphi), r1*sin(phi[i]+dphi), 0)

        l1 = gmsh.model.occ.addLine(p1, p2)
        l2 = gmsh.model.occ.addLine(p2, p3)
        push!(lines, l1)
        push!(lines, l2)
    end
    gmsh.model.occ.removeAllDuplicates()
    fat_loop_id = gmsh.model.occ.addCurveLoop(lines)
    # fat_id = gmsh.model.occ.addCircle(0, 0, 0, target_radius)
    # fat_loop_id = gmsh.model.occ.addCurveLoop([fat_id])
    fat_surface_id = gmsh.model.occ.addPlaneSurface([fat_loop_id])

    # create tumour
    tumour_id = gmsh.model.occ.addCircle(-0.2, -0.1, 0, tumour_radius)
    tumour_loop_id = gmsh.model.occ.addCurveLoop([tumour_id])
    tumour_surface_id = gmsh.model.occ.addPlaneSurface([tumour_loop_id])

    # cut out circles from each other

    outTags, outMap = gmsh.model.occ.fragment([(2, surface_id)], [(2, fat_surface_id), (2, tumour_surface_id)])
    gmsh.model.occ.synchronize()

    gmsh.model.addPhysicalGroup(2, [4], 2000, "freespace")
    gmsh.model.addPhysicalGroup(2, [5], 2001, "fat")
    gmsh.model.addPhysicalGroup(2, [3], 2002, "tumour")

    gmsh.model.addPhysicalGroup(1, [24], 1000, "pec")

    # gmsh.fltk.run()
    # x::Int32 = 1.5
end

function createCircle()
    circle_id = gmsh.model.occ.addCircle(0, 0, 0, 1)
    curve_loop_id = gmsh.model.occ.addCurveLoop([circle_id])
    surface_id = gmsh.model.occ.addPlaneSurface([curve_loop_id])

    p1 = gmsh.model.occ.addPoint(-0.25, -0.25, 0)
    p2 = gmsh.model.occ.addPoint(-0.25, 0.25, 0)
    p3 = gmsh.model.occ.addPoint(0.25, 0.25, 0)
    p4 = gmsh.model.occ.addPoint(0.25, -0.25, 0)

    l1 = gmsh.model.occ.addLine(p1, p2)
    l2 = gmsh.model.occ.addLine(p2, p3)
    l3 = gmsh.model.occ.addLine(p3, p4)

    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    gmsh.model.mesh.embed(1, [l1, l2, l3], 2, surface_id)

    gmsh.model.addPhysicalGroup(1, [l1, l2, l3], 2000, "pec")

    gmsh.model.addPhysicalGroup(2, [surface_id], 3000, "freespace")
    circle_boundary = gmsh.model.getBoundary((2, circle_id), true, false, false)
    circle_boundary = getfield.(circle_boundary, 2)
    gmsh.model.addPhysicalGroup(1, circle_boundary, 2001, "abc")

    # gmsh.fltk.run()
end

function createMesh(geoOrder::Integer, lc::Float64)
    points = gmsh.model.getEntities(0)
    gmsh.model.mesh.setSize(points, lc)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(geoOrder)
end

"""
Create mesh using lc per physical group

"""
function createMesh(geoOrder::Integer, lc_map::Dict{Int64, Float64})
    # SETUP MESH FIELDS based on lc_map
    fields = []

    lc_max = maximum([value for (key,value) in lc_map])

    physical_surfaces = gmsh.model.getPhysicalGroups(2)
    for (dim, tag) in physical_surfaces
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)

        lc = lc_map[tag]

        field = gmsh.model.mesh.field.add("Constant")
        gmsh.model.mesh.field.setNumbers(field, "SurfacesList", entities)
        gmsh.model.mesh.field.setNumber(field, "VIn", lc)
        gmsh.model.mesh.field.setNumber(field, "VOut", lc_max)

        push!(fields, field)
    end

    min_field = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(min_field)
    
    tags = gmsh.model.mesh.field.list()
    
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(geoOrder)
end

# reset gmsh
if gmsh.isInitialized() == 1
    gmsh.finalize()
end


λ0 = 0.75
f = c0 / λ0

println("λ0=$λ0 f=$f")

ω0 = 2*pi*f
k2 = ω0^2 * ϵ0 * μ0
k0 = sqrt(k2)

σfat = 0.02
σtumour = 0.05

ϵfat = 3 * ϵ0
ϵtumour = 12 * ϵ0

lc_freespace = λ0 / 5
lc_fat = lc_freespace / sqrt(3)
lc_tumour = lc_freespace / sqrt(12)

lc_map = Dict(
    2000 => lc_freespace,
    2001 => lc_fat,
    2002 => lc_tumour
)

eps_map = Dict(
    2000 => ϵ0,
    2001 => ϵfat,
    2002 => ϵtumour
)

sigma_map = Dict(
    2000 => 0.0,
    2001 => σfat,
    2002 => σtumour
)

gmsh.initialize()

# createCircle()
# createTETest()
createWeirdChamber()
createMesh(2, lc_map)
mesh = loadMesh()

basis_order = 2
triangleSpace = HCurlElement(2, basis_order)
globalSpace = HCurlSpace(mesh, triangleSpace)

N = numDOF(globalSpace)

eps_vector = zeros(ComplexF64, numElements(mesh.elements))
sigma_vector = zeros(ComplexF64, numElements(mesh.elements))

eps_inside = zeros(Float64, numElements(mesh.elements))
eps_outside = zeros(Float64, numElements(mesh.elements))

for (key, value) in mesh.physics_map
    (dim, tag) = key
    if dim == 2
        for elementTag in mesh.physics_map[key]
            idx = searchsortedfirst(mesh.elements.tags, elementTag)
            eps_vector[idx] = ω0^2 * μ0 * eps_map[tag]
            sigma_vector[idx] = 1im*ω0*μ0*sigma_map[tag]
            
            if tag != 2000
                eps_inside[idx] = eps_map[tag]
            else
                eps_outside[idx] = eps_map[tag]
            end
        end
    end
end

Minside = massMatrix(globalSpace, coefficients=eps_inside)
dropzeros!(Minside)
Moutside = massMatrix(globalSpace, coefficients=eps_outside)
dropzeros!(Moutside)

M = massMatrix(globalSpace, coefficients=eps_vector)
R = massMatrix(globalSpace, coefficients=sigma_vector)
S = curlMatrix(globalSpace)

pec_list = getPhysicsDOF(globalSpace, 1, 1000)

# pec_list = []
# Y0 = 1im * k0
# Kabc = integrateRobinBoundary(globalSpace, 1000, Y0)
# R += Kabc

K = S - M + R


nonpec_list = setdiff(1:N, pec_list)

Kf = K[nonpec_list, nonpec_list]
Mf = M[nonpec_list, nonpec_list]
Rf = R[nonpec_list, nonpec_list]
NP = size(Kf, 1)

A = [Kf spzeros(NP, NP); spzeros(NP, NP) sparse(I, NP, NP)]
B = [2*Mf - Rf Mf; sparse(I, NP, NP) spzeros(NP, NP)]

F = lu(A)

function bigmatvec(x)
    b = B*x
    return F \ b
end

λ, v, info = eigsolve(
    bigmatvec, 2*NP, 100, 
    EigSorter(f -> abs(f)*(real(1/f) <= -0.99)), 
    ishermitian=false, krylovdim=150)
println(info.converged)
# because we used shift-invert
λ = 1 ./ λ

E = [zeros(ComplexF64, numDOF(globalSpace)) for i in 1:length(v)]

for (i, mode) in enumerate(v)
    E[i][nonpec_list] = mode[1:NP]
end

energy_ratio = [real(mode' * Moutside * mode) / real(mode' * Minside * mode) for mode in E]
order = sortperm(energy_ratio)

λ = λ[order]
display(λ)
ω = ω0 * (1 .+ λ)
display(ω)

E = E[order]


fig = Figure()
# ax1 = Axis(fig[1,1], xlabel="idx", ylabel="Re{λ} [rad/s]", aspect=1)
# ax2 = Axis(fig[1,2], xlabel="idx", ylabel="Im{λ} [rad/s]", aspect=1)
ax1 = Axis(fig[1,1], xlabel="Re{λ} [rad/s]", ylabel="Im{λ} [rad/s]", aspect=1)
# scatter!(ax1, real.(ω))
# scatter!(ax2, imag.(ω))
scatter!(ax1, real.(ω), imag.(ω))
ylims!(0, 1e9)
display(fig)
save("eigenvalues_pec.png", fig)

saveVTK("eigs_pec", "E", globalSpace, E, levels=2)

# test case for point source
# remove PEC
# Kf = K[nonpec_list, nonpec_list]
# # create a test source term
# integral, tag = addPointSource(0,-0.75,0,basis_order)
# idx = searchsortedfirst(globalSpace.mesh.elements.tags, tag)
# @assert globalSpace.mesh.elements.tags[idx] == tag
# dof_idxs = globalSpace.dof_map[:,idx]

# b = zeros(eltype(Kf), N)
# E = zeros(eltype(Kf), N)
# b[dof_idxs] = integral

# bf = b[nonpec_list]

# F = lu(Kf)
# Ef = F \ bf

# E[nonpec_list] = Ef

# plotHighOrder(globalSpace, E)


gmsh.finalize()

