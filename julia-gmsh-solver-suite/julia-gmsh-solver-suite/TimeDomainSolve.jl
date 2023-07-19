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

# for gaussian source
include("Extras/Analytic.jl")

using Gmsh: gmsh
using SparseArrays
using StaticArrays
using LinearAlgebra

using CairoMakie
using ForwardDiff
using FFTW

function createWeirdChamber(; target_radius=0.4, tumour_radius=0.1)

    N = 11
    fibs = [1.0, 1.0]
    for i in 1:N-1
        push!(fibs, fibs[end] + fibs[end-1])
    end

    phi = 2 * pi * fibs ./ fibs[end]
    R0 = 2
    points = []
    for (i, angle) in enumerate(phi[begin:end-1])
        p = gmsh.model.occ.addPoint(R0 * cos(angle), R0 * sin(angle), 0)
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

function createTETest(; target_radius=0.4, tumour_radius=0.1)
    # create a cylindrical boundary
    circle_id = gmsh.model.occ.addCircle(0, 0, 0, 1)
    curve_loop_id = gmsh.model.occ.addCurveLoop([circle_id])
    surface_id = gmsh.model.occ.addPlaneSurface([curve_loop_id])

    # create fat
    fat_id = gmsh.model.occ.addCircle(0, 0, 0, target_radius)
    fat_loop_id = gmsh.model.occ.addCurveLoop([fat_id])
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

    gmsh.model.addPhysicalGroup(1, [curve_loop_id], 1000, "pec")

    gmsh.fltk.run()
end

function createMesh(geoOrder::Integer, lc_map::Dict{Int64,Float64})
    """
    Create mesh using lc per physical group
    """


    # SETUP MESH FIELDS based on lc_map
    fields = []

    lc_max = maximum([value for (key, value) in lc_map])

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

if gmsh.isInitialized() == 1
    gmsh.finalize()
end

# parameters for gaussian time derivative source
fc = 800e6 # centre frequency
ωc = fc * (2 * pi)
fs = 2e9
ts = 1 / fs
Δt = ts / 15
λ0 = c0 / fs

println("centre frequency = $(fc/1e6) [MHz]")
println("Freespace wavelength = $λ0 [m]")
println("Highest Perm wavelength = $(λ0/sqrt(12)) [m]")


# number of gaussian deviations
deviations = 6
limit = deviations / ωc
times = -limit:Δt:limit
println("Number of time steps needed: $(length(times))")
println("Nyquist frequency: $fs Δt: $Δt ")

magneticTimeFunction(t) = gaussDeriv(t, ωc)
currentTimeFunction(t) = μ0 * ForwardDiff.derivative(magneticTimeFunction, t)

λmesh = 0.25
# λmesh = 0.75

σair = 0.0
σfat = 0.02
σtumour = 0.05

ϵair = 1
ϵfat = 3
ϵtumour = 12

lc_freespace = λmesh / 5
lc_fat = lc_freespace / sqrt(3)
lc_tumour = lc_freespace / sqrt(12)

lc_map = Dict(
    2000 => lc_freespace,
    2001 => lc_fat,
    2002 => lc_tumour
)
ϵr_map = Dict(
    2000 => ϵair,
    2001 => ϵfat,
    2002 => ϵtumour
)

μr_map = Dict(
    2000 => 1.0,
    2001 => 1.0,
    2002 => 1.0
)

σ_map = Dict(
    2000 => σair,
    2001 => σfat,
    2002 => σtumour
)

gmsh.initialize()
createWeirdChamber()
createMesh(2, lc_map)
# gmsh.fltk.run()
# x::Int32 = 1.5

mesh = loadMesh()

basis_order = 2
triangle_space = HCurlElement(2, basis_order)
global_space = HCurlSpace(mesh, triangle_space)

N = numDOF(global_space)
pec_list = getPhysicsDOF(global_space, 1, 1000)
# pec_list = []
nonpec_list = setdiff(1:N, pec_list)

# integrate point source
evaluations, indexes = integratePointSource(global_space, -0.2, 0.2, 0.0, curl=true)

proj = SA_F64[0, 0, 1]
integrals = [dot(proj, eval) for eval in evaluations]
sourceIntegralsAllDOF = sparsevec(indexes, integrals, N)
source_integrals = Array(sourceIntegralsAllDOF[nonpec_list])

## setup FEM matrices

# for parameters
ele_ϵr = zeros(numElements(mesh.elements))
ele_σ = zeros(numElements(mesh.elements))
ele_μr_inv = zeros(numElements(mesh.elements))

for (key, value) in mesh.physics_map
    (dim, tag) = key
    if dim == 2
        for elementTag in mesh.physics_map[key]
            idx = searchsortedfirst(mesh.elements.tags, elementTag)
            ele_ϵr[idx] = ϵr_map[tag]
            ele_σ[idx] = σ_map[tag]
            ele_μr_inv[idx] = 1 / μr_map[tag]
        end
    end
end

T = massMatrix(global_space, coefficients=ele_ϵr)
R = massMatrix(global_space, coefficients=ele_σ)
# Rabc = integrateRobinBoundary(global_space, 1000, 1/η0)
# R += Rabc
S = curlMatrix(global_space, coefficients=ele_μr_inv)

# rescale time by c0
Β = 0.25
τ = Δt * c0
A = T + 0.5 * η0 * τ * R + Β * τ^2 * S
B = 2 * T - (1 - 2 * Β) * τ^2 * S
C = T - 0.5 * η0 * τ * R + Β * τ^2 * S

# remove pec
A = A[nonpec_list, nonpec_list]
B = B[nonpec_list, nonpec_list]
C = C[nonpec_list, nonpec_list]

F = cholesky(A)


Enext = zeros(size(A, 1))
Ecurr = zeros(size(A, 1))
Eprev = zeros(size(A, 1))
rhs = zeros(size(A, 1))

solutionCoeffs = Vector{Float64}[]
solutionTimes = Float64[]

# complexSolution = zeros(ComplexF64, size(T,1))

for (i, t) in enumerate(times)

    # compute contribution from sources

    # magnetic sources
    time_scalar = τ^2 * (Β * magneticTimeFunction(t + Δt) + (1 - 2 * Β) * magneticTimeFunction(t) + Β * magneticTimeFunction(t - Δt))

    # electric sources
    # time_scalar = τ^2 *(Β * f(t + Δt) + (1-2*Β) * f(t) + Β * f(t-Δt))

    global rhs = B * Ecurr - C * Eprev - source_integrals * time_scalar

    # compute Enext
    global Enext = F \ rhs

    # save solution
    if i % 4 == 1
        solution = zeros(size(T, 1))
        solution[nonpec_list] = Ecurr
        push!(solutionCoeffs, solution)
        push!(solutionTimes, t)

        # update complexSolution
        # complexSolution .+= (solution .* exp(1im * ωc * t))
    end

    # update Ecurr and Eprev
    global Eprev = Ecurr
    global Ecurr = Enext

    # output progress
    if i % 100 == 0
        println("Time $t step: $i/$(length(times))")
    end
end

# create complex solution directly
# k0 = ωc * sqrt(ϵ0*μ0)
# K = S - k0^2*T + 1im*ωc*μ0*R
# Kf = K[nonpec_list, nonpec_list]
# helmholtz_sol = Kf \ -source_integrals

# helmholtz_sol_full = copy(complexSolution)
# helmholtz_sol_full[nonpec_list] = helmholtz_sol

# helmholtz_sol_full ./= real(helmholtz_sol_full' * T * helmholtz_sol_full)
# complexSolution ./= real(complexSolution' * T * complexSolution)
# saveVTK("complex_field_fromtime", "E", global_space, [complexSolution, helmholtz_sol_full], levels=1)

# savePVD("time_series", global_space, solutionCoeffs, solutionTimes, levels=1)
saveVTK("time_series", "E", global_space, solutionCoeffs, levels=0, times=solutionTimes)
gmsh.finalize()