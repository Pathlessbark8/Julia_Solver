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

# for point source
include("Extras/Analytic.jl")

using Gmsh: gmsh
using SparseArrays
using StaticArrays
using LinearAlgebra
using MAT

function createTETest(;target_radius = 0.4, tumour_radius = 0.1)
    # instead of circle, use flat faceted thing?
    # phi = LinRange(0, 2*pi, 14)[begin:end-1]
    # points = []
    # for (x, y) in zip(cos.(phi), sin.(phi))
    #     p = gmsh.model.occ.addPoint(x, y, 0)
    #     push!(points, p)
    # end

    # # from points, create lines
    # lines = []
    # for iline in eachindex(points)
    #     l = gmsh.model.occ.addLine(points[iline], points[1 + (iline) % length(points)])
    #     push!(lines, l)
    # end

    # curve_loop_id = gmsh.model.occ.addCurveLoop(lines)
    # surface_id = gmsh.model.occ.addPlaneSurface([curve_loop_id])

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
end

function createMesh(geoOrder::Integer, lc_map::Dict{Int64, Float64})
"""
Create mesh using lc per physical group
"""


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

if gmsh.isInitialized() == 1
    gmsh.finalize()
end

f0 = c0 / 1.0
f1 = c0 / 0.15

λmesh = 0.25

frequencies = LinRange(f0, f1, 100)

println("$(frequencies[begin]), $(frequencies[2] - frequencies[1]), $(frequencies[end])")

σfat = 0.02
σtumour = 0.05

ϵfat = 3 * ϵ0
ϵtumour = 12 * ϵ0

lc_freespace = λmesh / 5
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
createTETest()
createMesh(2, lc_map)

mesh = loadMesh()

basis_order = 2
triangleSpace = HCurlElement(2, basis_order)
globalSpace = HCurlSpace(mesh, triangleSpace)

# to hold ω^2μϵ
eps_vector = zeros(ComplexF64, numElements(mesh.elements))
# to hold jωμσ
sigma_vector = zeros(ComplexF64, numElements(mesh.elements))

N = numDOF(globalSpace)
pec_list = getPhysicsDOF(globalSpace, 1, 1000)
nonpec_list = setdiff(1:N, pec_list)

S = curlMatrix(globalSpace)

nSources = 10
Rsource = 0.8
phi = LinRange(0, 2*pi, nSources+1)[begin:end-1]
x = Rsource*cos.(phi)
y = Rsource*sin.(phi)

I_ms = Int32[]
J_ms = Int32[]
V_ms = Float64[]
for i = 1:nSources
    evaluation, idxs = integratePointSource(globalSpace, x[i], y[i], 0.0, curl=true)
    
    proj = SA_F64[0,0,1]
    integral = dot.(Ref(proj), evaluation)
    
    for j in eachindex(integral)
        push!(I_ms, i)
        push!(J_ms, idxs[j])
        push!(V_ms, integral[j])
    end
end
Ms = sparse(I_ms, J_ms, V_ms, nSources, N)
display(Ms)
B = Array(Ms')
Bf = B[nonpec_list, :]
U = zeros(ComplexF64, N, nSources)

solutions = []

for f in frequencies
    ω = 2 * π * f
    for (key, value) in mesh.physics_map
        (dim, tag) = key
        if dim == 2
            for elementTag in mesh.physics_map[key]
                idx = searchsortedfirst(mesh.elements.tags, elementTag)
                eps_vector[idx] = ω^2 * μ0 * eps_map[tag]
                sigma_vector[idx] = 1im * ω * μ0 * sigma_map[tag]
            end
        end
    end
    
    M = massMatrix(globalSpace, coefficients=eps_vector)
    R = massMatrix(globalSpace, coefficients=sigma_vector)

    K = S - M + R
    Kf = K[nonpec_list, nonpec_list]

    @time Uf = Kf \ Bf
    U[nonpec_list, :] = Uf

    measure = Ms * U
    push!(solutions, reshape(measure, (1, size(measure)...)))

end

matrix = reduce(vcat, solutions)

matwrite("data.mat", 
    Dict("data" => matrix, 
    "frequencies" => collect(frequencies)))
gmsh.finalize()