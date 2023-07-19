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
using SuiteSparse
SuiteSparse.UMFPACK.umf_ctrl[8] = 0 # disable iterative refinement

using KrylovKit

const M_TO_MM::Float64 = 1e-3
const SHIFT::Float64 = 1e-2
const CHAMBER_HEIGHT::Float64 = 0.14
const CHAMBER_RADIUS::Float64 = 0.115
const FLANGE_RADIUS::Float64 = 0.14

# phantom shape
const PHANTOM_LENGTH::Float64 = 96.4 / 1000
const PHANTOM_RADIUS::Float64 = 50 / 1000

# antenna
const ANTENNTA_RADIUS::Float64 = 0.9*CHAMBER_RADIUS


function addCircPlane(x,y,z,r)
    line = gmsh.model.occ.addCircle(x,y,z,r)
    loop = gmsh.model.occ.addCurveLoop([line])
    surf = gmsh.model.occ.addPlaneSurface([loop])
    return surf
end

function createdSmoothChamber()
    # v = gmsh.model.occ.importShapes("Chamber_v5_LC_FlatInner.iges")
    # # convert to m from mm
    # gmsh.model.occ.dilate(v, 0, 0, 0, M_TO_MM, M_TO_MM, M_TO_MM)
    # gmsh.model.occ.translate(v, 0, 0, SHIFT)

    # quadratic spline
    # x = LinRange(0, CHAMBER_RADIUS, 10)
    # points = []
    # for ix in x
    #     a = CHAMBER_HEIGHT / CHAMBER_RADIUS^2
    #     z = a*ix^2 - CHAMBER_HEIGHT
    #     p = gmsh.model.occ.addPoint(ix, 0.0, z)
    #     push!(points, p)
    # end
    # spline = gmsh.model.occ.addSpline(points)

    p1 = gmsh.model.occ.addPoint(0, 0, -CHAMBER_HEIGHT)
    p2 = gmsh.model.occ.addPoint(CHAMBER_RADIUS, 0, -CHAMBER_HEIGHT)
    p3 = gmsh.model.occ.addPoint(CHAMBER_RADIUS, 0, 0)
    points = [p1, p2, p3]

    spline = gmsh.model.occ.addBSpline(points)

    revolve = gmsh.model.occ.revolve([(1,spline)], 0,0,-CHAMBER_HEIGHT, 0,0,1,2*pi)
    
    topLoop = gmsh.model.occ.addCurveLoop([3])
    topSurf = gmsh.model.occ.addPlaneSurface([topLoop])
    
    # create a volume
    shell = gmsh.model.occ.addSurfaceLoop([1,2])
    bottomVol = gmsh.model.occ.addVolume([shell])

    # remove extra surface
    gmsh.model.occ.remove([(2,1)])

    # top bubble
    eps = 0.0
    sphere = gmsh.model.occ.addSphere(0,0,-eps,FLANGE_RADIUS)
    cutTool = gmsh.model.occ.addCylinder(0,0,-eps,0,0,-FLANGE_RADIUS, FLANGE_RADIUS)
    halfSphere, _ = gmsh.model.occ.cut([(3,sphere)], [(3,cutTool)])

    # intersect
    new, old = gmsh.model.occ.fuse(halfSphere, [(3,bottomVol)])
    
    # we now have one volume

    # add breast phantom
    extraHeight = max(0, 2*PHANTOM_RADIUS - PHANTOM_LENGTH)
    longCyl = gmsh.model.occ.addCylinder(0,0,extraHeight,0,0, -(PHANTOM_LENGTH-PHANTOM_RADIUS) - extraHeight, PHANTOM_RADIUS)
    phantomSphere = gmsh.model.occ.addSphere(0,0,-(PHANTOM_LENGTH - PHANTOM_RADIUS), PHANTOM_RADIUS)
    phantom, _ = gmsh.model.occ.fuse([(3,longCyl)], [(3,phantomSphere)])

    if extraHeight > 0
        extraCyl = gmsh.model.occ.addCylinder(0,0,0,0,0,extraHeight,PHANTOM_RADIUS)
        phantom, _ = gmsh.model.occ.cut(phantom, [(3,extraCyl)])
    end
    
    # cut out breast phantom from vol
    newVol, _ = gmsh.model.occ.cut([(3,1)], [(3,2)], -1, true, false)
    
    # add tumour in phantom
    tum = gmsh.model.occ.addSphere(0.02, 0.01,-0.04, 0.015)
    phantom, _ = gmsh.model.occ.cut(phantom, [(3,tum)], -1, true, false)

    # add top pec plate
    outerLoop = gmsh.model.occ.addCurveLoop([17])
    innerLoop = gmsh.model.occ.addCurveLoop([7, 8])
    topPec = gmsh.model.occ.addPlaneSurface([outerLoop, innerLoop])
    
    newVol, map = gmsh.model.occ.fragment(newVol, [(2, topPec)])
    # gmsh.model.occ.removeAllDuplicates()

    
    gmsh.model.occ.synchronize()
    # add physical groups
    # volumetric
    gmsh.model.addPhysicalGroup(3, [4,5], 3000, "air")
    gmsh.model.addPhysicalGroup(3, [2], 3001, "phantom")
    gmsh.model.addPhysicalGroup(3, [3], 3002, "sphere")

    # surfaces
    gmsh.model.addPhysicalGroup(2, [11,13,14], 2000, "pec")
    gmsh.model.addPhysicalGroup(2, [12], 2001, "abc")
end

"""
Create mesh using lc per physical group
"""
function create3DMesh(
    geoOrder::Integer, 
    lc_map::Dict{Int64, Float64})

    # SETUP MESH FIELDS based on lc_map
    fields = []

    lc_max = maximum([value for (key,value) in lc_map])

    physical_surfaces = gmsh.model.getPhysicalGroups(3)

    for (dim, tag) in physical_surfaces
        if haskey(lc_map, tag)
            entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)

            lc = lc_map[tag]

            field = gmsh.model.mesh.field.add("Constant")
            gmsh.model.mesh.field.setNumbers(field, "VolumesList", entities)
            gmsh.model.mesh.field.setNumber(field, "VIn", lc)
            gmsh.model.mesh.field.setNumber(field, "VOut", lc_max)

            push!(fields, field)
        end
    end
    
    min_field = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(min_field)
    
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(geoOrder)
end

"""
Create points sampling -z half of unit sphere

https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
"""
function fibonacciHalfSphere(samples)
    points = Float64[]
    phi = pi * (sqrt(5) - 1)

    for i in 0:samples-1
        z = - *(i / (samples - 1))
        radius = sqrt(1-z^2)

        theta = phi*i

        x = cos(theta)*radius
        y = sin(theta)*radius

        append!(points, [x,y,z])
    end
    return reshape(points, (3, samples))
end

function createSources(n_sources, antenna_radius, chamber_radius)
    # potential antenna positions
    points = fibonacciHalfSphere(n_sources)
    points .*= antenna_radius
    points[3,:] .-= 0.1*chamber_radius

    
    dir = cross.(Ref(SA_F64[0,0,1]), reinterpret(SVector{3,Float64}, points))
    lengths = norm.(dir)
    dir = reinterpret(Float64, dir)

    mask = vec(lengths .≈ 0)
    dir[:,mask] .= [0,1,0]
    lengths[mask] .= 1.0
    dir ./= lengths

    return points, dir
end

if gmsh.isInitialized() == 1
    gmsh.finalize()
end

# physics

# going to fun at 1 GHz
f0 = 2.4e9
ω0 = 2*pi*f0
λ0 = c0 / f0
# 5 samples per wavelength (high order basis functions)
basis_order = 2
λ_mesh = λ0 / 5

lc_map = Dict(
    3000 => λ_mesh,
    3001 => λ_mesh / sqrt(3),
    3002 => λ_mesh / sqrt(12)
)

ϵr_map = Dict(
    3000 => 1.0, # air
    3001 => 3.0, # fat
    3002 => 12.0 # tumour
)

σ_map = Dict(
    3000 => 0.0, # air
    3001 => 0.1, # fat 
    3002 => 0.2 # tumour
)

# just so I dont forget
bc_map = Dict(
    2000 => PEC(E()),
    2001 => ABC()
)

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

createdSmoothChamber()
create3DMesh(2, lc_map)

mesh = loadMesh()


# setup finite element space
local_element = HCurlElement(3, basis_order)
global_space = HCurlSpace(mesh, local_element)

N = numDOF(global_space)
# create pec/dirichlet lists
pec_list = getPhysicsDOF(global_space, 2, 2000)
non_pec_list =  setdiff(1:N, pec_list)

epsr_vector = zeros(numElements(mesh.elements))
sigma_vector = zeros(numElements(mesh.elements))

eps_inside = zeros(Float64, numElements(mesh.elements))
eps_outside = zeros(Float64, numElements(mesh.elements))

for (key, value) in mesh.physics_map
    (dim, tag) = key
    if dim == 3
        for elementTag in mesh.physics_map[key]
            idx = searchsortedfirst(mesh.elements.tags, elementTag)
            epsr_vector[idx] = ϵr_map[tag]
            sigma_vector[idx] = σ_map[tag]

            if tag != 3000
                eps_inside[idx] = ϵr_map[tag]
            else
                eps_outside[idx] = ϵr_map[tag]
            end
        end
    end
end

# create matrices
Minside = massMatrix(global_space, coefficients=eps_inside)
dropzeros!(Minside)
Moutside = massMatrix(global_space, coefficients=eps_outside)
dropzeros!(Moutside)

M = massMatrix(global_space, coefficients=epsr_vector)
R = massMatrix(global_space, coefficients=sigma_vector)
S = curlMatrix(global_space)
Kabc = integrateRobinBoundary(global_space, 2001, 1.0)

M = M*(ω0^2 * ϵ0 * μ0)
R = R*(1im*ω0*μ0)
R += 1im*ω0*sqrt(ϵ0*μ0)*Kabc

K = S - M + R

Kf = K[non_pec_list, non_pec_list]
Mf = M[non_pec_list, non_pec_list]
Rf = R[non_pec_list, non_pec_list]
NP = size(Kf, 1)

A = [Kf spzeros(NP, NP); spzeros(NP, NP) sparse(I, NP, NP)]
B = [2*Mf - Rf Mf; sparse(I, NP, NP) spzeros(NP, NP)]

F = lu(A)

function bigmatvec(x)
    b = B*x
    return F \ b
end

λ, v, info = eigsolve(
    bigmatvec, 2*NP, 150, 
    EigSorter(f -> abs(f)*(real(1/f) <= -0.99)), 
    ishermitian=false, krylovdim=200)
# because we used shift-invert
λ = 1 ./ λ
ω = ω0 * (1 .+ λ)

modes = [zeros(ComplexF64, numDOF(global_space)) for i in 1:length(v)]
for (i, mode) in enumerate(v)
    modes[i][non_pec_list] = mode[1:NP]
end

energy_ratio = [real(mode' * Moutside * mode) / real(mode' * Minside * mode) for mode in modes]
order = sortperm(imag.(λ))
λ = λ[order]
ω = ω[order]
modes = modes[order]
display(ω)



# save to paraview
saveVTK(
    "3D_chamber",
    "eigs",
    global_space,
    modes,
    levels=2
)

# no longer need gmsh
gmsh.finalize()


    
