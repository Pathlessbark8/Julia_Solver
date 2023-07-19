using Revise

if isdefined(@__MODULE__, :EIL)
    
else
    include("EIL.jl")
end

using LinearAlgebra
using Plots
using Parameters
using NodesAndModes
using StaticArrays
using .EIL

using Gmsh:gmsh

gmsh.finalize()
setupGmsh() 

# ==================================================================================
# The purpose of this test script is to figure out:
# 1) how to compute the surface normal of a triangular element from the jacobians
# 
# The answer is:
# If we let the jacobian of the triangle be [a1 a2 a3] as three columns, then 
# a3 = cross(a1, a2) and is equal to the unit surface normal to the triangle.
#
# Note that this vector is also equal to the third column of transpose(inv(J)) (or 
# the third row of inv(J)) as the third column of J is normalized.
# ==================================================================================

#basic types for this unit test
tri_type = gmsh.model.mesh.getElementType("triangle", 1)

#-----------------------------------------------------------------------------------
# To answer question #1 we will start by creating a circular mesh with a normal
# zaxis = [alpha, beta, gamma]
#-----------------------------------------------------------------------------------

zaxis = [2, 0.5, 1]

circle_id = gmsh.model.occ.addCircle(0,0,0, 1, -1, 0, 2*pi, zaxis);
lc1 = 0.3
gmsh.model.occ.mesh.setSize(gmsh.model.occ.getEntities(0), lc1)
gmsh.model.occ.synchronize()
curve_loop_id = gmsh.model.geo.addCurveLoop([circle_id])
surface_id = gmsh.model.geo.addPlaneSurface([curve_loop_id]);

gmsh.model.geo.synchronize();
gmsh.model.addPhysicalGroup(2, [surface_id], 3000, "freespace")
circle_boundary = gmsh.model.getBoundary((2,circle_id))
circle_boundary = abs.(getfield.(circle_boundary, 2))
gmsh.model.addPhysicalGroup(1, circle_boundary, 2000, "pec")
gmsh.model.mesh.generate(2)
gmsh.write("PECCircleAngle.msh")
# gmsh.fltk.run()

#-----------------------------------------------------------------------------------
# Next we will get the jacobian of a single (flat) triangle element at one point
# within the element.
# The assertion below allows us to conclude that the third column of the jacobian
# matrix is equal to the surface normal to the triangle.
#-----------------------------------------------------------------------------------

tri_element_tags, tri_element_nodes = gmsh.model.mesh.getElementsByType(tri_type) 
tri_element_nodes = convert.(Int64, tri_element_nodes)
tri_element_nodes = reshape(tri_element_nodes, 3, length(tri_element_nodes) ÷ 3)
n_tri = length(tri_element_tags)

# we will just take the barycentre
local_point = [1/3, 1/3, 0]
tri_index = 1
jacobian, determinant, global_points = gmsh.model.mesh.getJacobian(tri_element_tags[tri_index], local_point[:])
jacobian = reshape(jacobian, 3, 3)

inv_jacobian_trans = inv(jacobian')

a1 = jacobian[:,1]
a2 = jacobian[:,2]
a3 = jacobian[:,3]

nhat = cross(a1, a2)
nhat = nhat/norm(nhat)

@assert(norm(a3 - nhat) < 1e-12) 
@assert(norm(inv_jacobian_trans[:,3] - a3) < 1e-12)

