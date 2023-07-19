using Revise

if isdefined(@__MODULE__, :EIL)
    
else
    include("src/EIL.jl")
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
# 1) What the jacobians of a line embedded in 3D space return from gmsh
#
# The answer to this question appears to be:
# If we let the Jacobian matrix be [a1, a2, a3] as three columns, then
# - a1 is directed along the line (as expected)
# - a2 is related to cross(a1, z-hat or x-hat) depending on the components of a1 
# - a3 = cross(a1, a2) completing a right-handed coordinate system along the line.
# 
# We also note that the vectors a2 and a3 are also unit vectors, which implies that they
# can be found in the corresponding columns of inv(J')
# ==================================================================================

#basic types for this unit test
line_type = gmsh.model.mesh.getElementType("line", 1)

#-----------------------------------------------------------------------------------
# To answer question #1 we will start by creating a line
#-----------------------------------------------------------------------------------

zaxis = [54, 4, 5]
zaxis_minus_z = zaxis - [0, 0, zaxis[3]]
p1 = [0, 0, 0]
p2 = p1 + zaxis
p1tag = 1
lc = 0.3
p1tag = gmsh.model.geo.addPoint(p1[1], p1[2], p1[3],lc)
p2tag = gmsh.model.geo.addPoint(p2[1], p2[2], p2[3],lc)
ltag = gmsh.model.geo.addLine(p1tag,p2tag)
gmsh.model.geo.synchronize();
gmsh.model.addPhysicalGroup(1, [ltag], 1000, "pec")
gmsh.model.mesh.generate(2)
gmsh.write("PECLineAngle.msh")
gmsh.fltk.run()

#-----------------------------------------------------------------------------------
# Next we will get the jacobian of a single (flat) triangle element at one point
# within the element
#-----------------------------------------------------------------------------------

line_element_tags, line_element_nodes = gmsh.model.mesh.getElementsByType(line_type) 
line_element_nodes = convert.(Int64, line_element_nodes)
println(line_element_nodes)
readline()
line_element_nodes = reshape(line_element_nodes, 2, length(line_element_nodes) ÷ 2)
n_line = length(line_element_tags)

# get the jacobian along the line
local_point = [0, 0, 0] # (0,0,0) is the barycenter
line_index = 3
jacobian, determinant, global_points = gmsh.model.mesh.getJacobian(line_element_tags[line_index], local_point[:])
jacobian = reshape(jacobian, 3, 3)
inv_jacobian_trans = inv(jacobian')

a1 = jacobian[:,1]

xhat = [1, 0, 0]
yhat = [0, 1, 0]
zhat = [0, 0, 1]

@assert(norm(a1/norm(a1) - (p2 - p1)/norm(p2-p1)) < 1e-12) #works for straight lines and tells us that a1 is tangential to the line

# According to the gmsh source code a1 is computed as follows:
# if the x or y components are the largest components in a1 (by magnitude) then a2 = a1 x zhat (normalized)
# otherwise a2 = a1 x xhat (normalized)


#a2 is either a1 x zhat or a1 x xhat with the decision based on if a1 has large x/y components compared to z
if (abs(a1[1]) >= abs(a1[2]) && abs(a1[1]) >= abs(a1[3])) || (abs(a1[1]) >= abs(a1[1]) && abs(a1[1]) >= abs(a1[3]))
    a1_cross_zhat = cross(a1, [0,0,1])
    a1_cross_zhat_normalized = a1_cross_zhat/norm(a1_cross_zhat)
    a2 = a1_cross_zhat_normalized
else
    a1_cross_xhat = cross(a1, [1,0,0])
    a1_cross_xhat_normalized = a1_cross_xhat/norm(a1_cross_xhat)
    a2 = a1_cross_xhat_normalized
end  

a3 = cross(a1, a2)/norm(cross(a1,a2)) #a1 x a2 = a3
    
@assert(norm(a1 - jacobian[:,1]) < 1e-12)
@assert(norm(a2 - jacobian[:,2]) < 1e-12)
@assert(norm(a3 - jacobian[:,3]) < 1e-12) 
@assert(norm(a2 - inv_jacobian_trans[:,2]) < 1e-12)
@assert(norm(a3 - inv_jacobian_trans[:,3]) < 1e-12)

println("jacobian[:,1] = ", jacobian[:,1])
println("jacobian[:,2] = ", jacobian[:,2])
println("jacobian[:,3] = ", jacobian[:,3])
println("a1 = ", a1)
println("a2 = ", a2)
println("a3 = ", a3)

