# -----------------------------------------------------------------------------
#
#  Gmsh Julia extended tutorial 3
#
#  Post-processing data import: list-based
#
# -----------------------------------------------------------------------------

import Gmsh:gmsh

gmsh.initialize(append!(["gmsh"], ARGS))

# Gmsh supports two types of post-processing data: "list-based" and
# "model-based". Both types of data are handled through the `view' interface.

# List-based views are completely independent from any model and any mesh: they
# are self-contained and simply contain lists of coordinates and values, element
# by element, for 3 types of fields (scalar "S", vector "V" and tensor "T") and
# several types of element shapes (point "P", line "L", triangle "T", quadrangle
# "Q", tetrahedron "S", hexahedron "H", prism "I" and pyramid "Y"). (See `x4.jl'
# for a tutorial on model-based views.)

# To create a list-based view one should first create a view:
t1 = gmsh.view.add("A list-based view")

# List-based data is then added by specifying the type as a 2 character string
# that combines a field type and an element shape (e.g. "ST" for a scalar field
# on triangles), the number of elements to be added, and the concatenated list
# of coordinates (e.g. 3 "x" coordinates, 3 "y" coordinates, 3 "z" coordinates
# for first order triangles) and values for each element (e.g. 3 values for
# first order scalar triangles, repeated for each step if there are several time
# steps).

# Let's create two triangles...
triangle1 = [0., 1., 1., # x coordinates of the 3 triangle nodes
             0., 0., 1., # y coordinates of the 3 triangle nodes
             0., 0., 0.] # z coordinates of the 3 triangle nodes

# ... and append values for 10 time steps
for step in 0:9
    push!(triangle1, 10., 11. - step, 12.)  # 3 node values for each step
end

# List-based data is just added by concatenating the data for all the triangles:
gmsh.view.addListData(t1, "SP", 1, triangle1)


if !("-nopopup" in ARGS)
    gmsh.fltk.run()
end

gmsh.finalize()
