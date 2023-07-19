import Gmsh:gmsh #Julia package available for API that was installed
import LinearAlgebra

#using PyPlot
using Plots
#using GLMakie

#if the code crashes and gmsh doesn't get finalized, problems ensue on next run
try
    gmsh.finalize()
catch

end

#start up the gmsh api and set a few of my own GUI display preferences
gmsh.initialize()
gmsh.option.setNumber("Geometry.Volumes", 1)
gmsh.option.setNumber("Geometry.Surfaces", 1)
gmsh.option.setNumber("Geometry.Normals", 30)


circle_id = gmsh.model.occ.addCircle(0,0,0, 1, -1, 0, 2*pi);
circle_id
lc1 = 0.008
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
gmsh.write("CircleMeshVeryFine.msh")
gmsh.fltk.run()