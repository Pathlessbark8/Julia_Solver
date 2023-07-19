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


#create the circle
circle_id = gmsh.model.occ.addCircle(0,0,0, 1, -1, 0, 2*pi);

#define the order of mesh and dist between points in the mesh
lc1 = 0.2
order = 2

#generate the mesh
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

groups = gmsh.model.getPhysicalGroups()
println(groups)
# set the order of the mesh
gmsh.model.mesh.setOrder(order)

# define meshing element
elementType = gmsh.model.mesh.getElementType("triangle", order)

#Check meshed elements properties
elementProperties=gmsh.model.mesh.getElementProperties(elementType)

#Get Edge Nodes in the mesh and create the edges
edgeNodes = gmsh.model.mesh.getElementEdgeNodes(elementType)
primaryEdgeNodes = gmsh.model.mesh.getElementEdgeNodes(elementType, -1, true)
gmsh.model.mesh.createEdges()

# println("primaryEdgeNodes: ", primaryEdgeNodes)

#Retrive the edge tags and Element Tags
edgeTags,_ = gmsh.model.mesh.getEdges(primaryEdgeNodes)
elementTags,_ = gmsh.model.mesh.getElementsByType(elementType)

#Perform Resizing of Arrays to Edge wise ordering
num_elements = length(elementTags)
edgeNodesPerElement = length(edgeNodes) ÷ num_elements
edgeNodes = reshape(edgeNodes, (edgeNodesPerElement ÷ 3, 3, num_elements))
edgeTags = reshape(edgeTags, (3, num_elements))
edges2Elements = Dict()

maxElementTag = gmsh.model.mesh.getMaxElementTag()

# println("Edge Nodes are ", edgeNodes)
# println("Max Element Tag: ", maxElementTag)
# tt= edgeTags 
# tt = convert.(Int64, tt)
# println("Edge Tags are ", tt )
#Create the edges2Elements dictionary
for i in 1:num_elements
    for j in 1:3
        if !haskey(edges2Elements, edgeTags[j,i] + maxElementTag)
            edges2Elements[edgeTags[j,i] ] = [elementTags[i]]
        else
            push!(edges2Elements[edgeTags[j,i] ],(elementTags[i]))
        end
    end
end

s = gmsh.model.addDiscreteEntity(1)


#Identify and assign a unique tag to each new edge
# maxElementTag = gmsh.model.mesh.getMaxElementTag()
uniqueEdgeTags = Set()
tagsForLines = []
EdgeNodesForLines = []
for i in eachindex(elementTags)
    for j in 1:3
        edge = edgeTags[j, i]
        if !in(edge, uniqueEdgeTags)
            push!(uniqueEdgeTags, edge)
            push!(tagsForLines, edge )
            append!(EdgeNodesForLines,edgeNodes[:,j,i])
        end
    end
end

#specficy the order and shape of the Elements being added
elementType1D = gmsh.model.mesh.getElementType("Line", order)

#Add the correspinding Elements to the mesh
gmsh.model.mesh.addElementsByType(s, elementType1D, tagsForLines,EdgeNodesForLines)


#Add the physical groups to the mesh
# uniqueEdgeTags = collect(uniqueEdgeTags)
# uniqueEdgeTags=convert.(Int64, uniqueEdgeTags)
# println(uniqueEdgeTags)
# println(typeof(uniqueEdgeTags))

# gmsh.model.geo.addPhysicalGroup(1, tagsForLines , 1000, "Edges")

# groups = gmsh.model.getPhysicalGroups()
# println("Groups are ", groups)
#Save the mesh and display it
gmsh.write("CircleMeshInitialDG.msh")

# gmsh.fltk.run()
gmsh.finalize()
