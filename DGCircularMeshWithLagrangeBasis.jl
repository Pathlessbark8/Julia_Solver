import Gmsh:gmsh #Julia package available for API that was installed
import LinearAlgebra
using CairoMakie

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

#define the setMeshOrder of mesh and dist between points in the mesh
lc1 = 1
setMeshOrder=2

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

#set the setMeshOrder of the mesh
gmsh.model.mesh.setOrder(setMeshOrder)

#define meshing element
elementType = gmsh.model.mesh.getElementType("triangle", setMeshOrder)

#Check meshed elements properties
elementProperties=gmsh.model.mesh.getElementProperties(elementType)

#Get Edge Nodes in the mesh and create the edges
edgeNodes = gmsh.model.mesh.getElementEdgeNodes(elementType)
primaryEdgeNodes = gmsh.model.mesh.getElementEdgeNodes(elementType, -1, true)
gmsh.model.mesh.createEdges()

#Retrive the edge tags and Element Tags
edgeTags,_ = gmsh.model.mesh.getEdges(primaryEdgeNodes)
elementTags,_ = gmsh.model.mesh.getElementsByType(elementType)

#Perform Resizing of Arrays to Edge wise ordering
num_elements = length(elementTags)
edgeNodesPerElement = length(edgeNodes) ÷ num_elements
edgeNodes = reshape(edgeNodes, (edgeNodesPerElement ÷ 3, 3, num_elements))
edgeTags = reshape(edgeTags, (3, num_elements))
edges2Elements = Dict()

#Create the edges2Elements dictionary
for i in 1:num_elements
    for j in 1:3
        if !haskey(edges2Elements, edgeTags[j,i])
            edges2Elements[edgeTags[j,i]] = elementTags[i]
        else
            edges2Elements[edgeTags[j,i]] = edges2Elements[edgeTags[j,i]]*(elementTags[i])
        end
    end
end
s = gmsh.model.addDiscreteEntity(1)


#Identify and assign a unique tag to each new edge
maxElementTag = gmsh.model.mesh.getMaxElementTag()
uniqueEdgeTags = Set()
tagsForLines = []
EdgeNodesForLines = []
for i in eachindex(elementTags)
    for j in 1:3
        edge = edgeTags[j, i]
        if !in(edge, uniqueEdgeTags)
            push!(uniqueEdgeTags, edge)
            push!(tagsForLines, edge + maxElementTag)
            append!(EdgeNodesForLines,edgeNodes[:,j,i])
        end
    end
end

#specficy the setMeshOrder and shape of the Elements being added
elementType1D = gmsh.model.mesh.getElementType("Line", setMeshOrder)

#Add the correspinding Elements to the mesh
gmsh.model.mesh.addElementsByType(s, elementType1D, tagsForLines,EdgeNodesForLines)

elementType1D = gmsh.model.mesh.getElementType("Line", setMeshOrder)

#Save the mesh and display it
gmsh.write("CircleMeshVeryFine.msh")

interpolationOrder = 3

function pp(label, v, mult)
    println(" * " * string(length(v) / mult) * " " * label * ": " * string(v))
end

# Iterate over all the element types present in the mesh:

# Retrieve properties for the given element type
elementName, dim, order, numNodes, numPrimNodes, localNodeCoord =
    gmsh.model.mesh.getElementProperties(elementType1D)
println("\n** " * elementName * " **\n")
println("numPrimNodes: " * string(numPrimNodes))
x= -1:0.5:1
localCoords = zeros(3,length(x))
localCoords[1,:] = x
reshape(localCoords, (3*length(x)))
localCoords = vec(localCoords)
display(localCoords)
# println("Local Coords are " * string(localCoords))
# localCoords, weights =
# gmsh.model.mesh.getIntegrationPoints(t, "Gauss" * string(interpolationOrder))
# println("Weight are " * string(weights))
# pp("integration points to integrate order " *
#    string(interpolationOrder) * " polynomials", localCoords, 3)
println("Local Coords are " * string(localCoords))
numComponents, basisFunctions, numOrientations = gmsh.model.mesh.getBasisFunctions(elementType1D, localCoords, "Lagrange$setMeshOrder")
println("Number of components are " * string(numComponents))

print(length(basisFunctions))
basisFunctions = reshape(basisFunctions, ((setMeshOrder+1),length(basisFunctions) ÷ (setMeshOrder+1)))
pp("basis functions at integration points", basisFunctions, 1)
f=Figure()
ax=Axis(f[-1,1])
println(typeof(x))
println(typeof(basisFunctions))
k = length(basisFunctions) ÷ length(x)
println("Size if basisFunctions is " * string(size(basisFunctions)))
for i in 1:k
    lines!(ax,x,basisFunctions[i,:])
end
f
sc = display(f);
# save("LagrangeBasiswithOrder$setMeshOrder.png",f)




# # Iterate over all the element types present in the mesh:
# elementTypes = gmsh.model.mesh.getElementTypes()

# for t in elementTypes
#     # Retrieve properties for the given element type
#     elementName, dim, order, numNodes, numPrimNodes, localNodeCoord =
#         gmsh.model.mesh.getElementProperties(t)
#     println("\n** " * elementName * " **\n")
#     println("numPrimNodes: " * string(numPrimNodes))


#     localCoords, weights =
#         gmsh.model.mesh.getIntegrationPoints(t, "Gauss" * string(interpolationOrder))
#     println("Weight are " * string(weights))
#     pp("integration points to integrate order " *
#        string(interpolationOrder) * " polynomials", localCoords, 3)

#     println("TPYE OF localCoords is " * string(typeof(localCoords)))

#     numComponents, basisFunctions, numOrientations =
#         # gmsh.model.mesh.getBasisFunctions(t, localCoords, "H1Legendre$setMeshOrder")
#         gmsh.model.mesh.getBasisFunctions(t, localCoords, "Lagrange")
#     println("Number of components are " * string(numComponents))
#     pp("basis functions at integration points", basisFunctions, 1)
#     # jacobians, determinants, coords =
#     #     gmsh.model.mesh.getJacobian(1, localCoords)
#     # pp("Jacobian determinants at integration points", determinants, 1)
# end

gmsh.finalize()
