module Meshes

include("AbstractMesh.jl")
export AbstractMesh, AbstractElement, Line, Triangle, Tetrahedron, Nodes, Elements

include("GmshMeshUtils.jl")
export createEdges, createFaces, ElementProperties, Properties

include("SimpleMesh.jl")
export SimpleMesh, SimpleMeshConnectivity, loadMesh, setupMeshConnectivity, setTracesForEdges

end
