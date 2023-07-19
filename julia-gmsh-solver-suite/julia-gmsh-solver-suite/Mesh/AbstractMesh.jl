#=-------------------------------------------------------------------
 Nodes
-------------------------------------------------------------------=#

struct Nodes
    tags::Vector
    coords::Matrix
end
numNodes(nodes::Nodes) = length(nodes.tags)

#=-------------------------------------------------------------------
 Elements
-------------------------------------------------------------------=#

abstract type AbstractElement end

struct Line <: AbstractElement end
struct Triangle <: AbstractElement end
struct Tetrahedron <: AbstractElement end

struct Elements
    tags::Vector
    elementNodes::Matrix
end
numElements(elements::Elements) = length(elements.tags)

#=-------------------------------------------------------------------
 Mesh
-------------------------------------------------------------------=#

abstract type AbstractMesh end