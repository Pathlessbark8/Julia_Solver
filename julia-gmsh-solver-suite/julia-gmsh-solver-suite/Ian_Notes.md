# On Ordering of Global Points

Indexing global points for plotting purposes has turned out to be a bit of a pain, likely because I'm making it more complicated than it needs to be.

Every node in the gmsh mesh file comes with a tag. We index these nodes as 1:n_nodes, and can look up the node using the tag via the mesh.nodes.tag2idx dictionary. I suspect that these tags are sequential from 1:n_nodes, but can't be sure. Regardless, if we are looking up a node in the list, then we really want the idx. Elements point to the tag and not the index.

In order to create elements and/or basis points we loop over physical groups, entities, and element types (in that order). In our initial testing we only have one physical group, one entity and one element type. We pull elements in this order. We also pull jacobians (global nodal basis points) in this same order. We store these points as one big array:

global_points = [points on element 1, points on elements 2, points on element 3] and so on. This array has 3 rows, one for each coordinate of the points. 

We can get the points on element i as follows:

ele_index = ((point_index-1) $\div$ 3) + 1

and then loop over the number of points on that element. 

To map these points to the appropriate gmsh index (or tag) we need to compare all points on a given element, to the quad points on that element.