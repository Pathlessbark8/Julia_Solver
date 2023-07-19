function createCube(lcar1)
	L = 1 # dimension of box
	R = 0.1 # radius of sphere in box
	h = 0.1 # spacing between sphere

	gmsh.model.occ.addBox(0, 0, 0, L, L, L, 1)
	gmsh.model.occ.addSphere(L/2 - R - h/2, L/2, L/2, R, 2)
	gmsh.model.occ.addSphere(L/2 + R + h/2, L/2, L/2, R, 3)

	gmsh.model.occ.cut([(3,1)], [(3,2), (3,3)], 4, false, false)

	gmsh.model.occ.synchronize()

	# get boundary for original cube, not the cut cube, then delete original cube
	boundary = gmsh.model.getBoundary([(3,1)])
	gmsh.model.removeEntities([(3,1)])

	# create physical groups for volumes
	gmsh.model.addPhysicalGroup(3, [4], 3000)
	gmsh.model.addPhysicalGroup(3, [2], 3001)
	gmsh.model.addPhysicalGroup(3, [3], 3002)

	# physical groups for external surfaces
	for (i,b) in enumerate(boundary)
		gmsh.model.addPhysicalGroup(b[1], [b[2]], 2000+i)
	end

	ov = gmsh.model.getEntities(0)
	gmsh.model.mesh.setSize(ov, lcar1)
	gmsh.model.mesh.generate(3)
	
end

function getElements(dim, group)
	entities = gmsh.model.getEntitiesForPhysicalGroup(dim, group)

	elementTags = Vector{UInt64}(undef, 0)
	elementNodeTags = Vector{UInt64}(undef, 0)

	for e in entities
		elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, e)
		
		# we expect only a single element type
		@assert size(elemTypes, 1) == 1
		
		append!(elementTags, elemTags[1])
		append!(elementNodeTags, elemNodeTags[1])
	end

	nElem = size(elementTags, 1)
	nNodes = size(elementNodeTags, 1)
	nodesPerElem = nNodes ÷ nElem

	elementNodeTags = reshape(elementNodeTags, (nodesPerElem, nElem))
	return elementTags, elementNodeTags
end

function createPartitions(N)
	gmsh.option.setNumber("Mesh.PartitionCreateTopology", 1)
	gmsh.option.setNumber("Mesh.PartitionCreateGhostCells", 0)
	gmsh.option.setNumber("Mesh.PartitionCreatePhysicals", 1)
	gmsh.model.mesh.partition(N)
end

function getPartitions()
	N = gmsh.model.getNumberOfPartitions()
	partitions = [Vector{Tuple}(undef, 0) for n=1:N]
	
	entities = gmsh.model.getEntities()
	
	for e in entities
		parts = gmsh.model.getPartitions(e[1], e[2])
		if length(parts) > 0
			for p in parts
				push!(partitions[p], e)
			end
		end
	end
	return partitions
end

