module EIL

#FEMMesh
export Mesh, setupGmsh, readMesh, setupMeshConnectivity, MeshConnectivity

#FEMConstitutives
export SimpleRealConstitutives, SimpleComplexConstitutives, ComplexConstitutives

#FEMBoundaryConditions
export  BoundaryCondition

#FEMBasis
export setupLocalBasis, LocalBasis, setupDOFMap

#FEMLocalMatrices
export elementMassAndStiffnessMatrices

include("FEMMesh.jl")
include("FEMConstitutives.jl")
include("FEMBasis.jl")
include("FEMBoundaryConditions.jl")
include("FEMLocalMatrices.jl")
# include("FEMGlobalMatrices.jl")
# include("FEMSources.jl")


end