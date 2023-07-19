#=-------------------------------------------------------------------
 Dependencies
-------------------------------------------------------------------=#

using LinearAlgebra
using StaticArrays
using SparseArrays
using Plots

using Revise

includet("Mesh/Mesh.jl")
using .Meshes
Revise.track(Meshes, "Mesh/AbstractMesh.jl")
Revise.track(Meshes, "Mesh/GmshMeshUtils.jl")
Revise.track(Meshes, "Mesh/SimpleMesh.jl")

includet("FunctionSpaces/FunctionSpaces.jl")
using .FunctionSpaces
Revise.track(FunctionSpaces, "FunctionSpaces/HOneElement.jl")
Revise.track(FunctionSpaces, "FunctionSpaces/HOneSpace.jl")
Revise.track(FunctionSpaces, "FunctionSpaces/HCurlElement.jl")
Revise.track(FunctionSpaces, "FunctionSpaces/HCurlSpace.jl")

includet("Physics/Physics.jl")
using .Physics
Revise.track(Physics, "Physics/AbstractPhysics.jl")
Revise.track(Physics, "Physics/EMPhysics.jl")


includet("Sources/Sources.jl")
using .Sources
Revise.track(Sources, "Sources/AbstractSources.jl")
Revise.track(Sources, "Sources/EMSources.jl")

includet("Solvers/Solvers.jl")
using .Solvers
Revise.track(Solvers, "Solvers/AbstractSolvers.jl")
Revise.track(Solvers, "Solvers/EMSolverUtilities.jl")
Revise.track(Solvers, "Solvers/EMTMFEMSolver.jl")

includet("Extras/Visualization.jl")
using .Visualization
Revise.track(Visualization, "Extras/MakieViz.jl")

using Gmsh: gmsh
using SparseArrays
using StaticArrays
using LinearAlgebra

include("Extras/Analytic.jl")

using GLMakie


include("src/EIL.jl")
using .EIL

if gmsh.isInitialized() == 1
    gmsh.finalize()
end
gmsh.initialize()

ϵ0 = 8.85418782 * 10-12 #permittivity of free space
μ0 = 1.25663706 * 10-6 #permeability of free space

function setInitialConditions(elementTags)
    element_relative_permeability = ones(Float64,length(elementTags))
    element_relative_permittivity = ones(Float64,length(elementTags))
    for i in 1:length(elementTags)
        element_relative_permeability[i] = μ0
        element_relative_permittivity[i] = ϵ0
    end
    return element_relative_permeability,element_relative_permittivity
end


function setUpLocalMassMatrixTM(mesh,order,edgeTags,elementTags,Elements2Nbhs)
    referenceElement = FunctionSpaces.HOneElement(2,order)
    num_elements=length(elementTags)
    # (num_components, num_functions, num_points, num_orientations) = size(referenceElement.gradient_functions)
    # basis_functions=zeros(Float64,num_functions,num_points,3,num_elements)
    # basis_gradient_functions=Array{Float64}(undef,num_components,num_functions,num_points,num_orientations,3,num_elements)   

    println("Here ",size(referenceElement.functions))
    (num_bf_components, num_functions, num_points, num_bf_orientations) = size(referenceElement.functions)
    (num_gbf_components, num_functions, num_points, num_gbf_orientations) = size(referenceElement.gradient_functions)
    
    local_functions = referenceElement.functions
    gradient_functions = referenceElement.gradient_functions
    gradient_functions = reinterpret(SVector{num_gbf_components,Float64}, gradient_functions)
    gradient_functions = reshape(gradient_functions, (num_functions, num_points, num_gbf_orientations))
    # for i in 1:num_elements
    #     for j in 1:3
    #         basis_functions[:,:,:,j,i]=referenceElement.functions[:,:,:,:]
    #         # basis_gradient_functions[:,:,:,:,j,i]=referenceElement.gradient_functions[:,:,:,:]
    #     end
    # end

    basis_functions = zeros(num_points, num_functions)
    gradient_basis_functions = zeros(SVector{num_gbf_components,Float64}, num_points, num_functions)

    # for i in 1:num_elements
    #     for j in 1:3
    #         jacobians,determinants,coord = gmsh.model.mesh.getJacobian(edgeTags[j,i],referenceElement.quad_points)
    #         jacobians = reinterpret(SMatrix{3,3,Float64,9}, jacobians)
    #         orient = 1 + gmsh.model.mesh.getBasisFunctionsOrientationForElement(edgeTags[j,i], "Lagrange")
    #         for k = 1:num_points, l = 1:num_functions
    #             basis_functions[k, l] = local_functions[:,l,k,orient]
    #             gradient_basis_functions[k, l] =  jacobians[k]' \ gradient_functions[l, k, orient]
    #         end
    #     end
    # end

    M = zeros(num_functions, num_functions, num_elements)
    Sx = zeros(num_functions, num_functions, num_elements)
    Sy = zeros(num_functions, num_functions, num_elements)

    GlobalXYCoordinates::Array{Float64} = Array{Float64}(undef,2,num_points*num_elements)
    println(typeof(GlobalXYCoordinates))
    # asd
    yhat=[0; 1.0; 0]
    for i in 1:num_elements
        jacobians,determinants,coord = gmsh.model.mesh.getJacobian(elementTags[i],referenceElement.quad_points)
        coord=reshape(coord,(3,num_points))
        GlobalXYCoordinates[1,(i-1)*num_points+1:i*num_points]=coord[1,:]
        GlobalXYCoordinates[2,(i-1)*num_points+1:i*num_points]=coord[2,:]
        jacobians = reinterpret(SMatrix{3,3,Float64,9}, jacobians)
        orient = 1 + gmsh.model.mesh.getBasisFunctionsOrientationForElement(elementTags[i], "Lagrange$order")
        for k = 1:num_points, l = 1:num_functions
            basis_functions[k, l] = local_functions[1,l,k,orient]
            gradient_basis_functions[k, l] =  jacobians[k]' \ gradient_functions[l, k, orient]
        end
        for k in 1:num_functions
            for l in 1:num_functions
                tempM=0
                tempSx=0
                tempSy=0 
                for m in 1:num_points
                    tempM += referenceElement.quad_weights[m] * (basis_functions[m, l] * basis_functions[m, k]) * determinants[m]
                    tempSx += referenceElement.quad_weights[m] * (gradient_basis_functions[m, l][1]'*basis_functions[m, k]) * determinants[m]
                    tempSy += referenceElement.quad_weights[m] * (gradient_basis_functions[m, l][2]'*basis_functions[m, k]) * determinants[m]
                end 
                M[k,l,i]=tempM
                Sx[k,l,i]=tempSx
                Sy[k,l,i]=tempSy     
            end
        end
    end

    display(M[:,:,1])
    element_relative_permeability,element_relative_permittivity = setInitialConditions(elementTags)

    α = 1
    ω = 2*π*1e9

    edge_referenceElemenet = HOneElement(1,order)
    edge_num_points = length(edge_referenceElemenet.quad_weights)
    
    basis_function_type = gmsh.model.mesh.getElementType("triangle", 1)
    edge_num_components, edge_local_functions, edge_num_orientations =
        gmsh.model.mesh.getBasisFunctions(basis_function_type, edge_referenceElemenet.quad_points, "Lagrange$order")

    edge_num_functions = length(edge_local_functions) ÷ edge_num_components ÷ edge_num_orientations ÷ edge_num_points
    edge_local_functions = reshape(edge_local_functions, (edge_num_components, edge_num_functions, edge_num_points, edge_num_orientations))

    @assert num_functions == edge_num_functions

    edge_basis_functions = zeros(SVector{floor(Int64, edge_num_components),Float64}, edge_num_points, edge_num_functions)

    Normals = zeros(3,edge_num_points,3,num_elements)
    AMinusGlobal = zeros(3*edge_num_functions,3*edge_num_functions,3,num_elements)    
    APlusGlobal = zeros(3*edge_num_functions,3*edge_num_functions,3,num_elements)

    println("Here 0")
    for i in 1:num_elements
        for j in 1:3
            currentNbh = Elements2Nbhs[elementTags[i]][j]
            if i==1
                println("Neighbour: ",currentNbh)
            end
            # if currentNbh == 0
            #     continue
            # end
            jacobians,determinants,coord = gmsh.model.mesh.getJacobian(elementTags[i],edge_referenceElemenet.quad_points)
            jacobians = reinterpret(SMatrix{3,3,Float64,9}, jacobians)
            orient = 1 + gmsh.model.mesh.getBasisFunctionsOrientationForElement(edgeTags[j,i], "Lagrange$order")
            for k = 1:edge_num_points, l = 1:edge_num_functions
                edge_basis_functions[k, l] = edge_local_functions[:,l,k,orient]
            end
            for k in 1:edge_num_points
                Normals[:,k,j,i]= jacobians[k]' \ yhat
            end
        end
    end
    for i in 1:num_elements
        for j in 1:3
            currentNbh = Elements2Nbhs[elementTags[i]][j]

            jacobians,determinants,coord = gmsh.model.mesh.getJacobian(elementTags[i],edge_referenceElemenet.quad_points)
            jacobians = reinterpret(SMatrix{3,3,Float64,9}, jacobians)
            orient = 1 + gmsh.model.mesh.getBasisFunctionsOrientationForElement(edgeTags[j,i], "Lagrange$order")
            for k = 1:edge_num_points, l = 1:edge_num_functions
                edge_basis_functions[k, l] = edge_local_functions[:,l,k,orient]
            end

            currentNbhIndx = findall(x -> x == currentNbh, elementTags)
            if currentNbhIndx != []
                currentNbhIndx = currentNbhIndx[1]
            end
            Zplus=sqrt(element_relative_permeability[i]/element_relative_permittivity[i])
            if currentNbh == 0
                Zminus=Zplus
            else
                Zminus=sqrt(element_relative_permeability[currentNbhIndx]/element_relative_permittivity[currentNbhIndx])
            end

            Yplus = 1/Zplus
            Yminus = 1/Zminus 

            reciprocalOfZSum = 1/(Zplus + Zminus)
            reciprocalOfYSum = 1/(Yplus + Yminus)

            AMinusLocal = ones(3*edge_num_functions,3*edge_num_functions)
            for k in 1:3*edge_num_functions
                for l in 1:3*edge_num_functions
                    for m in 1:edge_num_points
                        nx=Normals[1,m,j,i]
                        ny=Normals[2,m,j,i]
                        if(k<=edge_num_functions && l<=edge_num_functions)
                            AMinusLocal[k,l] += reciprocalOfZSum*α*edge_referenceElemenet.quad_weights[m]*(-nx^2-ny^2)*edge_basis_functions[m,k][1]'*edge_basis_functions[m,l][1]*determinants[m]
                        elseif(k<=edge_num_functions && l<=2*edge_num_functions)
                            AMinusLocal[k,l] += -Zplus*reciprocalOfZSum*edge_referenceElemenet.quad_weights[m]*(-ny)*edge_basis_functions[m,k][1]'*edge_basis_functions[m,l-edge_num_functions][1]*determinants[m]
                        elseif(k<=edge_num_functions && l<=3*edge_num_functions)
                            AMinusLocal[k,l] += -Zplus*reciprocalOfZSum*edge_referenceElemenet.quad_weights[m]*(nx)*edge_basis_functions[m,k][1]'*edge_basis_functions[m,l-2*edge_num_functions][1]*determinants[m]
                        elseif(k<=2*edge_num_functions && l<=edge_num_functions)
                            AMinusLocal[k,l] += Yplus*reciprocalOfYSum*edge_referenceElemenet.quad_weights[m]*(ny)*edge_basis_functions[m,k-edge_num_functions][1]'*edge_basis_functions[m,l][1]*determinants[m]
                        elseif(k<=2*edge_num_functions && l<=2*edge_num_functions)
                            AMinusLocal[k,l] += reciprocalOfYSum*α*edge_referenceElemenet.quad_weights[m]*(-ny^2)*edge_basis_functions[m,k-edge_num_functions][1]'*edge_basis_functions[m,l-edge_num_functions][1]*determinants[m]
                        elseif(k<=2*edge_num_functions && l<=3*edge_num_functions)
                            AMinusLocal[k,l] += reciprocalOfYSum*α*edge_referenceElemenet.quad_weights[m]*(nx*ny)*edge_basis_functions[m,k-edge_num_functions][1]'*edge_basis_functions[m,l-2*edge_num_functions][1]*determinants[m]
                        elseif(k<=3*edge_num_functions && l<=edge_num_functions)
                            AMinusLocal[k,l] += Yplus*reciprocalOfYSum*edge_referenceElemenet.quad_weights[m]*(-nx)*edge_basis_functions[m,k-2*edge_num_functions][1]'*edge_basis_functions[m,l][1]*determinants[m]
                        elseif(k<=3*edge_num_functions && l<=2*edge_num_functions)
                            AMinusLocal[k,l] += reciprocalOfYSum*α*edge_referenceElemenet.quad_weights[m]*(nx*ny)*edge_basis_functions[m,k-2*edge_num_functions][1]'*edge_basis_functions[m,l-edge_num_functions][1]*determinants[m]
                        elseif(k<=3*edge_num_functions && l<=3*edge_num_functions)
                            AMinusLocal[k,l] += reciprocalOfYSum*α*edge_referenceElemenet.quad_weights[m]*(-nx^2)*edge_basis_functions[m,k-2*edge_num_functions][1]'*edge_basis_functions[m,l-2*edge_num_functions][1]*determinants[m]
                        end
                    end
                end
            end
            AMinusGlobal[:,:,j,i] = AMinusLocal

            APlusLocal = ones(3*edge_num_functions,3*edge_num_functions)
            if currentNbh != 0
                for k in 1:3*edge_num_functions
                    for l in 1:3*edge_num_functions
                        for m in 1:edge_num_points
                            nx=Normals[1,m,j,currentNbhIndx]
                            ny=Normals[2,m,j,currentNbhIndx]
                            if(k<=edge_num_functions && l<=edge_num_functions) 
                                APlusLocal[k,l] += -reciprocalOfZSum*α*edge_referenceElemenet.quad_weights[m]*(-nx^2-ny^2)*edge_basis_functions[m,k][1]'*edge_basis_functions[m,l][1]*determinants[m]
                            elseif(k<=edge_num_functions && l<=2*edge_num_functions)
                                APlusLocal[k,l] += Zplus*reciprocalOfZSum*edge_referenceElemenet.quad_weights[m]*(-ny)*edge_basis_functions[m,k][1]'*edge_basis_functions[m,l-edge_num_functions][1]*determinants[m]
                            elseif(k<=edge_num_functions && l<=3*edge_num_functions)
                                APlusLocal[k,l] += Zplus*reciprocalOfZSum*edge_referenceElemenet.quad_weights[m]*(nx)*edge_basis_functions[m,k][1]'*edge_basis_functions[m,l-2*edge_num_functions][1]*determinants[m]
                            elseif(k<=2*edge_num_functions && l<=edge_num_functions)
                                APlusLocal[k,l] += -Yplus*reciprocalOfYSum*edge_referenceElemenet.quad_weights[m]*(ny)*edge_basis_functions[m,k-edge_num_functions][1]'*edge_basis_functions[m,l][1]*determinants[m]
                            elseif(k<=2*edge_num_functions && l<=2*edge_num_functions)
                                APlusLocal[k,l] += -reciprocalOfYSum*α*edge_referenceElemenet.quad_weights[m]*(-ny^2)*edge_basis_functions[m,k-edge_num_functions][1]'*edge_basis_functions[m,l-edge_num_functions][1]*determinants[m]
                            elseif(k<=2*edge_num_functions && l<=3*edge_num_functions)
                                APlusLocal[k,l] += -reciprocalOfYSum*α*edge_referenceElemenet.quad_weights[m]*(nx*ny)*edge_basis_functions[m,k-edge_num_functions][1]'*edge_basis_functions[m,l-2*edge_num_functions][1]*determinants[m]
                            elseif(k<=3*edge_num_functions && l<=edge_num_functions)
                                APlusLocal[k,l] += -Yplus*reciprocalOfYSum*edge_referenceElemenet.quad_weights[m]*(-nx)*edge_basis_functions[m,k-2*edge_num_functions][1]'*edge_basis_functions[m,l][1]*determinants[m]
                            elseif(k<=3*edge_num_functions && l<=2*edge_num_functions)
                                APlusLocal[k,l] += -reciprocalOfYSum*α*edge_referenceElemenet.quad_weights[m]*(nx*ny)*edge_basis_functions[m,k-2*edge_num_functions][1]'*edge_basis_functions[m,l-edge_num_functions][1]*determinants[m]
                            elseif(k<=3*edge_num_functions && l<=3*edge_num_functions)
                                APlusLocal[k,l] += -reciprocalOfYSum*α*edge_referenceElemenet.quad_weights[m]*(-nx^2)*edge_basis_functions[m,k-2*edge_num_functions][1]'*edge_basis_functions[m,l-2*edge_num_functions][1]*determinants[m]
                            end
                        end
                    end
                end
                APlusGlobal[:,:,j,i] = APlusLocal
            end
        end
    end

    display(APlusGlobal[:,:,2,1])
    ElementMassMatrixComplex = zeros(ComplexF64,3*num_functions,3*num_functions,num_elements) 
    for i in 1:num_elements
        ElementMassMatrixComplex[1:num_functions,1:num_functions,i] = element_relative_permittivity[i] .* M[:,:,i]
        ElementMassMatrixComplex[num_functions+1:2*num_functions,num_functions+1:2*num_functions,i] = element_relative_permeability[i] .* M[:,:,i]
        ElementMassMatrixComplex[2*num_functions+1:3*num_functions,2*num_functions+1:3*num_functions,i] = element_relative_permeability[i] .* M[:,:,i]
    end
    println("Here 1")
    # ElementMassMatrixComplex = reshape(Complex.(ElementMassMatrixSimple),3*num_functions,3*num_functions,num_elements)
    ElementMassMatrixComplex = ElementMassMatrixComplex * 1im *ω

    ElementStiffnessMatrix = zeros(3*num_functions,3*num_functions,num_elements)
    for i in 1:num_elements
        ElementStiffnessMatrix[1:num_functions,num_functions+1:2*num_functions,i] = Sy[:,:,i]
        ElementStiffnessMatrix[1:num_functions,2*num_functions+1:3*num_functions,i] = -Sx[:,:,i]
        ElementStiffnessMatrix[num_functions+1:2*num_functions,1:num_functions,i] = Sy[:,:,i]
        ElementStiffnessMatrix[2*num_functions+1:3*num_functions,1:num_functions,i] = -Sx[:,:,i]
    end

    # println(ElementStiffnessMatrix[:,:,1])
    println("Here 2")

    # GlobalMainMatrix = spzeros(ComplexF64,3*num_functions*num_elements,3*num_functions*num_elements)
    I = Int64[]
    J = Int64[]
    V = ComplexF64[]
    for i in 1:num_elements
        tempMatrix = ElementMassMatrixComplex[:,:,i] + ElementStiffnessMatrix[:,:,i]
        # for j in 1:3
        #     tempMatrix += AMinusGlobal[:,:,j,i]
        # end
        for j in 1:3*num_functions
            for k in 1:3*num_functions
                if tempMatrix[j,k] != 0
                    push!(I,j+(i-1)*3*num_functions)
                    push!(J,k+(i-1)*3*num_functions)
                    push!(V,tempMatrix[j,k])
                end
            end
        end
        # GlobalMainMatrix[(i-1)*3*num_functions+1:i*3*num_functions,(i-1)*3*num_functions+1:i*3*num_functions] = tempMatrix
    end

    println("Here 3")
    GlobalMainMatrix = sparse(I,J,V)
    display(GlobalMainMatrix)
    # for i in 1:num_elements
    #     for j in 1:3
    #         currentNbh = Elements2Nbhs[elementTags[i]][j]
    #         if currentNbh == 0
    #             continue
    #         end
    #         currentNbhIndx = searchsortedfirst(elementTags, currentNbh)
    #         for k in 1:3*num_functions
    #             for m in 1:3*num_functions
    #                 if APlusGlobal[k,m,j,i] != 0
    #                     push!(I,k+(i-1)*3*num_functions)
    #                     push!(J,m+(currentNbhIndx-1)*3*num_functions)
    #                     push!(V,APlusGlobal[k,m,j,i])
    #                 end
    #             end
    #         end
    #         # GlobalMainMatrix[(i-1)*3*num_functions+1:i*3*num_functions,(currentNbhIndx-1)*3*num_functions+1:currentNbhIndx*3*num_functions] = APlusGlobal[:,:,j,i]
    #     end
    # end
    GlobalMainMatrix = sparse(I,J,V)
    println("Here 4")
    display(GlobalMainMatrix)
    source::EMTMGaussianCurrentDensity = EMTMGaussianCurrentDensity(0.0,0.0,0.5,0.5,1)

    Jz = evaluateSource(source,GlobalXYCoordinates,1.0)
    # Jz[1:num_functions] .= 1
    # Jz[num_functions+1:end] .=0
    SolutionMatrixGlobal = zeros(3*num_functions*num_elements)

    I = Int64[]
    # J =[]
    V = ComplexF64[]
    for i in 1:num_elements
        sumJz=0
        _,determinants,_ = gmsh.model.mesh.getJacobian(elementTags[i],referenceElement.quad_points)
        JzMatrix = zeros(num_functions)
        for j in 1:num_functions
            for k in 1:num_points
                JzMatrix[j] += determinants[k]*referenceElement.quad_weights[k]*basis_functions[k,j][1]*Jz[(i-1)*num_points+k]
            end
        end
        # for j in 1:num_functions
        #     push!(I,j+3*(i-1)*num_functions)
        #     push!(V,JzMatrix[j])
        # end
        SolutionMatrixGlobal[(3*(i-1)*num_functions+1):(3*(i-1)*num_functions+num_functions)] = JzMatrix
    end
    # SolutionMatrixGlobal = sparsevec(I,V)
    println("Here 5")
    # SolutionMatrixGlobal .= 0
    # SolutionMatrixGlobal[end] = 1
    Bb = GlobalMainMatrix \ SolutionMatrixGlobal
    display(maximum(abs.(SolutionMatrixGlobal)))
    # XYZCoordinates = Array{Float64}(undef,4*num_elements*num_points)
    # for i in 1:num_elements*num_points
    #     XYZtemp=Array{Float64}(undef,4)
    #     XYZtemp[1]=GlobalXYCoordinates[1,i]
    #     XYZtemp[2]=GlobalXYCoordinates[2,i]
    #     XYZtemp[3]=Jz[i]
    #     XYZtemp[4]=Jz[i]
    #     XYZCoordinates[4*(i-1)+1:4*i]=XYZtemp
    # end
    # t1 = gmsh.view.add("Points")
    # len::Int64=length(XYZCoordinates)/4
    # gmsh.view.addListData(t1, "SP", len,XYZCoordinates)
    # gmsh.fltk.run()

    println("Here 6")
    plot_points = 3
    XYZTriangleCoordinates = Array{Float64}(undef,4*plot_points*num_elements)
    localTriangleCoordinates = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0]
    # RandomXYZCoordinates = zeros(3 * plot_points)
    # RandomXYZCoordinates[1:3:3*plot_points] = rand(1, plot_points)
    # for i in 1:plot_points
    #     RandomXYZCoordinates[3*(i-1)+2] = (1-RandomXYZCoordinates[3*(i-1)+1])*rand()
    # end
    # println(RandomXYZCoordinates[1:12])
    type_number = gmsh.model.mesh.getElementType("triangle", 1)
    for i in 1:num_elements
        _,_,coord= gmsh.model.mesh.getJacobian(elementTags[i],localTriangleCoordinates)
        local_num_components, local_functions, local_num_orientations = gmsh.model.mesh.getBasisFunctions(type_number, localTriangleCoordinates, "Lagrange$order")
        local_num_points::Int64 = length(localTriangleCoordinates)/3
        @assert local_num_points == plot_points
        local_num_functions = length(local_functions) ÷ local_num_components ÷ local_num_orientations ÷ local_num_points
        local_functions = reshape(local_functions, (local_num_components, local_num_functions, local_num_points, local_num_orientations))
        local_basis_functions = zeros(SVector{floor(Int64, local_num_components),Float64}, local_num_points, local_num_functions)
        orient = 1 + gmsh.model.mesh.getBasisFunctionsOrientationForElement(elementTags[i], "Lagrange$order")
        for k = 1:local_num_points, l = 1:local_num_functions
            local_basis_functions[k, l] = local_functions[:,l,k,orient]
        end
        
        
        XYZTriangleCoordinates[4*plot_points*(i-1)+ 1:4*plot_points*(i-1) + plot_points]=coord[1:3:3*plot_points]
        XYZTriangleCoordinates[4*plot_points*(i-1)+ plot_points + 1: 4*plot_points*(i-1)+2*plot_points]=coord[2:3:3*plot_points]
        XYZTriangleCoordinates[4*plot_points*(i-1)+ 2*plot_points + 1:4*plot_points*(i-1)+3*plot_points]=coord[3:3:3*plot_points]
        # for j in 1:plot_points
        #     XYZTriangleCoordinates[4*plot_points*(i-1)+j] = coord[3*(j-1)+1:3*j]
        # end
        basis_functions_per_element = zeros(local_num_functions)

        for j in 1:local_num_points
            EzValue = 0
            for k in 1:local_num_functions
                basis_functions_per_element[k] = local_basis_functions[j,k][1]
            end
            EzValue = Bb[3*(i-1)*num_functions+1:3*(i-1)*num_functions+num_functions]'*basis_functions_per_element
            # println(EzValue)
            XYZTriangleCoordinates[4*plot_points*(i-1)+3*plot_points+j] = imag(EzValue)
        end
    #     # XYZTriangleCoordinates[]
    #     # XYZTriangleCoordinates[12*(i-1)+10:12*i] .= EzValue
    end
    # println(XYZTriangleCoordinates[1:24])
    # println(Bb[1:12])
    t1 = gmsh.view.add("Triangles")
    len1::Int64=length(XYZTriangleCoordinates)/12
    gmsh.view.addListData(t1, "ST", len1,XYZTriangleCoordinates)
    gmsh.fltk.run()
end


# problem geometry
mesh_filename = "./CircleMeshInitialDG.msh"
# mesh_filename = "./CircleMeshVeryFine.msh"

setMeshOrder = 2
problem_type = EIL.EMTMType()
basis_type = EIL.NodalBasisType()

mesh = loadMesh(mesh_filename)

edgeConectivity,edgeTags,elementTags,Elements2Nbhs = setTracesForEdges(mesh_filename,setMeshOrder)

# println(Elements2Nbhs)
# gmsh.fltk.run()
# for i in eachindex(elementTags)
#     println(i)
# end

@time setUpLocalMassMatrixTM(mesh,setMeshOrder,edgeTags,elementTags,Elements2Nbhs)






# print(local_basis)
# mesh_connectivity = setupMeshConnectivity(mesh)
