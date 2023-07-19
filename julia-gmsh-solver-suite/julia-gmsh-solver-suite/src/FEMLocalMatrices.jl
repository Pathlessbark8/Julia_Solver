#=------------------------------------
Element Matrices Struct
-------------------------------------=#

abstract type AbstractElementMatrices end

struct ElementMatrices2DTM <: AbstractElementMatrices
    M::Array{Float64,3}                     #nbasis x nbasis x nelements
    Sx::Array{Float64,3}                    #nbasis x nbasis x nelements
    Sy::Array{Float64,3}                    #nbasis x nbasis x nelements
end

struct ElementMatrices2DTE <: AbstractElementMatrices
    Mx::Array{Float64,3}                    #nbasis x nbasis x nelements
    My::Array{Float64,3}                    #nbasis x nbasis x nelements
    S::Array{Float64,3}                     #nbasis x nbasis x nelements
end

struct ElementMatrices3D <: AbstractElementMatrices
    Mx::Array{Float64,3}                    #nbasis x nbasis x nelements
    My::Array{Float64,3}                    #nbasis x nbasis x nelements
    Mz::Array{Float64,3}                    #nbasis x nbasis x nelements
    Sx::Array{Float64,3}                     #nbasis x nbasis x nelements
    Sy::Array{Float64,3}                     #nbasis x nbasis x nelements
    Sz::Array{Float64,3}                     #nbasis x nbasis x nelements
end


#=------------------------------------
elementMassAndStiffnessMatrices() - 2D TM
-------------------------------------=#
function elementMassAndStiffnessMatrices(mesh::Mesh, element_type, basis::LocalBasis, problem_type::EMTMType)

    #I have a thought to make this tag-based so we can compute them in batches of elements - but can discuss that in the future
    M = Array{Float64, 3}(undef,  basis.n_basis, basis.n_basis, basis.n_elements) #nbasis x nbasis x nelements
    Sx = Array{Float64, 3}(undef, basis.n_basis, basis.n_basis, basis.n_elements) #nbasis x nbasis x nelements
    Sy = Array{Float64, 3}(undef, basis.n_basis, basis.n_basis, basis.n_elements) #nbasis x nbasis x nelements
    T = Array{Float64, 2}(undef, 3, basis.n_quad_points)

    W = diagm(basis.quad_weights)    #element-independent weight matrix
    dBx = Array{Float64, 2}(undef, basis.n_basis, basis.n_quad_points) #derivative of basis in x
    dBy = Array{Float64, 2}(undef, basis.n_basis, basis.n_quad_points) #derivative of basis in y
    
    #construct matrices for each element
    for iele in range(1,basis.n_elements)
        J = diagm(basis.J[:,iele]) #diagonal matrix of quad point jacobians for this element
        B = basis.basis[1,:,:,1]   #1st (only) component, n_basis x n_quad, 1st (only) orientation
        M[:,:,iele] = B*(J*W)*B'   #mass matrix is integral of basis functions times basis functions       
      
        #fancy way of computing the derivatives of basis functions at each quad point
        for ipoint in 1:basis.n_quad_points
            #jacobian^-T (at a point) * gradient of local basis = gradient in global coordinates (with a transpose or something)?
            T = basis.jacobians[:,:,ipoint,iele]'\basis.grad_basis[:,:,ipoint,1]
            dBx[:,ipoint] = T[1,:] #build these point by point
            dBy[:,ipoint] = T[2,:]
        end

        #once we have the derivatives of the basis functions at every point
        #the stiffness matrices are derivatives integrated with derivatives (second-order equations)
        #could add some additional options here if we don't want second order equations
        Sx[:,:,iele] = dBx*(J*W)*dBx'
        Sy[:,:,iele] = dBy*(J*W)*dBy'
    end

    #create the struct
    element_matrices = ElementMatrices2DTM(M, Sx, Sy)

    return element_matrices
end


#=------------------------------------
elementMassAndStiffnessMatrices() - 2D TE
-------------------------------------=#
function elementMassAndStiffnessMatrices(basis::LocalBasis, problem_type::EMTEType)

    #I have a thought to make this tag-based so we can compute them in batches of elements - but can discuss that in the future
    Mx = Array{Float64, 3}(undef,  basis.n_basis, basis.n_basis, basis.n_elements) #nbasis x nbasis x nelements
    My = Array{Float64, 3}(undef,  basis.n_basis, basis.n_basis, basis.n_elements) #nbasis x nbasis x nelements
    S = Array{Float64, 3}(undef, basis.n_basis, basis.n_basis, basis.n_elements) #nbasis x nbasis x nelements
    
    all_bases_at_point = Array{Float64, 2}(undef, 3, basis.n_basis)
    all_curl_bases_at_point = Array{Float64, 2}(undef, 3, basis.n_basis)

    W = diagm(basis.quad_weights)    #element-independent weight matrix
    Bx = Array{Float64, 2}(undef, basis.n_basis, basis.n_quad_points) #derivative of basis in x
    By = Array{Float64, 2}(undef, basis.n_basis, basis.n_quad_points) #derivative of basis in y
    curlBz = Array{Float64, 2}(undef, basis.n_basis, basis.n_quad_points)

    #construct matrices for each element
    for iele in range(1,basis.n_elements)
        basis_orientation = basis.element_basis_orientations[iele]
        basis_orientation += 1 #gmsh to julia offset

        J = diagm(basis.J[:,iele]) #diagonal matrix of quad point jacobians for this element

        #construct basis matrices for integration over quad points (using all basis functions at each quad point)
        for ipoint in 1:basis.n_quad_points
            all_bases_at_point = (basis.jacobians[:,:,ipoint,iele]')\basis.basis[:,:,ipoint,basis_orientation]                
            Bx[:, ipoint] = all_bases_at_point[1, :] #x component of all bases at this one point
            By[:, ipoint] = all_bases_at_point[2, :] #y component of all bases at this one point
            all_curl_bases_at_point = basis.curl_basis[:, :, ipoint, basis_orientation]
            curlBz[:, ipoint] = (1/basis.J[ipoint,iele])*all_curl_bases_at_point[3, :] #just the z component for 2D TE - division by jacobian needed from the math
        end
            
        Mx[:,:,iele] = Bx*(J*W)*Bx'   #mass matrix is integral of basis functions times basis functions       
        My[:,:,iele] = By*(J*W)*By'   #mass matrix is integral of basis functions times basis functions         
            
        S[:, :, iele] = curlBz*(J*W)*curlBz'   #stiffness matrix is integral of (curl psi) dot (curl psi), which only has a z component for 2D TE problems
            
    end

    #create the struct
    element_matrices = ElementMatrices2DTE(Mx, My, S)

    return element_matrices
end


#=------------------------------------
elementMassAndStiffnessMatrices() - 3D
-------------------------------------=#
function elementMassAndStiffnessMatrices(basis::LocalBasis, problem_type::EM3DType)

    #I have a thought to make this tag-based so we can compute them in batches of elements - but can discuss that in the future
    Mx = Array{Float64, 3}(undef,  basis.n_basis, basis.n_basis, basis.n_elements) #nbasis x nbasis x nelements
    My = Array{Float64, 3}(undef,  basis.n_basis, basis.n_basis, basis.n_elements) #nbasis x nbasis x nelements
    Mz = Array{Float64, 3}(undef,  basis.n_basis, basis.n_basis, basis.n_elements) #nbasis x nbasis x nelements
    Sx = Array{Float64, 3}(undef, basis.n_basis, basis.n_basis, basis.n_elements) #nbasis x nbasis x nelements
    Sy = Array{Float64, 3}(undef, basis.n_basis, basis.n_basis, basis.n_elements) #nbasis x nbasis x nelements
    Sz = Array{Float64, 3}(undef, basis.n_basis, basis.n_basis, basis.n_elements) #nbasis x nbasis x nelements
    
    all_bases_at_point = Array{Float64, 2}(undef, 3, basis.n_basis)
    all_curl_bases_at_point = Array{Float64, 2}(undef, 3, basis.n_basis)

    W = diagm(basis.quad_weights)    #element-independent weight matrix
    Bx = Array{Float64, 2}(undef, basis.n_basis, basis.n_quad_points) #derivative of basis in x
    By = Array{Float64, 2}(undef, basis.n_basis, basis.n_quad_points) #derivative of basis in y
    Bz = Array{Float64, 2}(undef, basis.n_basis, basis.n_quad_points) #derivative of basis in y
    curlBx = Array{Float64, 2}(undef, basis.n_basis, basis.n_quad_points)
    curlBy = Array{Float64, 2}(undef, basis.n_basis, basis.n_quad_points)
    curlBz = Array{Float64, 2}(undef, basis.n_basis, basis.n_quad_points)

    #construct matrices for each element
    for iele in range(1,basis.n_elements)
        basis_orientation = basis.element_basis_orientations[iele]
        basis_orientation += 1 #gmsh to julia offset

        J = diagm(basis.J[:,iele]) #diagonal matrix of quad point jacobians for this element

        #construct basis matrices for integration over quad points (using all basis functions at each quad point)
        for ipoint in 1:basis.n_quad_points
            all_bases_at_point = (basis.jacobians[:,:,ipoint,iele]')\basis.basis[:,:,ipoint,basis_orientation]                
            Bx[:, ipoint] = all_bases_at_point[1, :] #x component of all bases at this one point
            By[:, ipoint] = all_bases_at_point[2, :] #y component of all bases at this one point
            Bz[:, ipoint] = all_bases_at_point[3, :] #y component of all bases at this one point
            all_curl_bases_at_point = basis.curl_basis[:, :, ipoint, basis_orientation]
            curlBx[:, ipoint] = (1/basis.J[ipoint,iele])*all_curl_bases_at_point[1, :]
            curlBy[:, ipoint] = (1/basis.J[ipoint,iele])*all_curl_bases_at_point[2, :]
            curlBz[:, ipoint] = (1/basis.J[ipoint,iele])*all_curl_bases_at_point[3, :]
        end
            
        Mx[:,:,iele] = Bx*(J*W)*Bx'   #mass matrix is integral of basis functions times basis functions       
        My[:,:,iele] = By*(J*W)*By'   #mass matrix is integral of basis functions times basis functions         
        Mz[:,:,iele] = Bz*(J*W)*Bz'   #mass matrix is integral of basis functions times basis functions         
        
        Sx[:, :, iele] = curlBx*(J*W)*curlBx'   #stiffness matrix is integral of (curl psi) dot (curl psi), which only has a z component for 2D TE problems
        Sy[:, :, iele] = curlBy*(J*W)*curlBy'   #stiffness matrix is integral of (curl psi) dot (curl psi), which only has a z component for 2D TE problems
        Sz[:, :, iele] = curlBz*(J*W)*curlBz'   #stiffness matrix is integral of (curl psi) dot (curl psi), which only has a z component for 2D TE problems
            
    end

    #create the struct
    element_matrices = ElementMatrices3D(Mx, My, Mz, Sx, Sy, Sz)

    return element_matrices
end

