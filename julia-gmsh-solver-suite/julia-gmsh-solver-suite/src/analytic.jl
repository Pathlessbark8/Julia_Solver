function addPointSource(x, y, z, basisOrder; curl=true)
    tag, _, _, u, v, w = gmsh.model.mesh.getElementByCoordinates(x, y, z, 2)
    localCoord = Float64[u, v, w]
    if curl
        basisName = "CurlHcurlLegendre$basisOrder"
    else
        basisName = "HcurlLegendre$basisOrder"
    end

    # integrate line source
    type, _, _, _ = gmsh.model.mesh.getElement(tag)
    _, localFuncs, numOrientations = gmsh.model.mesh.getBasisFunctions(type, localCoord, basisName)
    numFunctions = length(localFuncs) ÷ numOrientations ÷ 3
    localFuncs = reshape(reinterpret(SVector{3,Float64}, localFuncs), (numFunctions, numOrientations))

    jac, det, _ = gmsh.model.mesh.getJacobian(tag, localCoord)
    jac = reinterpret(SMatrix{3, 3, Float64, 9}, jac)

    orient = 1 + gmsh.model.mesh.getBasisFunctionsOrientationForElement(tag, basisName)
    
    if curl
        globEvaluation = [(localFuncs[i,orient] / det[1])[3] for i=1:numFunctions]
    else
        globEvaluation = [(jac[1]' \ localFuncs[i,orient])[1] for i=1:numFunctions]
    end
    return globEvaluation, tag
end

function E(x,y,xp,yp)
    R = sqrt((x-xp)^2+(y-yp)^2)
    jkR = 1im*ω*sqrt(ϵ0*μ0)*R
    t = exp(-jkR)*(jkR+1) / R^3 / (4*π)
    E = -t*SA_F64[(y - yp); -(x - xp)]
end
