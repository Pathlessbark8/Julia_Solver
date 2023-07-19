function TECurlGreensFunction(x,y,xp,yp)
    R = sqrt((x-xp)^2+(y-yp)^2)
    jkR = 1im*ω*sqrt(ϵ0*μ0)*R
    t = exp(-jkR)*(jkR+1) / R^3 / (4*π)
    E = -t*SA_F64[(y - yp); -(x - xp)]
end

# normalized by dividing by sqrt(a)
@inline function gaussSin(t, a, ωc)
    return exp(-(a*t)^2) * sin(ωc*t)
end

# normalize by multiplying by ωc
@inline function gaussDeriv(t, ωc)
    return -2*t*ωc*exp(-(ωc*t)^2 / 4)
end

