using SpecialFunctions

abstract type AbstractEMSource <: AbstractSource end

# we need to settle on a time/frequency domain break down at some point soon.

struct EMTMElectricLineSource <: AbstractEMSource
    x::Float64      # x-position
    y::Float64      # y-position
    Jz::ComplexF64  # complex amplitude
    ε::ComplexF64   # homogeneous permittivity of medium 
    μ::ComplexF64   # homogeneous permeability of medium
end

struct EMTMPlanewaveSource <: AbstractEMSource
    khatx::Float64  # x-direction
    khaty::Float64  # y-direction
    A::ComplexF64   # complex amplitude
    ε::ComplexF64   # homogeneous permittivity of medium 
    μ::ComplexF64   # homogeneous permeability of medium
end

struct EM3DDipole <: AbstractEMSource 
    
end

struct EMTMGaussianCurrentDensity <: AbstractEMSource
    μx::Float64     # x-directed mean
    μy::Float64     # y-directed mean
    σx::Float64     # x-directed standard deviation
    σy::Float64     # y-directed standard deviation
    Jz::Float64     # complex amplitude scaling
end

function evaluateFieldsFromSource(source::EMTMElectricLineSource, points::Array{Float64,2}, frequency::Float64)
    
    ω = 2*pi*frequency
    
    rx = points[1,:] .- source.x
    ry = points[2,:] .- source.y
    r = sqrt.(rx .* rx + ry .* ry)
    
    k = ω*sqrt(source.ε * source.μ)
    η = sqrt(source.μ / source.ε)
    
    Ez = -source.Jz*(1im/4.0)*besselh.(0, 2, k*r)
    Hx = -source.Jz*(1.0/(4.0*η))*(ry./r).*besselh.(1, 2, k*r)
    Hy = source.Jz*(1.0/(4.0*η))*(rx./r).*besselh.(1, 2, k*r)

    return Ez, Hx, Hy
end

function evaluateFieldsFromSource(source::EMTMPlanewaveSource, points::Array{Float64, 2}, frequency::Float64)
    ω = 2*pi*frequency
    k = ω*sqrt(source.ε * source.μ)
    kx = source.khatx*k
    ky = source.khaty*k
    kdotr = kx .* points[1,:] + ky .* points[2, :]

    Ez = source.exp(-1im*kdotr)
    Hx = 0 .* Ez  # will need to update these appropriately at some point
    Hy = 0 .* Ez

    return Ez, Hx, Hy
end

function evaluateSource(source::EMTMGaussianCurrentDensity, points::Array{Float64,2}, frequency::Float64)

    Jz = (1 / (sqrt(2*π)*source.σx) * 1 / (sqrt(2*π)*source.σy)) .* exp.(-(points[1,:] .- source.μx).^2 ./ (2*source.σx^2)) .* exp.(-(points[2,:] .- source.μy).^2 ./ (2*source.σy^2)) 
    #Mx = 0 .* Jz
    #My = 0 .* Jz
    return Jz #, Mx, My
end
