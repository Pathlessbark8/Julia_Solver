module Sources

include("AbstractSources.jl")

export AbstractSource

include("EMSources.jl")

export EMAbstractSource, EMTMElectricLineSource, EMTMPlanewaveSource, EMTMGaussianCurrentDensity, evaluateSource

end