# This file will generate structs for various parameters used in JoukowskiMain.jl
module generateStructs
export foilParams, vortexParams, sysParams

mutable struct foilParams
    x
    y
    theta
    U
    V
    Omega
    gammac
    rc
    ua
    uf
end

mutable struct vortexParams
    maxNumOfVortices::Int
    timeOfVortexSheddingEvent
    vortexFlag::Int
    mvortex
    mvortex0
end

mutable struct sysParams
    n # Dimension of configuration manifold (3 for SE(2))
    N # Dimension of the entire space given that maxNumOfVortices are permitted in the flow
    Pinitial # total initial momentum in the system
    npv::Int # number of initial vortices in flow
end

end
