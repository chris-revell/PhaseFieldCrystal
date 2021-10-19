#
#  Visualise.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 30/06/2021.
#
#
# Function to visualise results as animated gif; second function to import existing results from file

module Visualise

# Import Julia packages
using Plots
using ColorSchemes
using Printf
using JLD2
using UnPack
using DifferentialEquations

function visualise(sol, freeEnergies, params, path)

    @unpack nGrid, lSpace, h, r, ϕ₀, a, δt, tMax = params

    mat1 = zeros(nGrid^2)

    # # Plot phase field
    ENV["GKSwstype"]=100
    ENV["JULIA_GR_PROVIDER"] = "GR"
    anim = @animate for (i,u) in enumerate(sol.u)
        uInternal = reshape(u,(nGrid,nGrid))
        heatmap(uInternal,title="t=$(@sprintf("%.2f", sol.t[i]))",clims=(-1,1),aspect_ratio=:equal,border=:none,show=false,color=:hawaii)
    end
    @info "Saving animated gif"
    gif(anim,"$(path[1:end-5])_u.gif",fps=10)

    plot(sol.t,freeEnergies)
    savefig("$(path[1:end-5])_freeEnergyVsTime.png")

    return nothing

end

# Function to import data needed for visualisation from file
# For example, run visualise(importData("output/folder/data.jld2")...)
function importData(path)

    data = load(path)
    @unpack sol, freeEnergies, params = data

    return sol, freeEnergies, params, path

end

export visualise, importData

end
