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
using CairoMakie
using ColorSchemes
using Printf
using UnPack
using DifferentialEquations

function visualise(sol, freeEnergies, params, path)

    @unpack nX, nY, lSpace, r, ϕ0, a, δt, tMax = params
    println("test")
    uInternal = Observable(rand(nX,nY))

    fig1 = Figure()
    ga1 = fig1[1,1] = GridLayout()
    ax1 = Axis(ga1[1,1],aspect=DataAspect())
    heatmap!(ax1,uInternal,colorrange=(-1.0, 1.0),colormap=:viridis)
    Colorbar(ga1[1, 2],colormap=:viridis,vertical=true)
    hidedecorations!(ax1)
    ax1.title = "t=0.0"
    ax1.yreversed = true

    tSteps = range(1,length(sol.t),step=1)

    record(fig1,"$(path[1:end-5])_u.gif",tSteps; framerate=10) do i
        ax1.title = "t=$(sol.t[i])"
        uInternal[] = transpose(reshape(sol.u[i],(nY,nX)))
        uInternal[] = uInternal[]
    end

    fig2 = Figure()
    ax2 = Axis(fig2[1,1])
    lines!(ax2,sol.t,freeEnergies)
    ax2.xlabel = "Time"
    ax2.ylabel = "Free Energy"
    save("$(path[1:end-5])_freeEnergyVsTime.png",fig2)

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
