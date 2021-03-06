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
using JLD2

function visualise(sol, freeEnergies, params, path)

    @unpack nX, nY, lX, r, ϕ0, a, δt, tMax = params

    uInternal = Observable(rand(nX,nY))

    fig1 = Figure()
    ga1 = fig1[1,1] = GridLayout()
    ax1 = CairoMakie.Axis(ga1[1,1],aspect=DataAspect())
    heatmap!(ax1,uInternal,colorrange=(-1.0, 1.0),colormap=:bwr)
    #Colorbar(ga1[1, 2],colormap=:viridis,vertical=true)
    hidedecorations!(ax1)
    hidespines!(ax1)
    ax1.title = "t=0.0"
    ax1.yreversed = true

    tSteps = range(1,length(sol.t),step=1)

    record(fig1,"$(path[1:end-5])_u.mp4",tSteps; framerate=10) do i
        ax1.title = "t=$(@sprintf("%.2f", sol.t[i]))"
        uInternal[] = transpose(reshape(sol.u[i],(nY,nX)))
        uInternal[] = uInternal[]
    end

    fig2 = Figure()
    ax2 = CairoMakie.Axis(fig2[1,1])
    lines!(ax2,sol.t,freeEnergies)
    ax2.xlabel = "Time"
    ax2.ylabel = "Free Energy"
    save("$(path[1:end-5])_freeEnergyVsTime.png",fig2)

    fig3 = Figure(figure_padding=0,resolution=(1371,1049))
    ax3 = CairoMakie.Axis(fig3[1,1],aspect=DataAspect())
    ax3.yreversed = true
    heatmap!(ax3,transpose(reshape(sol.u[end],(nY,nX))),colorrange=(-1.0, 1.0),colormap=:bwr)
    hidedecorations!(ax3)
    hidespines!(ax3)
    resize_to_layout!(fig3)
    save("$(path[1:end-5])_finalState.png",fig3)

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
