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
using DrWatson

function visualise(u, t, ϕ0, r, m, nX, nY, lX, a, δt, tMax, subFolder, fileName)
    
    fig1 = Figure(figure_padding=0,resolution=(1000,1000),fontsize=64)
    ax1 = CairoMakie.Axis(fig1[1,1],aspect=DataAspect())
    uInternal1 = Observable(zeros(nX,nY))
    uInternal2 = Observable(zeros(nX,nY))
    heatmap!(ax1,uInternal1,colorrange=(-1.0, 1.0),colormap=(:bwr,0.5))
    heatmap!(ax1,uInternal2,colorrange=(-1.0, 1.0),colormap=(:bwr,0.5))
    hidedecorations!(ax1)
    hidespines!(ax1)
    ax1.title = "t=0.0"
    ax1.yreversed = true
    resize_to_layout!(fig1)
    tSteps = range(1,length(t),step=1)
    record(fig1,"$subFolder/$(fileName)_u.mp4",tSteps; framerate=10) do i
        display(i)
        ax1.title = "t=$(@sprintf("%.2f", t[i]))"
        uInternal1[] = transpose(reshape(u[i][1:nX*nY],(nY,nX)))
        uInternal1[] = uInternal1[]
        uInternal2[] = transpose(reshape(u[i][1+nX*nY:end],(nY,nX)))
        uInternal2[] = uInternal2[]
        save("$subFolder/$(fileName)$i.png",fig1)
    end

    # if freeEnergyFlag==1
    #     fig2 = Figure(figure_padding=0)
    #     ax2 = CairoMakie.Axis(fig2[1,1])
    #     lines!(ax2,t,freeEnergies)
    #     ax2.xlabel = "Time"
    #     ax2.ylabel = "Free Energy"
    #     safesave("$subFolder/$(fileName)_freeEnergyVsTime.png",fig2)
    # end

    # fig3 = Figure(figure_padding=0,resolution=(1000,1000))
    # ax3 = CairoMakie.Axis(fig3[1,1],aspect=DataAspect())
    # ax3.yreversed = true
    # heatmap!(ax3,transpose(reshape(u[end][1:nX*nY],(nY,nX))),colorrange=(-1.0, 1.0),colormap=(:bwr,0.5))
    # heatmap!(ax3,transpose(reshape(u[end][1+nX*nY:end],(nY,nX))),colorrange=(-1.0, 1.0),colormap=(:bwr,0.5))
    # hidedecorations!(ax3)
    # hidespines!(ax3)
    # resize_to_layout!(fig3)
    # # display("$subFolder/$(fileName[1:end-5])_finalState.png")
    # save("$subFolder/$(fileName)_finalState.png",fig3)

    return nothing

end

# Function to import data needed for visualisation from file
# For example, run visualise(importData("output/folder/data.jld2")...)
function importData(path)

    data = load(path)
    @unpack u, t, ϕ0, r, m, λ, nX, nY, lX, a, δt, tMax = data

    pathparts = splitpath(path)

    return u, t, ϕ0, r, m, nX, nY, lX, a, δt, tMax, joinpath(pathparts[1:end-1]), pathparts[end][1:end-5], 0

end

export visualise, importData

end
