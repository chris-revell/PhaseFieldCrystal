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
using Colors

function filterNonFibril(x)
    if x<0.0
        return x
    else
        return 0.0
    end
end

function visualise(u, t, nX, nY, subFolder, fileName)

    fig1 = Figure(figure_padding=0, resolution=(1000, 1000), fontsize=64)
    ax1 = CairoMakie.Axis(fig1[1, 1], aspect=DataAspect())
    uTmp = zeros(nX, nY)
    uInternal1 = Observable(zeros(nX, nY))
    uInternal2 = Observable(zeros(nX, nY))
    heatmap!(ax1, uInternal1, colorrange=(-1.0, 1.0), colormap=(:bwr, 1.0))#)ColorScheme([RGBA{Float64}(0.0, 0.0, 0.0, i) for i in 0.0:0.001:1.0]))
    # heatmap!(ax1, uInternal2, colorrange=(-1.0, 1.0), colormap=(:bwr, 0.5))
    hidedecorations!(ax1)
    hidespines!(ax1)
    ax1.title = "t=0.0"
    ax1.yreversed = true
    resize_to_layout!(fig1)
    tSteps = range(1, length(t), step=1)
    mov = VideoStream(fig1, framerate=10)

    for i=1:length(t)
        display(i)
        ax1.title = "t=$(@sprintf("%.2f", t[i]))"
        # uInternal1[] = min.(transpose(reshape(u[i][1:nX*nY], (nY, nX))),transpose(reshape(u[i][1+nX*nY:end], (nY, nX))))
        uTmp .= filterNonFibril.(transpose(reshape(u[i][1:nX*nY], (nY, nX)))) .- filterNonFibril.(transpose(reshape(u[i][1+nX*nY:end], (nY, nX))))
        uInternal1[] .= uTmp # transpose(reshape(u[i][1:nX*nY], (nY, nX)))
        uInternal1[] = uInternal1[]
        # uTmp = 
        # uInternal2[] = -1.0.*filterNonFibril.(uTmp)
        # uInternal2[] = uInternal2[]

        recordframe!(mov)
        # save("$subFolder/$(fileName)$i.png", fig1)
    end
    safesave(datadir("sims", subFolder, "$fileName.mp4"), mov)

    # if freeEnergyFlag==1
    #     fig2 = Figure(figure_padding=0)
    #     ax2 = CairoMakie.Axis(fig2[1,1])
    #     lines!(ax2,t,freeEnergies)
    #     ax2.xlabel = "Time"
    #     ax2.ylabel = "Free Energy"
    #     safesave("$subFolder/$(fileName)_freeEnergyVsTime.png",fig2)
    # end

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
