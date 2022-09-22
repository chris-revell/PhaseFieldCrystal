using DrWatson
@quickactivate
using UnPack
using CairoMakie
using ColorSchemes
using DifferentialEquations
using JLD2
using LaTeXStrings

folderPath = "PhaseSpace"

# rs            = [0.30,0.35,0.40,0.45,0.50,0.55]
# ms            = [0.1,0.2,0.3,0.4,0.5]
# ϕ0s           = [-0.420,-0.425,-0.430,-0.435,-0.440,-0.445]

rs            = [0.5,0.6,0.7,0.8,0.9,1.0]
ms            = [0.1,0.3,0.5,0.7,0.9,1.1]
ϕ0s           = [-0.36,-0.38,-0.39,-0.40,-0.42,-0.44]

set_theme!(figure_padding=0,fontsize=32)

runs = [f for f in readdir(datadir("fromCSF",folderPath,"processed")) if f[end-4:end]==".jld2"]

for m in ms
    filteredRuns = [f for f in runs if occursin("m=$m",f)]
    
    fig = Figure(resolution=(6000,6000))
    axesDict = Dict()
    for (i,f) in enumerate(filteredRuns)
        
        data = load(datadir("fromCSF",folderPath,f))
        @unpack u, params = data
        @unpack ϕ0, r, m, nX, nY, lX, a, δt, tMax = params
    
        ax = CairoMakie.Axis(fig[i%6+1,i÷6+1],aspect=DataAspect())
        ax.yreversed = true
        heatmap!(ax,transpose(reshape(u[end],(nY,nX))),colorrange=(-1.0, 1.0),colormap=:bwr)
        hidedecorations!(ax)
        hidespines!(ax)
        # ax.title = "phi_0=$(ϕ0), r= $r"
        axesDict[(r,ϕ0)] = ax
    end

    for (i,r) in enumerate(rs)
        for (j,ϕ₀) in enumerate(ϕ0s)
            if (r,ϕ₀) in keys(axes)
                fig[7-i,7-j][1,1] = axesDict[(r,ϕ₀)]
                Colorbar(fig[7-i,7-j][1,2], limits=(-1,1),colormap=:bwr)#, size = 25)
                Label(fig[7-i,7-j,Bottom()],L"r=%$r, \phi_0=%$ϕ₀",textsize=64)
            end
        end
    end

    colsize!(fig.layout,1,Aspect(1.0,1.0))
    colsize!(fig.layout,2,Aspect(1.0,1.0))
    colsize!(fig.layout,3,Aspect(1.0,1.0))
    colsize!(fig.layout,4,Aspect(1.0,1.0))
    colsize!(fig.layout,5,Aspect(1.0,1.0))
    colsize!(fig.layout,6,Aspect(1.0,1.0))
    resize_to_layout!(fig)

    display(fig)
    save(datadir("fromCSF",folderPath,"plots","phasespace_m=$m.png"),fig)
end

    