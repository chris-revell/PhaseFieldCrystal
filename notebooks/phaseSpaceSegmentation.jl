using DrWatson
@quickactivate
using UnPack
using CairoMakie
using ColorSchemes
using DifferentialEquations
using JLD2
using LaTeXStrings

folderPath = "PhaseSpace"

rs            = [0.5,0.6,0.7,0.8,0.9,1.0]
ms            = [0.1,0.3,0.5,0.7,0.9,1.1]
ϕ0s           = [-0.36,-0.38,-0.39,-0.40,-0.42,-0.44]

set_theme!(figure_padding=0,fontsize=32)

runs = [f for f in readdir(datadir("fromCSF",folderPath)) if f[end-4:end]==".jld2"]

for m in ms
    filteredRuns = [f for f in runs if occursin("m=$m",f)]
    
    fig = Figure(resolution=(6000,6000))
    axes = Dict()
    for (i,f) in enumerate(filteredRuns)
        
        data = load(datadir("fromCSF",folderPath,f))
        @unpack sol, params = data
        @unpack ϕ0, r, m, nX, nY, lX, a, δt, tMax = params
    
        ax = CairoMakie.Axis(fig[i%6+1,i÷6+1],aspect=DataAspect())
        ax.yreversed = true
        heatmap!(ax,transpose(reshape(sol.u[end],(nY,nX))),colorrange=(-1.0, 1.0),colormap=:bwr)
        hidedecorations!(ax)
        hidespines!(ax)
        # ax.title = "phi_0=$(ϕ0), r= $r"
        axes[(r,ϕ0)] = ax
    end

    for (i,r) in enumerate(rs)
        for (j,ϕ₀) in enumerate(ϕ0s)
            if (r,ϕ₀) in keys(axes)
                fig[i,j] = axes[(r,ϕ₀)]
            end
        end
    end

    Label(fig[1,1,Left()],L"r=0.5",textsize=96)
    Label(fig[2,1,Left()],L"r=0.6",textsize=96)
    Label(fig[3,1,Left()],L"r=0.7",textsize=96)
    Label(fig[4,1,Left()],L"r=0.8",textsize=96)
    Label(fig[5,1,Left()],L"r=0.9",textsize=96)
    Label(fig[6,1,Left()],L"r=1.0",textsize=96)

    Label(fig[6,1,Bottom()],L"\phi_0=-0.36",textsize=96)
    Label(fig[6,2,Bottom()],L"\phi_0=-0.38",textsize=96)
    Label(fig[6,3,Bottom()],L"\phi_0=-0.39",textsize=96)
    Label(fig[6,4,Bottom()],L"\phi_0=-0.40",textsize=96)
    Label(fig[6,5,Bottom()],L"\phi_0=-0.42",textsize=96)
    Label(fig[6,6,Bottom()],L"\phi_0=-0.44",textsize=96)

    colsize!(fig.layout,1,Aspect(1.0,1.0))
    colsize!(fig.layout,2,Aspect(1.0,1.0))
    colsize!(fig.layout,3,Aspect(1.0,1.0))
    colsize!(fig.layout,4,Aspect(1.0,1.0))
    colsize!(fig.layout,5,Aspect(1.0,1.0))
    colsize!(fig.layout,6,Aspect(1.0,1.0))
    resize_to_layout!(fig)

    display(fig)
    save(datadir("fromCSF","PhaseSpace","plots","phasespace_m=$m.png"),fig)
end

    