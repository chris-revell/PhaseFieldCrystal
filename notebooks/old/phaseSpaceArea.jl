using DrWatson
@quickactivate
using UnPack
using CairoMakie
using ColorSchemes
using DifferentialEquations
using JLD2
using LaTeXStrings
using GeometryBasics
using PerceptualColourMaps

folderPath = "PhaseSpace2"

rs            = [0.30,0.35,0.40,0.45,0.50,0.55]
ms            = [0.1,0.2,0.3,0.4,0.5]
ϕ0s           = [-0.420,-0.425,-0.430,-0.435,-0.440,-0.445]

# rs            = [0.5,0.6,0.7,0.8,0.9,1.0]
# ms            = [0.1,0.3,0.5,0.7,0.9,1.1]
# ϕ0s           = [-0.36,-0.38,-0.39,-0.40,-0.42,-0.44]

set_theme!(figure_padding=0,fontsize=14)

runs = [f for f in readdir(datadir("fromCSF",folderPath,"processed")) if f[end-4:end]==".jld2"]

# for m in ms
m=0.1
    filteredRuns = [f for f in runs if occursin("m=$m",f)]

    results = Dict()
    # rResults = Vector{Float64}([])
    # ϕ0Results = Vector{Float64}([])
    # sumResults = Vector{Float64}([])
    for (i,f) in enumerate(filteredRuns)
        data = load(datadir("fromCSF",folderPath,"processed",f))
        @unpack u, params = data
        @unpack ϕ0, r, m, nX, nY, lX, a, δt, tMax = params
        # push!(rResults,r)
        # push!(ϕ0Results,ϕ0)
        # push!(sumResults,sum(u[end][u[end].>0.5]))
        results[(ϕ0,r)] = sum(u[end][u[end].>0.0])./size(u[end])[1]
    end
# display(sumResults)
    fig = Figure(resolution=(600,600))
    ax = Axis(fig[1,1])    
    scatter!(ax,Point2.(first.(collect(results))),markersize=25,color=last.(collect(results)),colormap=(cmap("L2")),colorrange=(minimum(last.(collect(results))),maximum(last.(collect(results)))))
    Colorbar(fig[1,2],limits=(minimum(last.(collect(results))),maximum(last.(collect(results)))),colormap=(cmap("L2")))
    display(fig)
    # save(datadir("fromCSF",folderPath,"plots","sum_m=$m.png"),fig)
# end

    