using DrWatson
@quickactivate
using UnPack
using CairoMakie; set_theme!(figure_padding=0,fontsize=32)
using ColorSchemes
using DifferentialEquations
using JLD2
using LaTeXStrings
using DataFrames
using PerceptualColourMaps

folderPath = "PhaseSpace"

results = collect_results!(datadir("fromCSF",folderPath,"processed"); subfolders = true)

fig = Figure(resolution=(6000,6000))
axes = Dict()

df = DataFrame(t=Float64[], ϕ0=Float64[], r=Float64[], m=Float64[], nX=Int64[], nY=Int64[], lX=Float64[], a=Float64[], δt=Float64[], tMax=Float64[])
us = []
for result in eachrow(results)
    params = result[:params]
    @unpack ϕ0, r, m, nX, nY, lX, a, δt, tMax = params
    push!(df, Dict(:t=>tMax, :ϕ0=>ϕ0, :r=>r, :m=>m, :nX=>nX, :nY=>nY, :lX=>lX, :a=>a, :δt=>δt, :tMax=>tMax))
    push!(us,result[:u])
end


for (i,result) in enumerate(eachrow(subset(df, :m => m -> m.== 0.1)))
    ax = CairoMakie.Axis(fig[i%6+1,i÷6+1],aspect=DataAspect())
    ax.yreversed = true
    heatmap!(ax,transpose(reshape(us[i][end],(nY,nX))),colorrange=(-1.0, 1.0),colormap=cmap("D1"))
    hidedecorations!(ax)
    hidespines!(ax)    
    axes[(result[:r],result[:ϕ0])] = ax
end 

sortedrs  = unique!(sort(first.(keys(axes))))
sortedϕ0s = unique!(sort(last.(keys(axes))))

for (i,r) in enumerate(sortedrs)
    for (j,ϕ0) in enumerate(sortedϕ0s)
        fig[length(sortedrs)-i,length(sortedϕ0s)-j][1,1] = axes[(r,ϕ0)]
        # Colorbar(fig[length(sortedrs)-i,length(sortedϕ0s)-j][1,2], limits=(-1,1),colormap=:bwr)#, size = 25)
        Label(fig[length(sortedrs)-i,length(sortedϕ0s)-j,Bottom()],L"r=%$r, \phi_0=%$ϕ0",textsize=64)
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


