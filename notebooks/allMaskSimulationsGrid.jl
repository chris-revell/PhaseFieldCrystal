using DrWatson; @quickactivate
using Images
using CairoMakie
using Colors
using GeometryBasics
using LinearAlgebra
using FromFile
using FileIO
using ConcaveHull
using VoronoiCells
using GR: delaunay
using CSV
using DataFrames

@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions

# Collate results as a dataframe 
results = collect_results!(datadir("fromCSF","allMasks"); subfolders = true)

runs = [f for f in readdir(datadir("exp_pro","masks","ok")) if f[end-3:end]==".png"]

fig = Figure(resolution=(6000,6000),fontsize=64)

axes = Dict()
sizes = Dict()

for (i,r) in enumerate(runs)

    subsetResults = subset(results, :path => m -> occursin.(r[1:end-4],m))

    maskIn = load(datadir("exp_pro","masksCompressed",r[1:end-4],"$(r[1:end-4]).png"))
    maskImage = fill(RGBA(1,1,1,1),size(maskIn))
    for i=1:size(maskIn)[1]
        for j=1:size(maskIn)[2]
            if maskIn[i,j] < 0.5
                maskImage[i,j] = RGBA(0.0,0.0,0.0,1.0)
            else
                maskImage[i,j] = RGBA(0.0,0.0,0.0,0.0)
            end
        end
    end    

    ax = CairoMakie.Axis(fig[(i-1)%6+1,(i-1)÷6+1],aspect=DataAspect())
    uMat = reshape(subsetResults[1,:u][end],(subsetResults[1,:nY],subsetResults[1,:nX]))
    heatmap!(ax,rotr90(uMat),colorrange=(-1.0, 1.0),colormap=:bwr)
    image!(ax,rotr90(maskImage))
    hidedecorations!(ax)
    hidespines!(ax)
    Label(fig[(i-1)%6+1,(i-1)÷6+1, BottomLeft()], "$i", valign = :bottom, font = "TeX Gyre Heros Bold", padding = (0, 10, 10, 0))
    axes[r] = ax
    sizes[r] = size(maskImage)
end 

xMax = maximum(last.(values(sizes)))
yMax = maximum(first.(values(sizes)))

for r in runs 
    xlims!(axes[r],(0,xMax))
    ylims!(axes[r],(0,yMax))
end

for i=1:6
    colsize!(fig.layout, i, Aspect(1, 1.0))
end

# colgap!(fig.layout, 1, -350)
# colgap!(fig.layout, 2, -300)
# colgap!(fig.layout, 3, -300)
# colgap!(fig.layout, 4, -250)
# colgap!(fig.layout, 5, -50)

# rowgap!(fig.layout, 1, -200)
# rowgap!(fig.layout, 2, -200)
# rowgap!(fig.layout, 3, -300)


colgap!(fig.layout, 1, -500)
colgap!(fig.layout, 2, -400)
colgap!(fig.layout, 3, -400)
colgap!(fig.layout, 4, -300)
colgap!(fig.layout, 5, -50)

rowgap!(fig.layout, 1, -200)
rowgap!(fig.layout, 2, -300)
rowgap!(fig.layout, 3, -400)
rowgap!(fig.layout, 4, -100)

resize_to_layout!(fig)
display(fig)

save(datadir("fromCSF","allMasks","grid.png"),fig)



