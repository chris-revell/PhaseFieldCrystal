using UnPack
using JLD2
using DrWatson
using DifferentialEquations
using LinearAlgebra
using SparseArrays
using FromFile: @from
using CairoMakie
using Makie.GeometryBasics
using VoronoiCells
using ConcaveHull
using GR: delaunay

# Import local files
@from "$(srcdir("CreateLaplacian.jl"))" using CreateLaplacian
@from "$(srcdir("CreateDivAlphaGrad.jl"))" using CreateDivAlphaGrad

function neighbourColours(x)
    if x == 6
        return (:white, 0.0)
    elseif x == 5
        return (:red, 1.0)
    elseif x == 7
        return (:blue, 1.0)
    else
        return (:grey, 1.0)
    end
end

runs = [f for f in readdir(datadir("exp_pro","emCentroidsInteractive")) if occursin(".jld2",f)]

data = load(datadir("exp_pro","emCentroidsInteractive",runs[1]))

@unpack centroidLocations = data

xs = [x[1] for x in centroidLocations]
ys = [x[2] for x in centroidLocations]
n, tri = delaunay(xs, ys)


scalingFactor = maximum(abs.([xs ys])) / (1 - 3eps(Float64))
shiftedCentroidLocations = centroidLocations ./ scalingFactor
# shiftedCentroidLocations .+= Point2(0, 1)

# Count neighbours of each centroid in the triangulation 
nNeighbours = [length(findall(x -> x == i, tri)) for i = 1:length(shiftedCentroidLocations)]

# Concave hull to identify boundary fibrils 
hull = concave_hull(centroidLocations, 1)
# Indices of fibrils within the hull 
hullInds = [findall(x -> Point2(x...) == v, centroidLocations)[1] for v in hull.vertices]
hullPoints = [Point2(Float64(h[1]),Float64(h[2])) for h in hull.vertices]

# Voronoi tessellation of centroid positions within (0,0) (1,1) box
tess = voronoicells(shiftedCentroidLocations, Rectangle(Point2(0, 0), Point2(1, 1)))
tessAreas = voronoiarea(tess)



fig = Figure(resolution=(500,500))
ax = Axis(fig[1,1],aspect=DataAspect())
for (j, c) in enumerate(tess.Cells)
    # if j ∉ hullInds
        vertices = [v .* scalingFactor for v in c]
        poly!(ax, vertices, color=(j ∉ hullInds ? neighbourColours(nNeighbours[j]) : :black), strokecolor=(:black, 1.0), strokewidth=1.0)
    # end
end
for i=1:n
    poly!(ax,centroidLocations[tri[i,:]],color=(:white,0.0),strokecolor=(:orange,1.0),strokewidth=1.0)
end
scatter!(ax,centroidLocations)

xlims!(ax,0.9*minimum(first.(centroidLocations)),1.05*maximum(first.(centroidLocations)))
ylims!(ax,0.9*minimum(last.(centroidLocations)),1.05*maximum(last.(centroidLocations)))
hidedecorations!(ax)
hidespines!(ax)

limits = ax.finallimits[]
boundingRectangle = [
    Point2(0.0,0.0), 
    Point2(700.0, 0.0), 
    Point2(700.0,700.0), 
    Point2(0.0,700.0)
]

p = Polygon(
    boundingRectangle,
    [hullPoints]
)
poly!(ax,p,color=:blue)


display(fig)