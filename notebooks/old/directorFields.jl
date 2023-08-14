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


save("test.png",fig)

# Create A matrix mapping voronoi cells to voronoi edges

using DelaunayTriangulation

ptsArray = zeros(Float64,(2,length(centroidLocations)))
for (i,point) in enumerate(centroidLocations)
    ptsArray[1,i] = point[1]
    ptsArray[2,i] = point[2]
end
tri2 = triangulate(ptsArray)
vorn = voronoi(tri2)

vertexPoints = Point2.(vorn.polygon_points)

fig2 = Figure(resolution=(750,750)); ax2=Axis(fig2[1,1],aspect=DataAspect())
xlims!(ax2,0.9*minimum(first.(centroidLocations)),1.05*maximum(first.(centroidLocations)))
ylims!(ax2,0.9*minimum(last.(centroidLocations)),1.05*maximum(last.(centroidLocations)))
hidedecorations!(ax2)
hidespines!(ax2)


scatter!(ax2, centroidLocations,color=:red)
annotations!(ax2, string.(collect(1:length(centroidLocations))), centroidLocations, color=:red)
scatter!(ax2, vertexPoints, color=:blue)
annotations!(ax2, string.(collect(1:length(vertexPoints))), vertexPoints, color=:blue)

display(fig2)


fig3 = Figure(resolution=(750,750)); ax3=Axis(fig3[1,1],aspect=DataAspect())
xlims!(ax3,0.9*minimum(first.(centroidLocations)),1.05*maximum(first.(centroidLocations)))
ylims!(ax3,0.9*minimum(last.(centroidLocations)),1.05*maximum(last.(centroidLocations)))
hidedecorations!(ax3)
hidespines!(ax3)

pairsFull = collect(keys(tri2.adjacent.adjacent))
pairs = unique([(min(p...),max(p...)) for p in pairs])
for p in get_edges(tri2) #pairs
    if min(p...) > 0
        lines!(ax3,centroidLocations[[p...]],color=:black)
    else
        lines!(ax3,[centroidLocations[max(p...)],Point2(0,0)],color=:black)
    end
    # @show a = [centroidLocations[p[1]]], centroidLocations[p[2]]
end

display(fig3)


using SparseArrays
A = spzeros(length(pairs),length(vorn.polygon_points))
for p in pairs 
    A[p[1],p[2]] = 1
    A[p[2],p[1]] = -1
end

B = spzeros(length(centroidLocations),length(pairs))
for c=1:length(centroidLocations)
    verticesOfCell = vorn.polygons[c]
    for e=2:length(verticesOfCell)
        edge = findall(x->x==(e,e-1)||x==(e-1,e), pairs)[1]
        if A[edge,verticesOfCell[e]] == 1
            B[c,edge] = 1
        else
            B[c,edge] = -1
        end
    end
end