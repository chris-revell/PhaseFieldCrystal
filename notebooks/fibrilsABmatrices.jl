using UnPack
using JLD2
using DrWatson
using DifferentialEquations
using LinearAlgebra
using SparseArrays
using FromFile: @from
using CairoMakie
using Makie.GeometryBasics
using DelaunayTriangulation
using ConcaveHull

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

ptsArray = zeros(Float64,(2,length(centroidLocations)))
for (i,point) in enumerate(centroidLocations)
    ptsArray[1,i] = point[1]
    ptsArray[2,i] = point[2]
end

triangulation_unconstrained = triangulate(ptsArray)
# tessellation_unconstrained = voronoi(triangulation_unconstrained)
tessellation_constrained = voronoi(triangulation_unconstrained,true)

R = [get_polygon_point(tessellation_constrained,i) for i=1:num_polygon_vertices(tessellation_constrained)]


orderedPairs = unique([(min(p...),max(p...)) for p in keys(tessellation_constrained.adjacent.adjacent)])

nVerts = length(R)
nEdges = length(orderedPairs)
nCells = length(centroidLocations)

# Construct A matrix mapping tessellation edges to tessellation vertices 
A = spzeros(nEdges,nVerts)
for (edgeIndex,vertices) in enumerate(orderedPairs)
    A[edgeIndex,vertices[1]] = 1
    A[edgeIndex,vertices[2]] = -1
end


# NB get_polygon(tessellation_constrained,x) or tessellation_constrained.polygons[x] return indices of vertices around cell x ordered anti-clockwise, with first and last element the same

# Construct B matrix mapping voronoi cell around each fibril to surrounding edges between vertices in tessellation
B = spzeros(nCells,nEdges)
for c=1:nCells
    for i=2:length(tessellation_constrained.polygons[c]) # iterating around polygon anti-clockwise 
        vertexLeading = tessellation_constrained.polygons[c][i]  # Leading with respect to anticlockwise direction around cell
        vertexTrailing = tessellation_constrained.polygons[c][i-1]                
        # Find index of edge connecting these vertices 
        edge = (findall(x->x!=0,A[:,vertexLeading])∩findall(x->x!=0,A[:,vertexTrailing]))[1]        
        if A[edge,vertexLeading] > 0
            B[c,edge] = 1
        else
            B[c,edge] = -1
        end
    end
end

fig, ax, sc = voronoiplot(tessellation_constrained, strokecolor=:red, markersize=9)
# fig = Figure(resolution=(750,750)); ax=Axis(fig[1,1],aspect=DataAspect())
# triplot!(ax, triangulation_unconstrained, strokewidth=1, strokecolor=(:black, 0.4))
# annotations!(ax,string.(collect(1:length(centroidLocations))),centroidLocations,color=:red)
display(fig)