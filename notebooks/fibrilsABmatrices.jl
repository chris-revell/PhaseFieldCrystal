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

# boundingRectangle = [
#     Point2(0.0,0.0), 
#     Point2(700.0, 0.0), 
#     Point2(700.0,700.0), 
#     Point2(0.0,700.0)
# ]
# p = Polygon(
#     boundingRectangle,
#     [hullPoints]
# )
# poly!(ax,p,color=:blue)

# Find Delaunay triangulation and Voronoi tessellation
pointsArray = zeros(Float64,(2,length(centroidLocations)))
for (i,point) in enumerate(centroidLocations)
    pointsArray[1,i] = point[1]
    pointsArray[2,i] = point[2]
end
triangulation = triangulate(pointsArray)
tessellation = voronoi(triangulation, true)

# Convert vertices of Voronoi tessellation to Point2 type
vertexPoints = Point2.(tessellation.polygon_points)


# Find all edges in Delaunay triangulation
# NB get_edges(triangulation) returns the number of edges including ghost get_edges
# Exclude ghost edges by removing pairs with negative vertex indices
# Use collect to convert set to vector
pairsInTriangulation = collect(filter(x->min(x...)>0,get_edges(triangulation)))
# for p in pairsInTriangulation
#     lines!(ax,centroidLocations[[p...]],color=:black)
# end

edgesInTessellation = Tuple{Int64, Int64}[]
for p in pairsInTriangulation
    # Find tessellation vertices shared by points at both ends of pair p in triangulation
    # Note that the list of pairs does not include ghost edges in the triangulation so there should be no ghost points here...
    tessellationVertices = tessellation.polygons[p[1]]∩tessellation.polygons[p[2]]
    push!(edgesInTessellation,Tuple(tessellationVertices))
end
# Add convex hull as edges 
# Indices in convex hull object are correctly ordered anti-clockwise already, 
# and start and end with the same index
# Note, need to relabel indices for use in A, B matrices, 
# so their indices are all > max index of vertices in tessellation 
for i=2:length(triangulation.convex_hull.indices)
    edgePair = (num_polygon_vertices(tessellation)+i,num_polygon_vertices(tessellation)+i-1)
    push!(edgesInTessellation,edgePair)
end

# R matrix of vertices for A, B matrices includes vertices in tessellation (excluding ghost vertices) 
# and vertices in convex hull of triangulation
R = tessellation.polygon_points
# for i=1:length(triangulation.convex_hull.indices)-1
#     push!(R,Tuple(triangulation.convex_hull.points[:,i]))
# end

# Construct A matrix mapping tessellation edges to tessellation vertices 
A = spzeros(length(edgesInTessellation),num_polygon_vertices(tessellation))
for (edgeIndex,vertices) in enumerate(edgesInTessellation)
    A[edgeIndex,vertices[1]] = 1
    A[edgeIndex,vertices[2]] = -1
end

# B = spzeros(length(centroidLocations),length(pairs))
# for c=1:length(centroidLocations)
#     verticesOfCell = tessellation.polygons[c]
#     for e=2:length(verticesOfCell)
#         edge = findall(x->x==(e,e-1)||x==(e-1,e), pairs)[1]
#         if A[edge,verticesOfCell[e]] == 1
#             B[c,edge] = 1
#         else
#             B[c,edge] = -1
#         end
#     end
# end


fig = Figure(resolution=(750,750)); ax=Axis(fig[1,1],aspect=DataAspect())
xlims!(ax,0.9*minimum(first.(centroidLocations)),1.05*maximum(first.(centroidLocations)))
ylims!(ax,0.9*minimum(last.(centroidLocations)),1.05*maximum(last.(centroidLocations)))
hidedecorations!(ax)
hidespines!(ax)
# scatter!(ax,first.(R),last.(R))
# scatter!(ax,centroidLocations, color=:red)
voronoiplot!(ax, tessellation, strokecolor=:red, strokewidth=0.2, show_generators=false)
triplot!(ax, triangulation, strokewidth=0.2, strokecolor=(:black, 0.4), show_convex_hull=true)
display(fig)