using JLD2
using DrWatson
using LinearAlgebra
using SparseArrays
using FromFile: @from
using CairoMakie
using DelaunayTriangulation
using StaticArrays
using UnPack

# # Import local files
@from "$(srcdir("CreateLaplacian.jl"))" using CreateLaplacian
@from "$(srcdir("CreateDivAlphaGrad.jl"))" using CreateDivAlphaGrad

runs = [f for f in readdir(datadir("exp_pro", "emCentroidsInteractive")) if occursin(".jld2", f)]

data = load(datadir("exp_pro", "emCentroidsInteractive", runs[1]))

@unpack centroidLocations = data

ptsArray = zeros(Float64, (2, length(centroidLocations)))
for (i, point) in enumerate(centroidLocations)
    ptsArray[1, i] = point[1]
    ptsArray[2, i] = point[2]
end

triangulation_unconstrained = triangulate(ptsArray)
tessellation_constrained = voronoi(triangulation_unconstrained, true)

#Exclude points outside constraining boundary
usableVertices = Int64[]
for a in values(tessellation_constrained.polygons)
    push!(usableVertices,a...)
end
sort!(unique!(usableVertices))
outerVertices = setdiff(collect(1:num_polygon_vertices(tessellation_constrained)),usableVertices)

# Map vertex indices in tessellation to vertex indices in incidence matrices (after excluding outer vertices)
vertexIndexingMap = Dict(usableVertices.=>collect(1:length(usableVertices)))

R = SVector.(tessellation_constrained.polygon_points[usableVertices])

# Find pairs of vertices connected by edges in tessellation 
# Use incidence matrix indexing for vertices, and exclude outer vertices 
pairs = [(vertexIndexingMap[p[1]],vertexIndexingMap[p[2]]) for p in keys(tessellation_constrained.adjacent.adjacent) if p[1]∈usableVertices && p[2]∈usableVertices]
# Ensure lowest index is first in tuple, and remove duplicates 
orderedPairs = unique([(min(p...), max(p...)) for p in pairs])

nVerts = length(R)
nEdges = length(orderedPairs)
nCells = length(centroidLocations)

# Construct A matrix mapping tessellation edges to tessellation vertices 
A = spzeros(Int64, nEdges, nVerts)
for (edgeIndex, vertices) in enumerate(orderedPairs)
    A[edgeIndex, vertices[1]] = 1
    A[edgeIndex, vertices[2]] = -1
end


# NB get_polygon(tessellation_constrained,x) or tessellation_constrained.polygons[x] return indices of vertices around cell x ordered anti-clockwise, with first and last element the same

# Construct B matrix mapping voronoi cell around each fibril to surrounding edges between vertices in tessellation
# NB assume ϵᵢ is a clockwise rotation so cell orientation is into page. 
B = spzeros(Int64, nCells, nEdges)
for c = 1:nCells
    for i = 2:length(tessellation_constrained.polygons[c])
        vertexLeading = vertexIndexingMap[tessellation_constrained.polygons[c][i-1]]  # Leading with respect to *clockwise* direction around cell
        vertexTrailing = vertexIndexingMap[tessellation_constrained.polygons[c][i]]
        # Find index of edge connecting these vertices 
        edge = (findall(x -> x != 0, A[:, vertexLeading])∩findall(x -> x != 0, A[:, vertexTrailing]))[1]
        if A[edge, vertexLeading] > 0
            B[c, edge] = 1
        else
            B[c, edge] = -1
        end
    end
end