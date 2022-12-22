using DataFrames
using CSV
using CairoMakie
using GeometryBasics
using LinearAlgebra
using FileIO
using ConcaveHull
using GR: delaunay
using Statistics
using StatsBase
using DelimitedFiles
using Dates

runs = [r for r in Vector(readdlm(datadir("exp_pro", "filesToUse.txt"))[:, 1]) if !(occursin("mp13ko", r) || occursin("18tailT_4800X_HUI_0007_0", r) || occursin("18tailT_4800X_HUI_0008_0", r))]

# spacingData = DataFrame(CSV.File(datadir("exp_pro","emCentroidMeasurements","spacingData.csv")))
# defectCounts = DataFrame(CSV.File(datadir("exp_pro","emCentroidMeasurements","defectCounts.csv")))

lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro", "lengthMeasurements", "lengthMeasurements.csv")))
croppedLX = DataFrame(CSV.File(datadir("exp_pro", "lengthMeasurements", "croppedLX.csv")))

pairLengthsDict = Dict()
nNeighboursDict = Dict()

voronoiSizeThresh = 2.0

function bothPeripheralTest(label1,label2,hullInds,tessAreas,areaThresh)
    if label1 ∈ hullInds || tessAreas[label1] > areaThresh
        if label2 ∈ hullInds || tessAreas[label2] > areaThresh
            return true # both peripheral 
        else 
            return false # only 1 peripheral 
        end
    else
        return false # 1 or neither peripheral 
    end
end


for (i, r) in enumerate(runs)

    subsetLengths = subset(lengthMeasurements, :File => m -> occursin.(r[1:end-6], m))
    lengthPerPixelEM = (subsetLengths[!, :length]./subsetLengths[!, :Pixels])[1]

    centroidData = load(datadir("exp_pro", "emCentroidsInteractive", "$(r[1:end-4]).jld2"))
    @unpack centroidLocations = centroidData

    imageIn = load(datadir("exp_pro", "cropped", r))
    imSize = size(imageIn)

    # Put centroid locations into a format for tessellation and triangulation 
    xs = [x[1] for x in centroidLocations]
    ys = [x[2] for x in centroidLocations]
    scalingFactor = maximum(imSize) / (1 - 3eps(Float64))
    shiftedCentroidLocations = centroidLocations ./ scalingFactor

    # Delaunay triangulation of centroid locations using function from GR
    n, tri = delaunay(xs, ys)
    # Count neighbours of each centroid in the triangulation 
    nNeighbours = [length(findall(x -> x == i, tri)) for i = 1:length(centroidLocations)]
    # Concave hull to identify boundary fibrils 
    hull = concave_hull(centroidLocations, 1)
    # Indices of fibrils within the hull 
    hullInds = sort([findall(x -> Point2(x...) == v, centroidLocations)[1] for v in hull.vertices])

    # Voronoi tessellation of centroid positions within (0,0) (1,1) box
    tess = voronoicells(shiftedCentroidLocations, Rectangle(Point2(0, 0), Point2(1, 1)))

    tessAreas = voronoiarea(tess)
    tessAreasFiltered = [tessAreas[a] for a in 1:length(tessAreas) if a ∉ hullInds]
    # display(tessAreasFiltered)        
    meanArea = mean(tessAreasFiltered)

    # Store number of neighbours without boundary elements
    nNeighboursDict[r] = nNeighbours[filter(x -> (x∉hullInds && tessAreas[x]<voronoiSizeThresh*meanArea), eachindex(nNeighbours))]

    pairs = Tuple[]
    for t in eachrow(tri)
        for i = 1:3
            labelSizeOrder = sort(t)
            if (labelSizeOrder[1], labelSizeOrder[2]) ∉ pairs
                if bothPeripheralTest(labelSizeOrder[1], labelSizeOrder[2], hullInds, tessAreas, voronoiSizeThresh * meanArea)
                    nothing
                else
                    push!(pairs, (labelSizeOrder[1], labelSizeOrder[2]))
                end
            end
            if (labelSizeOrder[1], labelSizeOrder[3]) ∉ pairs
                if bothPeripheralTest(labelSizeOrder[1], labelSizeOrder[3], hullInds, tessAreas, voronoiSizeThresh * meanArea)
                    nothing
                else
                    push!(pairs, (labelSizeOrder[1], labelSizeOrder[3]))
                end
            end
            if (labelSizeOrder[2], labelSizeOrder[3]) ∉ pairs
                if bothPeripheralTest(labelSizeOrder[2], labelSizeOrder[3], hullInds, tessAreas, voronoiSizeThresh * meanArea)
                    nothing
                else
                    push!(pairs, (labelSizeOrder[2], labelSizeOrder[3]))
                end
            end
        end
    end

    lengths = norm.(centroidLocations[first.(pairs)] - centroidLocations[last.(pairs)])
    lengths .*= lengthPerPixelEM * 1000.0
    pairLengthsDict[r] = lengths
end

fig1 = CairoMakie.Figure(resolution=(1000, 1000), fontsize=32)
ax1 = CairoMakie.Axis(fig1[1, 1])
for r in runs
    points = Point2[]
    hNorm = normalize(fit(Histogram, pairLengthsDict[r], 0:10:160), mode=:pdf)
    edgeVec = collect(hNorm.edges[1])
    for i = 1:length(hNorm.weights)
        push!(points, Point2(mean(edgeVec[i:i+1]), hNorm.weights[i]))
    end
    lines!(ax1, points)
end
ax1.xlabel = "Edge length/nm"
ax1.ylabel = "Density"
save(datadir("exp_pro", "emCentroidMeasurements", "emLengthHistogram$(Dates.format(Dates.now(), "yyyy-mm-dd-HH-MM-SS")).png"), fig1)


fig2 = CairoMakie.Figure(resolution=(1000, 1000), fontsize=32)
ax2 = CairoMakie.Axis(fig2[1, 1])
for r in runs
    # lines!(ax2,nNeighboursDict[r])
    points = Point2[]
    hNorm = normalize(fit(Histogram, nNeighboursDict[r], 0.5:11.5), mode=:pdf)
    edgeVec = collect(hNorm.edges[1])
    for i = 1:length(hNorm.weights)
        push!(points, Point2(mean(edgeVec[i:i+1]), hNorm.weights[i]))
    end
    lines!(ax2, points)
end
ax2.xlabel = "Neighbour count"
ax2.ylabel = "Density"
save(datadir("exp_pro", "emCentroidMeasurements", "emNeighbourHistogram$(Dates.format(Dates.now(), "yyyy-mm-dd-HH-MM-SS")).png"), fig2)