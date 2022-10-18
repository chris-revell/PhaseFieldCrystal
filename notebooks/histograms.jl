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

runs = [f for f in readdir(datadir("exp_pro","masks","ok")) if f[end-3:end]==".png"]

# spacingData = DataFrame(CSV.File(datadir("exp_pro","emCentroidMeasurements","spacingData.csv")))
# defectCounts = DataFrame(CSV.File(datadir("exp_pro","emCentroidMeasurements","defectCounts.csv")))

lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","lengthMeasurements.csv")))
lengthPerPixel = lengthMeasurements[!,:length]./lengthMeasurements[!,:Pixels]
lengthPerPixelDict = Dict()
for r in runs 
    subsetLengths = subset(lengthMeasurements, :File => m -> occursin.(r[1:end-6],m))
    lengthPerPixelDict[r] = (subsetLengths[!,:length]./subsetLengths[!,:Pixels])[1]
end

imLengths = Dict()
defectCountsDict = Dict()

for (i,r) in enumerate(runs)
    centroidData = load(datadir("exp_pro","emCentroidsInteractive","$(r[1:end-4]).jld2"))
    @unpack centroidLocations = centroidData
    # Put centroid locations into a format for tessellation and triangulation 
    xs = [x[1] for x in centroidLocations]
    ys = [x[2] for x in centroidLocations]
    scalingFactor = maximum(abs.([xs ys]))/(1-3eps(Float64))
    shiftedCentroidLocations = centroidLocations./scalingFactor
    # Delaunay triangulation of centroid locations using function from GR
    n, tri = delaunay(xs,ys)
    # Count neighbours of each centroid in the triangulation 
    nNeighbours = [length(findall(x->x==i,tri)) for i=1:length(shiftedCentroidLocations)]
    # Concave hull to identify boundary fibrils 
    hull = concave_hull(shiftedCentroidLocations,1)
    # Indices of fibrils within the hull 
    hullInds = sort([findall(x->Point2(x...)==v,shiftedCentroidLocations)[1] for v in hull.vertices])
    pairs = Tuple[]
    for t in eachrow(tri)
        for i=1:3
            labelSizeOrder = sort(t)
            if (labelSizeOrder[1],labelSizeOrder[2]) ∉ pairs 
                if labelSizeOrder[1] ∈ hullInds && labelSizeOrder[2] ∈ hullInds
                    nothing 
                else 
                    push!(pairs,(labelSizeOrder[1],labelSizeOrder[2]))
                end
            end 
            if (labelSizeOrder[1],labelSizeOrder[3]) ∉ pairs 
                if labelSizeOrder[1] ∈ hullInds && labelSizeOrder[3] ∈ hullInds
                    nothing 
                else 
                    push!(pairs,(labelSizeOrder[1],labelSizeOrder[3]))
                end
            end 
            if (labelSizeOrder[2],labelSizeOrder[3]) ∉ pairs 
                if labelSizeOrder[2] ∈ hullInds && labelSizeOrder[3] ∈ hullInds
                    nothing 
                else 
                    push!(pairs,(labelSizeOrder[2],labelSizeOrder[3]))
                end
            end
        end
    end 

    localCountDict = Dict()
    for i in nNeighbours
        if i ∉ hullInds
            if i ∈ keys(localCountDict)
                localCountDict[i] += 1
            else
                localCountDict[i] = 1
            end
        end
    end
    pointsVec = Point2[]
    for k in keys(localCountDict)
        push!(pointsVec,Point2(k,localCountDict[k]))
    end
    defectCountsDict[r] = pointsVec

    # dfLocal = DataFrame(defectCountsDict)
    # dfLocal[!,:file] = [r]
    # defectCountsDataFrame = vcat(defectCountsDataFrame,dfLocal,cols = :union)

    lengths = norm.(centroidLocations[first.(pairs)]-centroidLocations[last.(pairs)])
    lengths .*= lengthPerPixelDict[r]*1000.0
    imLengths[r] = lengths
end 

fig = CairoMakie.Figure(resolution=(2000,1000))
ax1 = CairoMakie.Axis(fig[1,1])
for r in runs
    points = Point2[]
    hNorm = normalize(fit(Histogram, imLengths[r]), mode=:pdf)
    edgeVec = collect(hNorm.edges[1])    
    for i=1:length(hNorm.weights)
        push!(points,Point2(mean(edgeVec[i:i+1]),hNorm.weights[i]))
    end
    lines!(ax1,points)
end
ax1.xlabel = L"Edge~length/nm"
ax1.ylabel = L"Density"

ax2 = CairoMakie.Axis(fig[1,2])
for r in runs
    lines!(ax2,defectCountsDict[r])
end

display(fig)