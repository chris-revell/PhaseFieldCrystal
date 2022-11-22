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
using DelimitedFiles

@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions

function neighbourColours(x)
    if x==6
        return (:white,0.0)
    elseif x==5 
        return (:red,0.75)
    elseif x==7
        return (:blue,0.75)
    else 
        return (:grey,0.75)
    end
end

mkpath(datadir("exp_pro","emCentroidNeighbours"))

runs = [r for r in Vector(readdlm(datadir("exp_pro","filesToUse.txt"))[:,1]) if !occursin("mmpko",r)]

croppedLX = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","croppedLX.csv")))

fig = Figure(resolution=(6000,6000),backgroundcolor=:white,fontsize=64)

axes = Dict()
sizes = Dict()

# Container to store defect counts 
defectCountsDataFrame = DataFrame()

for (k,r) in enumerate(runs)

    maskData = load(datadir("exp_pro","masks",r[1:end-4],"$(r[1:end-4]).jld2"))
    @unpack newIndexMap, lX, h = maskData

    centroidData = load(datadir("exp_pro","emCentroidsInteractive","$(r[1:end-4]).jld2"))
    @unpack centroidLocations = centroidData

    imageIn = load(datadir("exp_pro","cropped",r))
    imSize = size(imageIn)
    grayImage = Gray.(imageIn)

    croppedLXRow  = filter(:file => f->f==r, croppedLX)
    lengthPerPixel = croppedLXRow[1,:lX]/size(imageIn)[2]

    # Put centroid locations into a format for tessellation and triangulation 
    xs = [x[1] for x in centroidLocations]
    ys = [x[2] for x in centroidLocations]
    scalingFactor = maximum(size(newIndexMap))/(1-3eps(Float64))
    shiftedCentroidLocations = centroidLocations./scalingFactor

    # Delaunay triangulation of centroid locations using function from GR
    n, tri = delaunay(xs,ys)
    # Count neighbours of each centroid in the triangulation 
    nNeighbours = [length(findall(x->x==i,tri)) for i=1:length(shiftedCentroidLocations)]

    # Concave hull to identify boundary fibrils 
    hull = concave_hull(shiftedCentroidLocations)
    # Indices of fibrils within the hull 
    hullInds = sort([findall(x->Point2(x...)==v,shiftedCentroidLocations)[1] for v in hull.vertices])

    defectCountsDict = Dict()
    for i in eachindex(nNeighbours)
        if i ∉ hullInds
            if "$(nNeighbours[i])" ∈ keys(defectCountsDict)
                defectCountsDict["$(nNeighbours[i])"] += 1
            else
                defectCountsDict["$(nNeighbours[i])"] = 1
            end
        end
    end
    dfLocal = DataFrame(defectCountsDict)
    dfLocal[!,:file] = [r]
    defectCountsDataFrame = vcat(defectCountsDataFrame,dfLocal,cols = :union)

    runDefectProportion = 1-defectCountsDict["6"]/(length(nNeighbours)-length(hullInds))
    
    # Voronoi tessellation of centroid positions within (0,0) (1,1) box
    tess = voronoicells(shiftedCentroidLocations, Rectangle(Point2(0, 0), Point2(1, 1)))

    ax2 = CairoMakie.Axis(fig[(k-1)%6+1,(k-1)÷6+1],aspect=DataAspect(), backgroundcolor=:white)
    image!(ax2,rotr90(grayImage))    
    for (i,c) in enumerate(tess.Cells)
        if i ∉ hullInds
            vertices = [v.*scalingFactor for v in c]
            poly!(ax2, vertices, color=neighbourColours(nNeighbours[i]),strokecolor=(:black,1.0),strokewidth=1.0)
        # else
        #     vertices = [v.*scalingFactor for v in c]
        #     poly!(ax2, vertices, color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1.0)
        end
    end
    hideydecorations!(ax2)
    ax2.xticks = [0,size(grayImage)[2]]
    ax2.xtickformat = x -> string.(round.((Float64.(x).*lengthPerPixel),digits=2))
    # ax2.yticks = []
    # ax2.ytickformat = y -> string.(round.((Float64.(y).*lengthPerPixel),digits=2))
    hidespines!(ax2)
    
    Label(fig[(k-1)%6+1,(k-1)÷6+1, Bottom()], "$k, d=$(round(runDefectProportion,digits=2))", padding = (0, 10, 10, 0), color=:black)

    axes[r] = ax2
    sizes[r] = size(grayImage)
end 

# lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","lengthMeasurements.csv")))
# lengthPerPixel = lengthMeasurements[!,:length]./lengthMeasurements[!,:Pixels]
# lengthPerPixelDict = Dict()
# lengthDict = Dict()
# for r in runs 
#     subsetLengths = subset(lengthMeasurements, :File => m -> occursin.(r[1:end-6],m))
#     lengthPerPixelDict[r] = (subsetLengths[!,:length]./subsetLengths[!,:Pixels])[1]
#     lengthDict[r] = lengthPerPixelDict[r].*sizes[r]
# end

# xMax = maximum(first.(values(lengthDict)))
# yMax = maximum(last.(values(lengthDict)))

# for r in runs 
#     xlims!(axes[r],(0,xMax)./lengthPerPixelDict[r])
#     ylims!(axes[r],(0,yMax)./lengthPerPixelDict[r])
# end

# colgap!(fig.layout, 1, -700)
# colgap!(fig.layout, 2, -600)
# colgap!(fig.layout, 3, -600)
# colgap!(fig.layout, 4, -500)
# colgap!(fig.layout, 5, -100)

# rowgap!(fig.layout, 1, -100)
# rowgap!(fig.layout, 2, -100)
# rowgap!(fig.layout, 3, -200)

resize_to_layout!(fig)
display(fig)

# CSV.write(datadir("exp_pro", "emCentroidMeasurements", "defectCounts.csv"), defectCountsDataFrame)

save(datadir("exp_pro","emCentroidNeighbours","emNeighboursGridWithDefectProportion.png"),fig)



