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
using LaTeXStrings
using StatsBase
using Statistics

@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions

function binariseSimulation!(uij)
    if uij < -0.5
        return Gray(1.0)
    else
        return Gray(0.0)
    end
end 

# mask = "17tailT_4800X_HUI_0001_0"
runs = Vector(readdlm(datadir("exp_pro","filesToUse.txt"))[:,1])

lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","lengthMeasurements.csv")))
lengthPerPixel = lengthMeasurements[!,:length]./lengthMeasurements[!,:Pixels]
lengthPerPixelDict = Dict()
for r in runs 
    subsetLengths = subset(lengthMeasurements, :File => m -> occursin.(r[1:end-6],m))
    lengthPerPixelDict[r] = (subsetLengths[!,:length]./subsetLengths[!,:Pixels])[1]
end

croppedLX = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","croppedLX.csv")))

for r in runs
    mask = r[1:end-4]

    maskIn = load(datadir("exp_pro","masksCompressed",mask,"$mask.png"))
    nX = size(maskIn)[2]
    imageIn = load(datadir("exp_pro","cropped","$mask.png"))

    subsetCroppedLX = filter(:file => m -> occursin.(mask,m), croppedLX)
    lengthPerPixelMask = subsetCroppedLX[1,:lX]/nX

    # Collate results as a dataframe 
    results = collect_results(datadir("fromCSF","allMasksPhasespaces2",mask); subfolders = false)

    fig = Figure(resolution=(2000,1000),fontsize=32)
    ax1 = CairoMakie.Axis(fig[1,2])
    ax2 = CairoMakie.Axis(fig[1,3])
    ax1.title = "Spacing distribution"
    ax1.xlabel = "Edge length/nm"
    ax1.ylabel = "Density"
    ax2.title = "Neighbour distribution"
    ax2.xlabel = "Neighbour count"
    ax2.ylabel = "Density"
    ax3 = CairoMakie.Axis(fig[1,1],aspect=DataAspect())
    hidedecorations!(ax3)
    hidespines!(ax3)
    image!(rotr90(imageIn))
    # ax3 = CairoMakie.Axis(fig[2,2])
    # ax4 = CairoMakie.Axis(fig[2,3])
    # ax3.title = "Em spacing"
    # ax3.xlabel = "Edge length/nm"
    # ax3.ylabel = "Density"
    # ax4.title = "EM neighbour distribution"
    # ax4.xlabel = "Neighbour count"
    # ax4.ylabel = "Density"

    # ylims!(ax1,(0,1))
    # ylims!(ax2,(0,1))
    # ylims!(ax3,(0,1))
    # ylims!(ax4,(0,1))

    # ax3 = CairoMakie.Axis(fig[2,2])
    # ax4 = CairoMakie.Axis(fig[2,3])

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
    lengths = norm.(centroidLocations[first.(pairs)]-centroidLocations[last.(pairs)])
    lengths .*= lengthPerPixelDict[r]*1000.0
    
    points = Point2[]
    hNorm = normalize(fit(Histogram, lengths, 0:10:160), mode=:pdf)
    edgeVec = collect(hNorm.edges[1])    
    for i=1:length(hNorm.weights)
        push!(points,Point2(mean(edgeVec[i:i+1]),hNorm.weights[i]))
    end
    lines!(ax1,points,linestyle="--",color=:black,label="EM data")
    points = Point2[]
    hNorm = normalize(fit(Histogram, nNeighbours[filter(x->x ∉ hullInds, eachindex(nNeighbours))], 0.5:11.5), mode=:pdf)
    edgeVec = collect(hNorm.edges[1])    
    for i=1:length(hNorm.weights)
        push!(points,Point2(mean(edgeVec[i:i+1]),hNorm.weights[i]))
    end
    lines!(ax2,points,linestyle="--",color=:black,label="EM data")
        
    # Loop to process data from each run 
    for i=1:nrow(results)
        
        # Convert simulation result to a 2D matrix
        uMat = reshape(results[i,:u][end],(results[i,:nY],results[i,:nX]))
        # Binarise grayscale image
        uImg = binariseSimulation!.(uMat)
        # Segment binarised image
        seg = fast_scanning(uImg, 0.01)    
        # Find centre of mass positions of all fibril segments. 
        centroidLocations = Point2{Float64}[]
        maxSize = 500
        for k in seg.segment_labels
            pixels = findall(x->x==k,seg.image_indexmap)
            com = Tuple(sum(pixels))./length(pixels)        
            # Exclude the one remaining segment above size of maxSize, representing the system background
            if length(pixels)<maxSize
                push!(centroidLocations,Point2{Float64}(com[2],-com[1]))
            end
        end

        # Put centroid locations into a format for tessellation and triangulation 
        xs = [x[1] for x in centroidLocations]
        ys = [x[2] for x in centroidLocations]
        scalingFactor = maximum(size(uMat))/(1-3eps(Float64))
        shiftedCentroidLocations = centroidLocations./scalingFactor
        shiftedCentroidLocations .+= Point2(0,1)

        # Delaunay triangulation of centroid locations using function from GR
        n, tri = delaunay(xs,ys)
        # Count neighbours of each centroid in the triangulation 
        nNeighbours = [length(findall(x->x==i,tri)) for i=1:length(shiftedCentroidLocations)]

        # Concave hull to identify boundary fibrils 
        hull = concave_hull(shiftedCentroidLocations,1)
        # Indices of fibrils within the hull 
        hullInds = sort([findall(x->Point2(x...)==v,shiftedCentroidLocations)[1] for v in hull.vertices])

        # Voronoi tessellation of centroid positions within (0,0) (1,1) box
        tess = voronoicells(shiftedCentroidLocations, Rectangle(Point2(0, 0), Point2(1, 1)))

        defectCountsDict = Dict()
        for j in eachindex(nNeighbours)
            if j ∉ hullInds
                if nNeighbours[j] ∈ keys(defectCountsDict)
                    defectCountsDict[nNeighbours[j]] += 1
                else
                    defectCountsDict[nNeighbours[j]] = 1
                end
            end
        end
        
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
        
        lengths = norm.(centroidLocations[first.(pairs)]-centroidLocations[last.(pairs)])
        lengths .*= lengthPerPixelMask*1000.0
        
        points = Point2[]
        hNorm = normalize(fit(Histogram, lengths, 0:10:160), mode=:pdf)
        edgeVec = collect(hNorm.edges[1])    
        for i=1:length(hNorm.weights)
            push!(points,Point2(mean(edgeVec[i:i+1]),hNorm.weights[i]))
        end
        lines!(ax1,points)

        points = Point2[]
        hNorm = normalize(fit(Histogram, nNeighbours[filter(x->x ∉ hullInds, eachindex(nNeighbours))], 0.5:11.5), mode=:pdf)
        edgeVec = collect(hNorm.edges[1])    
        for i=1:length(hNorm.weights)
            push!(points,Point2(mean(edgeVec[i:i+1]),hNorm.weights[i]))
        end
        lines!(ax2,points)

        display(mean(lengths))
    end 

    axislegend(ax1, position = :rt)
    axislegend(ax2, position = :rt)

    resize_to_layout!(fig)

    save(datadir("fromCSF","allMasksPhasespaces2",mask,"$(mask)Histograms.png"),fig)
end


