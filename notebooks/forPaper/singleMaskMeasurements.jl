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
using Statistics
using StatsBase
using Dates

@from "$(srcdir("ColourFunctions.jl"))" using ColourFunctions

function binariseSimulation!(uij)
    if uij < -0.5
        return 1.0
    else
        return 0.0
    end
end 

voronoiSizeThresh = 1.3

function bothPeripheralTest(label1, label2, hullInds, tessAreas, areaThresh)
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

# mask = "17tailT_4800X_HUI_0001_0"
runs = [r for r in Vector(readdlm(datadir("exp_pro", "filesToUse.txt"))[:, 1]) if !(occursin("mp13ko", r) || occursin("18tailT_4800X_HUI_0007_0", r) || occursin("18tailT_4800X_HUI_0008_0", r))]


lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","lengthMeasurements.csv")))
croppedLX = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","croppedLX.csv")))

for r in runs

    mask = r[1:end-4]

    # Collate simulation results as a dataframe 
    results = collect_results(datadir("fromCSF", "allMasksPhasespaceSeparateLengths", mask); subfolders=false)

    subsetCroppedLX = subset(croppedLX, :file => m -> m .== r)
    lengthPerPixelSim = subsetCroppedLX[1, :lX] ./ results[1, :nX]
    display(lengthPerPixelSim)

    subsetLengths = subset(lengthMeasurements, :File => m -> occursin.(r[1:end-6], m))
    lengthPerPixelEM = (subsetLengths[!, :length]./subsetLengths[!, :Pixels])[1]
    display(lengthPerPixelEM)

    centroidData = load(datadir("exp_pro", "emCentroidsInteractive", "$mask.jld2"))
    @unpack centroidLocations = centroidData
    centroidLocationsEM = centroidLocations
    imageIn = load(datadir("exp_pro", "cropped", r))
    imSize = size(imageIn)
    # Put centroid locations into a format for tessellation and triangulation 
    xsEM = [x[1] for x in centroidLocationsEM]
    ysEM = [x[2] for x in centroidLocationsEM]
    scalingFactorEM = maximum(imSize) / (1 - 3eps(Float64))
    shiftedCentroidLocationsEM = centroidLocationsEM ./ scalingFactorEM

    # Delaunay triangulation of centroid locations using function from GR
    nEM, triEM = delaunay(xsEM, ysEM)
    # Count neighbours of each centroid in the triangulation 
    nNeighboursEM = [length(findall(x -> x == i, triEM)) for i = 1:length(centroidLocationsEM)]
    # Concave hull to identify boundary fibrils 
    hullEM = concave_hull(centroidLocationsEM, 1)
    # Indices of fibrils within the hull 
    hullIndsEM = sort([findall(x -> Point2(x...) == v, centroidLocationsEM)[1] for v in hullEM.vertices])

    # Voronoi tessellation of centroid positions within (0,0) (1,1) box
    tessEM = voronoicells(shiftedCentroidLocationsEM, Rectangle(Point2(0, 0), Point2(1, 1)))
    tessAreasEM = voronoiarea(tessEM)
    tessAreasFilteredEM = [tessAreasEM[a] for a in 1:length(tessAreasEM) if a ∉ hullIndsEM]
    # display(tessAreasFiltered)        
    meanAreaEM = mean(tessAreasFilteredEM)

    # Store number of neighbours without boundary elements
    nNeighboursFilteredEM = nNeighboursEM[filter(x -> (x ∉ hullIndsEM && tessAreasEM[x] < voronoiSizeThresh * meanAreaEM), eachindex(nNeighboursEM))]

    pairsEM = Tuple[]
    for t in eachrow(triEM)
        for i = 1:3
            labelSizeOrder = sort(t)
            if (labelSizeOrder[1], labelSizeOrder[2]) ∉ pairsEM
                if bothPeripheralTest(labelSizeOrder[1], labelSizeOrder[2], hullIndsEM, tessAreasEM, voronoiSizeThresh * meanAreaEM)
                    nothing
                else
                    push!(pairsEM, (labelSizeOrder[1], labelSizeOrder[2]))
                end
            end
            if (labelSizeOrder[1], labelSizeOrder[3]) ∉ pairsEM
                if bothPeripheralTest(labelSizeOrder[1], labelSizeOrder[3], hullIndsEM, tessAreasEM, voronoiSizeThresh * meanAreaEM)
                    nothing
                else
                    push!(pairsEM, (labelSizeOrder[1], labelSizeOrder[3]))
                end
            end
            if (labelSizeOrder[2], labelSizeOrder[3]) ∉ pairsEM
                if bothPeripheralTest(labelSizeOrder[2], labelSizeOrder[3], hullIndsEM, tessAreasEM, voronoiSizeThresh * meanAreaEM)
                    nothing
                else
                    push!(pairsEM, (labelSizeOrder[2], labelSizeOrder[3]))
                end
            end
        end
    end

    lengthsEM = norm.(centroidLocationsEM[first.(pairsEM)] .- centroidLocationsEM[last.(pairsEM)])
    lengthsEM .*= lengthPerPixelEM*1000.0

    pairLengthsDict = Dict()
    nNeighboursDict = Dict()

    # Loop to process data from each run 
    for i = 1:nrow(results)

        # Convert simulation result to a 2D matrix
        uMat = reshape(results[i, :u][end], (results[i, :nY], results[i, :nX]))
        # Binarise grayscale image 
        uImg = binariseSimulation!.(uMat)
        # Convert matrix to a grayscale image
        uGray = Gray.(uImg)
        # Segment binarised image
        seg = fast_scanning(uGray, 0.01)
        # display(map(i->get_random_color(i), seg.image_indexmap))
        # Find centre of mass positions of all fibril segments. 
        centroidLocationsSim = Point2{Float64}[]
        maxSize = 500
        # seg = prune_segments(seg, i->(segment_pixel_count(seg,i)<20), (i,j)->(-segment_pixel_count(seg,j)))
        seg = prune_segments(seg, i -> (segment_pixel_count(seg, i) > maxSize), (i, j) -> (-segment_pixel_count(seg, j)))
        for k in seg.segment_labels
            pixels = findall(x -> x == k, seg.image_indexmap)
            com = Tuple(sum(pixels)) ./ length(pixels)
            # Exclude the one remaining segment above size of maxSize, representing the system background
            if length(pixels) < maxSize
                push!(centroidLocationsSim, Point2{Float64}(com[2], -com[1]))
            end
        end


        # Put centroid locations into a format for tessellation and triangulation 
        xsSim = [x[1] for x in centroidLocationsSim]
        ysSim = [x[2] for x in centroidLocationsSim]
        scalingFactorSim = maximum(size(uMat))/(1 - 3eps(Float64))
        shiftedCentroidLocationsSim = centroidLocationsSim./scalingFactorSim
        shiftedCentroidLocationsSim .+= Point2(0,1)

        # Delaunay triangulation of centroid locations using function from GR
        nSim, triSim = delaunay(xsSim, ysSim)

        # Concave hull to identify boundary fibrils 
        hullSim = concave_hull(centroidLocationsSim, 1)
        # Indices of fibrils within the hull 
        hullIndsSim = sort([findall(x -> Point2(x...) == v, centroidLocationsSim)[1] for v in hullSim.vertices])
        # Voronoi tessellation of centroid positions within (0,0) (1,1) box
        tessSim = voronoicells(shiftedCentroidLocationsSim, Rectangle(Point2(0, 0), Point2(1, 1)))
        tessAreasSim = voronoiarea(tessSim)
        tessAreasFilteredSim = [tessAreasSim[a] for a in 1:length(tessAreasSim) if a ∉ hullIndsSim]
        # display(tessAreasFiltered)        
        meanAreaSim = mean(tessAreasFilteredSim)


        # Count neighbours of each centroid in the triangulation 
        nNeighboursSim = [length(findall(x -> x == i, triSim)) for i = 1:length(centroidLocationsSim)]
        nNeighboursFilteredSim = nNeighboursSim[filter(x -> (x ∉ hullIndsSim && tessAreasSim[x] < voronoiSizeThresh*meanAreaSim), eachindex(nNeighboursSim))]
        
        nNeighboursDict[(results[i, :r], results[i, :ϕ0])] = nNeighboursFilteredSim

        pairsSim = Tuple[]
        for t in eachrow(triSim)
            for i = 1:3
                labelSizeOrder = sort(t)
                if (labelSizeOrder[1], labelSizeOrder[2]) ∉ pairsSim
                    if bothPeripheralTest(labelSizeOrder[1], labelSizeOrder[2], hullIndsSim, tessAreasSim, voronoiSizeThresh*meanAreaSim)
                        nothing
                    else
                        push!(pairsSim, (labelSizeOrder[1], labelSizeOrder[2]))
                    end
                end
                if (labelSizeOrder[1], labelSizeOrder[3]) ∉ pairsSim
                    if bothPeripheralTest(labelSizeOrder[1], labelSizeOrder[3], hullIndsSim, tessAreasSim, voronoiSizeThresh*meanAreaSim)
                        nothing
                    else
                        push!(pairsSim, (labelSizeOrder[1], labelSizeOrder[3]))
                    end
                end
                if (labelSizeOrder[2], labelSizeOrder[3]) ∉ pairsSim
                    if bothPeripheralTest(labelSizeOrder[2], labelSizeOrder[3], hullIndsSim, tessAreasSim, voronoiSizeThresh*meanAreaSim)
                        nothing
                    else
                        push!(pairsSim, (labelSizeOrder[2], labelSizeOrder[3]))
                    end
                end
            end
        end

        lengthsSim = norm.(centroidLocationsSim[first.(pairsSim)] .- centroidLocationsSim[last.(pairsSim)])
        lengthsSim .*= lengthPerPixelSim*1000.0
        pairLengthsDict[(results[i, :r], results[i, :ϕ0])] = lengthsSim

    end


    fig1 = CairoMakie.Figure(resolution=(500, 500), fontsize=32)
    ax1 = CairoMakie.Axis(fig1[1, 1])
    for r in keys(pairLengthsDict)
        points = Point2[]
        hNorm = normalize(fit(Histogram, pairLengthsDict[r], 0:10:160), mode=:pdf)
        edgeVec = collect(hNorm.edges[1])
        for i = 1:length(hNorm.weights)
            push!(points, Point2(mean(edgeVec[i:i+1]), hNorm.weights[i]))
        end
        lines!(ax1, points)
        # display(median(pairLengthsDict[r]))
    end
    points = Point2[]
    hNorm = normalize(fit(Histogram, lengthsEM, 0:10:160), mode=:pdf)
    edgeVec = collect(hNorm.edges[1])
    for i = 1:length(hNorm.weights)
        push!(points, Point2(mean(edgeVec[i:i+1]), hNorm.weights[i]))
    end
    lines!(ax1, points, linestyle="--", color=:black, width=4)
    # display(median(lengthsEM))
    # vlines!(ax1,mean(lengthsEM))
    ax1.xlabel = "Edge length/nm"
    ax1.ylabel = "Density"
    save(datadir("fromCSF", "allMasksPhasespaceSeparateLengths", mask, "$(mask)Length.png"), fig1)


    fig2 = CairoMakie.Figure(resolution=(500, 500), fontsize=32)
    ax2 = CairoMakie.Axis(fig2[1, 1])
    for r in keys(nNeighboursDict)
        # lines!(ax2,nNeighboursDict[r])
        points = Point2[]
        hNorm = normalize(fit(Histogram, nNeighboursDict[r], 0.5:11.5), mode=:pdf)
        edgeVec = collect(hNorm.edges[1])
        for i = 1:length(hNorm.weights)
            push!(points, Point2(mean(edgeVec[i:i+1]), hNorm.weights[i]))
        end
        lines!(ax2, points)
    end
    points = Point2[]
    hNorm = normalize(fit(Histogram, nNeighboursFilteredEM, 0.5:11.5), mode=:pdf)
    edgeVec = collect(hNorm.edges[1])
    for i = 1:length(hNorm.weights)
        push!(points, Point2(mean(edgeVec[i:i+1]), hNorm.weights[i]))
    end
    lines!(ax2, points, linestyle="--", color=:black, width=4)
    ax2.xlabel = "Neighbour count"
    ax2.ylabel = "Density"
    save(datadir("fromCSF", "allMasksPhasespaceSeparateLengths", mask, "$(mask)Neighbour.png"), fig2)

end


