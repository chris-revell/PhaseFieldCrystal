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

@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions

function neighbourColours(x)
    if x==6
        return (:white,0.0)
    elseif x==5 
        return (:red,1.0)
    elseif x==7
        return (:blue,1.0)
    else 
        return (:black,1.0)
    end
end

function binariseSimulation!(uij)
    if uij > 0.5
        return 1.0
    else
        return 0.0
    end
end 

# mask = "17tailT_4800X_HUI_0001_0"
runs = Vector(readdlm(datadir("exp_pro","filesToUse.txt"))[:,1])

defectCountsDataFrame = DataFrame()

for r in runs
    mask = r[1:end-4]

    # Collate results as a dataframe 
    results = collect_results(datadir("fromCSF","allMasksPhasespaces",mask); subfolders = false)

    # Loop to process data from each run 
    for i=1:nrow(results)
        
        # Convert simulation result to a 2D matrix
        uMat = reshape(results[i,:u][end],(results[i,:nY],results[i,:nX]))
        # Binarise grayscale image
        uImg = binariseSimulation!.(uMat)
        # Convert matrix to a grayscale image
        uGray = Gray.(uImg)
        # Segment binarised image
        seg = fast_scanning(uGray, 0.01)    
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
        # display(map(i->get_random_color(i), seg.image_indexmap))

        # Put centroid locations into a format for tessellation and triangulation 
        xs = [x[1] for x in centroidLocations]
        ys = [x[2] for x in centroidLocations]
        scalingFactor = maximum(abs.([xs ys]))/(1-3eps(Float64))
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
                if "$(nNeighbours[j])" ∈ keys(defectCountsDict)
                    defectCountsDict["$(nNeighbours[j])"] += 1
                else
                    defectCountsDict["$(nNeighbours[j])"] = 1
                end
            end
        end
        dfLocal = DataFrame(defectCountsDict)
        dfLocal[!,:file] = [r]
        dfLocal[!,:ϕ0]   = [results[i,:ϕ0]]
        dfLocal[!,:r]   = [results[i,:r]]
        defectCountsDataFrame = vcat(defectCountsDataFrame,dfLocal,cols = :union)

        runDefectProportion = 1-defectCountsDict["6"]/(length(nNeighbours)-length(hullInds))



    end 

  
end


