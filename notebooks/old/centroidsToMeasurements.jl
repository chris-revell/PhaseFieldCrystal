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
using Statistics
using StatsBase

@from "$(projectdir("src","ColourFunctions.jl"))" using ColourFunctions

mkpath(datadir("exp_pro","emCentroidMeasurements"))

runs = Vector(readdlm(datadir("exp_pro","filesToUse.txt"))[:,1])

fig1 = Figure(resolution=(6000,6000),backgroundcolor=:white,fontsize=64)

lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","lengthMeasurements.csv")))
lengthPerPixel = lengthMeasurements[!,:length]./lengthMeasurements[!,:Pixels]
lengthPerPixelDict = Dict()
for r in runs 
    subsetLengths = subset(lengthMeasurements, :File => m -> occursin.(r[1:end-6],m))
    lengthPerPixelDict[r] = (subsetLengths[!,:length]./subsetLengths[!,:Pixels])[1]
end

axes = Dict()
sizes = Dict()
imLengths = Dict()

for (i,r) in enumerate(runs)

    imageIn = load(datadir("exp_pro","cropped",r))
    
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

    ax2 = CairoMakie.Axis(fig1[mod(i-1,6)+1,(i-1)÷6+1],aspect=DataAspect(),backgroundcolor=:white)
    image!(ax2,rotr90(imageIn))
    #scatter!(ax2,centroidLocations,color=(:orange,1.0),markersize=ceil(Int64,10000/imSize[1]))
    scatter!(ax2,centroidLocations,color=(:orange,1.0),markersize=10)
    for p in pairs
        lines!(ax2,Point.(centroidLocations[[p[1],p[2]]]),width=2,color=:orange)
    end
    hidedecorations!(ax2)
    hidespines!(ax2)
    Label(fig1[mod(i-1,6)+1,(i-1)÷6+1, Bottom()], L"\mu=%$(round(mean(lengths),digits=1))nm,~\sigma=%$(round(std(lengths),digits=1))", valign = :bottom, font = "TeX Gyre Heros Bold", padding = (0, 10, 10, 0))
    
    axes[r] = ax2
    sizes[r] = size(imageIn)
    imLengths[r] = lengths
end 

lengthDict = Dict()
for r in runs 
    lengthDict[r] = lengthPerPixelDict[r].*sizes[r]
end

xMax = maximum(first.(values(lengthDict)))
yMax = maximum(last.(values(lengthDict)))

for r in runs 
    xlims!(axes[r],(0,xMax)./lengthPerPixelDict[r])
    ylims!(axes[r],(0,yMax)./lengthPerPixelDict[r])
end

resize_to_layout!(fig1)
display(fig1)
save(datadir("exp_pro","emCentroidMeasurements","grid.png"),fig1)


means = Float64[]
stds = Float64[]
for r in runs
    push!(means,mean(imLengths[r]))
    push!(stds,std(imLengths[r]))
end 
spacingData = DataFrame(file=runs, mean=means, std=stds)
CSV.write(datadir("exp_pro", "emCentroidMeasurements", "spacingData.csv"), spacingData)

fig3 = Figure(resolution=(1000,1000),backgroundcolor=:white,fontsize=32)
axErrorBars = CairoMakie.Axis(fig3[1,1])
errorbars!(axErrorBars,Point2.(collect(1:length(runs)), means),stds./2, color=:black)
scatter!(axErrorBars,Point2.(collect(1:length(runs)), means), color=:black)
axErrorBars.xlabel="File"
axErrorBars.ylabel="Mean spacing/nm"
axErrorBars.xticks = 1:length(runs)
axErrorBars.xtickformat = x -> runs[Int.(x)]
axErrorBars.xticklabelrotation = -π/2
resize_to_layout!(fig3)
display(fig3)
save(datadir("exp_pro","emCentroidMeasurements","means.png"),fig3)

fig2 = Figure(resolution=(1000,1000),backgroundcolor=:white,fontsize=32)
axHist = CairoMakie.Axis(fig2[1,1])
for r in runs
    points = Point2[]
    # h = fit(Histogram, imLengths[r])
    hNorm = normalize(fit(Histogram, imLengths[r]), mode=:pdf)
    display(hNorm.weights)
    edgeVec = collect(hNorm.edges[1])    
    for i=1:length(hNorm.weights)
        push!(points,Point2(mean(edgeVec[i:i+1]),hNorm.weights[i]))
    end
    lines!(axHist,points)
end
axHist.xlabel = L"Edge~length/nm"
axHist.ylabel = L"Density"
display(fig2)
save(datadir("exp_pro","emCentroidMeasurements","lengthsHistogram.png"),fig2)


