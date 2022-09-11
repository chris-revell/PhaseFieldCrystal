using DrWatson
@quickactivate
using UnPack
using CairoMakie; set_theme!(figure_padding=0,fontsize=32)
using ColorSchemes
using DifferentialEquations
using JLD2
using LaTeXStrings
using DataFrames
using PerceptualColourMaps
using Images
using ImageBinarization
using ImageSegmentation

folderPath = "PhaseSpace3"

results = collect_results!(datadir("fromCSF",folderPath); subfolders = true)

fig = Figure(resolution=(6000,6000))
axes = Dict()

minSize = 2000
maxSize = 10000

for (i,result) in enumerate(eachrow(subset(results, :m => m -> m.== 0.1)))
    ax = CairoMakie.Axis(fig[i%6+1,i÷6+1],aspect=DataAspect())
    ax.yreversed = true
    uMat = transpose(reshape(result[:u],(result[:nY],result[:nX])))
    heatmap!(ax,uMat,colorrange=(-1.0, 1.0),colormap=cmap("D1"))
    hidedecorations!(ax)
    hidespines!(ax)    
    axes[(result[:r],result[:ϕ0])] = ax

    # uImg = Gray.(uMat)

    # uBinarized = binarize(uImg,Otsu())

    # seg = felzenszwalb(uBinarized, 10)    
    # seg2 = prune_segments(seg, i->(segment_pixel_count(seg,i)<minSize), (i,j)->(segment_pixel_count(seg,j)))
    # seg3 = prune_segments(seg2, i->(segment_pixel_count(seg2,i)>maxSize), (i,j)->(segment_pixel_count(seg2,j)))
    # # segmentedImage3 = map(i->maskColour(i,seg3), labels_map(seg3))

    # centroidLocations = Point2{Float64}[]
    # for k in seg3.segment_labels
    #     pixels = zeros(2)
    #     count = 0
    #     for i=1:size(image)[1]
    #         for j=1:size(image)[2]
    #             if seg3.image_indexmap[i,j] == k
    #                 pixels .+= [j,-i]
    #                 count += 1
    #             end
    #         end
    #     end
    #     if count<1000
    #         push!(centroidLocations,Point2{Float64}(pixels./count...))
    #     end
    # end
    # centroidLocations .+= Point2(0.0,size(iImg)[1]*1.0)
    # xs = [x[1] for x in centroidLocations]
    # ys = [x[2] for x in centroidLocations]
    # scalingFactor = maximum([xs ys])/(1-3eps(Float64))
    # shiftedCentroidLocations = centroidLocations./scalingFactor

    # xlims!(ax,0,size(image)[2]/scalingFactor)
    # ylims!(ax,0,size(image)[1]/scalingFactor)
    # scatter!(ax,shiftedCentroidLocations,color=:green,markersize=10)
    # # lims=(minimum(areas),maximum(areas))
    # for i=1:n
    #     poly!(ax,shiftedCentroidLocations[tri[i,:]],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1.0)
    # end


end 

sortedrs  = unique!(sort(first.(keys(axes))))
sortedϕ0s = unique!(sort(last.(keys(axes))))

for (i,r) in enumerate(sortedrs)
    for (j,ϕ0) in enumerate(sortedϕ0s)
        fig[length(sortedrs)+1-i,length(sortedϕ0s)+1-j][1,1] = axes[(r,ϕ0)]
        # Colorbar(fig[length(sortedrs)-i,length(sortedϕ0s)-j][1,2], limits=(-1,1),colormap=:bwr)#, size = 25)
        Label(fig[length(sortedrs)+1-i,length(sortedϕ0s)+1-j,Bottom()],L"r=%$r, \phi_0=%$ϕ0",textsize=64)
    end
end

for i=1:length(sortedϕ0s)
    colsize!(fig.layout,i,Aspect(1.0,1.0))
end
resize_to_layout!(fig.layout)

display(fig)





save(datadir("fromCSF",folderPath,"phasespace_m=0.1.png"),fig)


