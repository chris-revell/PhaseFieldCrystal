using DrWatson; @quickactivate
using DataFrames
using CairoMakie
using Colors

runs = [f for f in readdir(datadir("fromCSF","lambdaTest")) if f!=".DS_Store"]

for mask in runs

    maskIn = load(datadir("exp_pro","masksCompressed",mask,"$mask.png"))
    maskImage = fill(RGBA(1,1,1,1),size(maskIn))
    for i=1:size(maskIn)[1]
        for j=1:size(maskIn)[2]
            if maskIn[i,j] < 0.5
                maskImage[i,j] = RGBA(0.0,0.0,0.0,1.0)
            else
                maskImage[i,j] = RGBA(0.0,0.0,0.0,0.0)
            end
        end
    end    

    res = collect_results(datadir("fromCSF","lambdaTest",mask))

    fig = Figure(resolution=(3000,1000),fontsize=64)

    axesDict = Dict()
    for (i,r) in enumerate(eachrow(subset(res, :r => r -> r .== 0.75)))
        ax = Axis(fig[1,i],aspect=DataAspect())
        uMat = reshape(r[:u][end],(r[:nY],r[:nX]))
        heatmap!(ax,rotr90(uMat),colorrange=(-1.0, 1.0),colormap=:bwr)
        hidedecorations!(ax)
        hidespines!(ax)
        # Map parameters to axis in axes dictionary 
        axesDict[r[:λ]] = ax
        image!(ax,rotr90(maskImage))
    end
    sortedλs  = unique!(sort(first.(keys(axesDict))))    
    for (i,λ) in enumerate(sortedλs)
        fig[1,i] = axesDict[λ]        
    end
    Label(fig[1,1,Left()],"r=0.75\n ϕ₀=0.42")

    axesDict = Dict()
    for (i,r) in enumerate(eachrow(subset(res, :r => r -> r .== 0.55)))
        ax = Axis(fig[1,i],aspect=DataAspect())
        uMat = reshape(r[:u][end],(r[:nY],r[:nX]))
        heatmap!(ax,rotr90(uMat),colorrange=(-1.0, 1.0),colormap=:bwr)
        hidedecorations!(ax)
        hidespines!(ax)
        # Map parameters to axis in axes dictionary 
        axesDict[r[:λ]] = ax
        image!(ax,rotr90(maskImage))        
    end
    sortedλs  = unique!(sort(first.(keys(axesDict))))    
    for (i,λ) in enumerate(sortedλs)
        fig[2,i] = axesDict[λ]      
        Label(fig[2,i,Bottom()],"λ=$(λ)")  
    end
    Label(fig[2,1,Left()],"r=0.55\n ϕ₀=0.45")
    

    resize_to_layout!(fig)

    save(datadir("fromCSF","lambdaTest",mask,"$(mask)_lambdaTest.png"),fig)

    # display(fig)

end