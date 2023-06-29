using DrWatson; @quickactivate
using DataFrames
using CairoMakie
using Colors

runs = [f for f in readdir(datadir("fromCSF","lambdaTest")) if f!=".DS_Store"]

mOrder = [0.1,0.2,0.3,0.4]
λOrder = [10.0,30.0,50.0,70.0,90.0]

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

    fig = Figure(resolution=(2000,2000),fontsize=32)

    axesDict = Dict()
    # for (i,r) in enumerate(eachrow(subset(res, :m => m -> m .== 0.1)))
    for (i,r) in enumerate(eachrow(res))
        @show i
        ax = Axis(fig[findfirst(x->x==r[:m],mOrder),findfirst(x->x==r[:λ],λOrder)],aspect=DataAspect())
        uMat = reshape(r[:u][end],(r[:nY],r[:nX]))
        heatmap!(ax,rotr90(uMat),colorrange=(-1.0, 1.0),colormap=:bwr)
        hidedecorations!(ax)
        hidespines!(ax)
        # Map parameters to axis in axes dictionary 
        axesDict[(r[:m],r[:λ])] = ax
        image!(ax,rotr90(maskImage))
        Label(fig[findfirst(x->x==r[:m],mOrder),findfirst(x->x==r[:λ],λOrder),Bottom()],"m=$(r[:m]) λ=$(r[:λ])")
    end

    resize_to_layout!(fig)

    save(datadir("fromCSF","lambdaTest",mask,"$(mask)_lambda-m-phasespace.png"),fig)

    display(fig)

end