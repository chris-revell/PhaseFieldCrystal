using DrWatson; @quickactivate
using DataFrames
using CairoMakie
using Colors
using FromFile
using GeometryBasics
using NumericalIntegration
using Printf

fig = Figure(fontsize=32,resolution=(800,500))
ax = Axis(fig[1,1])
ylims!(ax,(0,1))
r = "17tailT_4800X_HUI_0002_0"

function filterFunction(r,ϕ0)
    r==0.8 && ϕ0==0.43
end

maskIn = load(datadir("exp_pro","masksCompressed",r,"$r.png"))
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

results = collect_results(datadir("sims","timeResolution",r))
# results = filter([:r, :ϕ0] => filterFunction, results)

fig1 = Figure(figure_padding=0,resolution=(1000,1000),fontsize=64)
ax1 = CairoMakie.Axis(fig1[1,1],aspect=DataAspect())
nX = results[1,:nX]
nY = results[1,:nY]
uInternal = Observable(zeros(nX,nY))
heatmap!(ax1,uInternal,colorrange=(-1.0, 1.0),colormap=:bwr)
hidedecorations!(ax1)
hidespines!(ax1)
ax1.title = "t=0.0"
ax1.yreversed = true
resize_to_layout!(fig1)
image!(ax1,transpose(maskImage))

for (i,t) in enumerate(results[1,:t])
    uMat = reshape(results[1,:u][i],(nY,nX))
    ax1.title = "t=$(@sprintf("%.2f", t))"
    uInternal[] = transpose(uMat)
    uInternal[] = uInternal[]
    save(datadir("sims","timeResolution",r,"$i.png"),fig1)
end

