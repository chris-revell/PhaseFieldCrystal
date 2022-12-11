using DrWatson; @quickactivate
using DataFrames
using CairoMakie
using Colors
using FromFile
using GeometryBasics
using NumericalIntegration

fig = Figure(fontsize=32,resolution=(800,500))
ax = Axis(fig[1,1])
ylims!(ax,(0,1))
r = "17tailT_4800X_HUI_0012_0"

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

results = collect_results(datadir("fromCSF","allMasksPhasespaceSeparateLengths",r))
subsetResults = filter([:r, :ϕ0] => filterFunction, results)
points = Point2[]
# Grid for numerical integration 
h = 1.0#subsetResults[1,:h]
nY = subsetResults[1,:nY]
nX = subsetResults[1,:nX]
ptsX = range(0, stop=h*nX, length=nX)
ptsY = range(0, stop=h*nY, length=nY)
for (i,t) in enumerate(subsetResults[1,:t])
    uMat = reshape(subsetResults[1,:u][i],(nY,nX))
    ϕSum = integrate((ptsY,ptsX),uMat)
    aSum = integrate((ptsY,ptsX),Float64.(maskIn))
    push!(points,Point2(t,ϕSum/aSum))
end
lines!(ax,points)


function filterFunction(r,ϕ0)
    r==0.8 && ϕ0==0.4
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

subsetResults = filter([:r, :ϕ0] => filterFunction, results)
points = Point2[]
# Grid for numerical integration 
h = 1.0#subsetResults[1,:h]
nY = subsetResults[1,:nY]
nX = subsetResults[1,:nX]
ptsX = range(0, stop=h*nX, length=nX)
ptsY = range(0, stop=h*nY, length=nY)
for (i,t) in enumerate(subsetResults[1,:t])
    uMat = reshape(subsetResults[1,:u][i],(nY,nX))
    ϕSum = integrate((ptsY,ptsX),uMat)
    aSum = integrate((ptsY,ptsX),Float64.(maskIn))
    push!(points,Point2(t,ϕSum/aSum))
end
lines!(ax,points)

function filterFunction(r,ϕ0)
    r==0.8 && ϕ0==0.45
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

subsetResults = filter([:r, :ϕ0] => filterFunction, results)
points = Point2[]
# Grid for numerical integration 
h = 1.0#subsetResults[1,:h]
nY = subsetResults[1,:nY]
nX = subsetResults[1,:nX]
ptsX = range(0, stop=h*nX, length=nX)
ptsY = range(0, stop=h*nY, length=nY)
for (i,t) in enumerate(subsetResults[1,:t])
    uMat = reshape(subsetResults[1,:u][i],(nY,nX))
    ϕSum = integrate((ptsY,ptsX),uMat)
    aSum = integrate((ptsY,ptsX),Float64.(maskIn))
    push!(points,Point2(t,ϕSum/aSum))
end
lines!(ax,points)





ylims!(ax,(0,0.6))
xlims!(ax,(0,1000))

ax.xlabel = "Time"
ax.ylabel = L"ϕ̄"
display(fig)
save(datadir("fromCSF","phiConservation.png"),fig)