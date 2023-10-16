using DrWatson
using JLD2
using CairoMakie
using UnPack


function filterNonFibril(x)
    if x<0.0
        return x
    else
        return 0.0
    end
end

# Scan data folder
runs = [f for f in readdir(datadir("sims", "examples")) if occursin(".jld2", f)]

# Figure canvas
fig = Figure(resolution=(1800,900))

# Dictionary to store parameter, axis pairs
axesDict = Dict()

# Create a new axis for each parameter set 
for r in runs
    ax = Axis(fig[1,1], aspect=DataAspect())
    hidedecorations!(ax); hidespines!(ax)
    ax.yreversed = true
    data = load(datadir("sims", "examples", r))
    @unpack ϕ0_1, ϕ0_2, r1, r2, q2, u, nX, nY = data
    uTmp = filterNonFibril.(transpose(reshape(u[end][1:nX*nY], (nY, nX)))) .- filterNonFibril.(transpose(reshape(u[end][1+nX*nY:end], (nY, nX))))
    heatmap!(ax, uTmp, colorrange=(-1.0, 1.0), colormap=(:bwr, 1.0))    
    # Store parameter set=>axis pair in dictionary
    axesDict[(ϕ0_2,r2,q2)] = ax
end

# Sort parameter set keys from dictionary
allKeys = keys(axesDict)
sortedKeys = sort(collect(allKeys), by=x->(x[3]*1000+x[2]*100+x[3]))

# Rearrange axes within figure according to sorted parameter sets 
for i=1:2
    for j=1:3
        for k=1:3
            fig[j,k+3*(i-1)] = axesDict[sortedKeys[(i-1)*9+(j-1)*3+k]]
            Label(
                fig[j,k+3*(i-1), Bottom()], 
                L" \phi_{0,2}=%$(sortedKeys[(i-1)*9+(j-1)*3+k][1]), r_{2}=%$(sortedKeys[(i-1)*9+(j-1)*3+k][2]), q_{2}=%$(sortedKeys[(i-1)*9+(j-1)*3+k][3])",
                )
        end
    end
end

display(fig)
save(datadir("sims", "examples", "bimodalPFCparameterSpace.png"), fig)