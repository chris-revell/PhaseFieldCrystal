#
#  Visualise.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 30/06/2021.
#
#
# Function to visualise results as animated gif; second function to import existing results from file

module Visualise

# Import Julia packages
using Plots
using ColorSchemes
using Printf
using JLD2
using UnPack
using DifferentialEquations

function visualise(sol,∇²,nGrid,freeEnergies,folderName)

    mat1 = zeros(nGrid^2)

    # # Plot phase field
    ENV["GKSwstype"]=100
    ENV["JULIA_GR_PROVIDER"] = "GR"
    anim = @animate for (i,u) in enumerate(sol.u)
        uInternal = reshape(u,(nGrid,nGrid))
        heatmap(uInternal,title="t=$(@sprintf("%.2f", sol.t[i]))",clims=(-1,1),aspect_ratio=:equal,border=:none,show=false,color=:hawaii)
    end
    #withenv("DYLD_LIBRARY_PATH"=>Plots.FFMPEG.FFMPEG_jll.LIBPATH[]) do
    gif(anim,"$folderName/anim_u.gif",fps=10)
    #end

    # # Plot 2nd derivative of phase field
    # anim2 = @animate for (i,u) in enumerate(sol.u)
    #     mat1 .= ∇²*u
    #     uInternal = reshape(mat1,(nGrid,nGrid))
    #     heatmap(uInternal,title="t=$(@sprintf("%.2f", sol.t[i]))",aspect_ratio=:equal,border=:none,show=false,color=:roma)
    # end
    # gif(anim2,"$folderName/anim_del2u.gif",fps=10)

    # # Plot 4th derivative of phase field
    # anim3 = @animate for (i,u) in enumerate(sol.u)
    #     mat1 .= ∇²*u
    #     mat1 .= ∇²*mat1
    #     uInternal = reshape(mat1,(nGrid,nGrid))
    #     heatmap(uInternal,title="t=$(@sprintf("%.2f", sol.t[i]))",aspect_ratio=:equal,border=:none,show=false,color=:cork)
    # end
    # gif(anim3,"$folderName/anim_del4u.gif",fps=10)

    plot(sol.t,freeEnergies)
    savefig("$folderName/freeEnergyVsTime.png")

    return nothing

end

# Function to import data needed for visualisation from file
# For example, run visualise(importData("output/folder/data.jld2")...)
function importData(path)

    data = load("$path/data.jld2")
    @unpack sol,∇²,nGrid,freeEnergies = data

    return sol,∇²,nGrid,freeEnergies,path

end

export visualise, importData

end
