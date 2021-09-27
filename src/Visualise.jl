 #
#  Visualise.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 30/06/2021.
#
#
# Visualise results as animated gif

module Visualise

# Import Julia packages
using Plots
#ENV["GKSwstype"]="nul"
ENV["GKSwstype"]=100
ENV["JULIA_GR_PROVIDER"] = "GR"
using ColorSchemes
using Printf
using JLD2
using DifferentialEquations

# Import local modules


function visualise(sol,∇²,mat1,nGrid,nGhosts,h,freeEnergies,folderName)

    # Plot phase field
    anim = @animate for (i,u) in enumerate(sol.u)
        uInternal = reshape(u,(nGrid+nGhosts,nGrid+nGhosts))
        heatmap(uInternal[4:nGrid+3,4:nGrid+3],title="t=$(@sprintf("%.2f", sol.t[i]))",clims=(-1,1),aspect_ratio=:equal,border=:none,show=false,color=:hawaii)
    end
    gif(anim,"$folderName/anim_u.gif",fps=10)

    # Plot 2nd derivative of phase field

    anim2 = @animate for (i,u) in enumerate(sol.u)
        mat1 .= ∇²*u
        uInternal = reshape(mat1,(nGrid+nGhosts,nGrid+nGhosts))
        heatmap(uInternal[4:nGrid+3,4:nGrid+3],title="t=$(@sprintf("%.2f", sol.t[i]))",aspect_ratio=:equal,border=:none,show=false,color=:roma)
    end
    gif(anim2,"$folderName/anim_del2u.gif",fps=10)

    # Plot 4th derivative of phase field

    anim3 = @animate for (i,u) in enumerate(sol.u)
        mat1 .= ∇²*u
        mat1 .= ∇²*mat1
        uInternal = reshape(mat1,(nGrid+nGhosts,nGrid+nGhosts))
        heatmap(uInternal[4:nGrid+3,4:nGrid+3],title="t=$(@sprintf("%.2f", sol.t[i]))",aspect_ratio=:equal,border=:none,show=false,color=:cork)
    end
    gif(anim3,"$folderName/anim_del4u.gif",fps=10)

    plot(sol.t,freeEnergies)
    savefig("$folderName/freeEnergyVsTime.png")

    return nothing

end

# Function to import data needed for visualisation from file
# For example, run visualise(importData("output/folder/data.jld2")...,output/folder/)
function importData(path)

    data = load(path)

    sol          = data["sol"]
    nGrid        = data["nGrid"]
    
    h            = data["h"]
    folderName   = data["folderName"]
    freeEnergies = data["freeEnergies"]
    imageMask    = data["imageMask"]

    return sol,nGrid,h,freeEnergies,imageMask

end

export visualise, importData

end
