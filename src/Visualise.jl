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
using Plots; ENV["GKSwstype"]="nul" #ENV["GKSwstype"]=100 #ENV["JULIA_GR_PROVIDER"] = "GR"
using ColorSchemes
using Printf
using JLD2
using DifferentialEquations

# Import local modules


function visualise(sol,laplacianMatrix,mat1,N,h,freeEnergies,folderName)

    uInternal = zeros(N+6,N+6)

    # Plot phase field
    anim = @animate for (i,u) in enumerate(sol.u)
        uInternal .= reshape(u,(N+6,N+6))
        heatmap(uInternal[4:N+3,4:N+3],title="t=$(@sprintf("%.2f", sol.t[i])), Free energy=$(@sprintf("%.2f", freeEnergies[i]))",clims=(-1,1),aspect_ratio=:equal,border=:none,show=false,color=:hawaii)
    end
    gif(anim,"$folderName/anim_u.gif",fps=10)

    # Plot 2nd derivative of phase field

    anim2 = @animate for (i,u) in enumerate(sol.u)
        mat1 .= laplacianMatrix*u
        uInternal .= reshape(mat1,(N+6,N+6))
        heatmap(uInternal[4:N+3,4:N+3],title="t=$(@sprintf("%.2f", sol.t[i])), Free energy=$(@sprintf("%.2f", freeEnergies[i]))",aspect_ratio=:equal,border=:none,show=false,color=:roma)
    end
    gif(anim2,"$folderName/anim_del2u.gif",fps=10)

    # Plot 4th derivative of phase field

    anim3 = @animate for (i,u) in enumerate(sol.u)
        mat1 .= laplacianMatrix*u
        mat1 .= laplacianMatrix*mat1
        uInternal .= reshape(mat1,(N+6,N+6))
        heatmap(uInternal[4:N+3,4:N+3],title="t=$(@sprintf("%.2f", sol.t[i])), Free energy=$(@sprintf("%.2f", freeEnergies[i]))",aspect_ratio=:equal,border=:none,show=false,color=:cork)
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
    N            = data["N"]
    h            = data["h"]
    folderName   = data["folderName"]
    freeEnergies = data["freeEnergies"]
    imageMask    = data["imageMask"]

    return sol,N,h,freeEnergies,imageMask

end

export visualise, importData

end
