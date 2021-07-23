 #
#  Visualise.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 30/06/2021.
#
#
# Visualise results as animated gif

module Visualise

using Plots
using ColorSchemes
using Printf
using JLD2

using Laplacian

function visualise(sol,N,h,folderName,freeEnergies,imageMask)

    # Plot phase field
    anim = @animate for (i,u) in enumerate(sol.u)
        heatmap(u[4:N+3,4:N+3].*imageMask,title="t=$(@sprintf("%.2f", sol.t[i])), Free energy=$(@sprintf("%.2f", freeEnergies[i]))",clims=(-1,1),aspect_ratio=:equal,border=:none,show=false,color=:hawaii)
    end
    gif(anim,"output/$folderName/anim_u.gif",fps=5)

    # Plot 2nd derivative of phase field
    part1 = zeros(N+6,N+6)
    anim2 = @animate for (i,u) in enumerate(sol.u)
        ∇²!(part1, u, N, h, 0)
        heatmap(part1[4:N+3,4:N+3].*imageMask,title="t=$(@sprintf("%.2f", sol.t[i])), Free energy=$(@sprintf("%.2f", freeEnergies[i]))",aspect_ratio=:equal,border=:none,show=false,color=:roma)
    end
    gif(anim2,"output/$folderName/anim_del2u.gif",fps=5)

    # Plot 4th derivative of phase field
    part2 = zeros(N+6,N+6)
    anim3 = @animate for (i,u) in enumerate(sol.u)
        ∇²!(part1, u, N, h, 0)
        ∇²!(part2, part1, N, h, 1)
        heatmap(part2[4:N+3,4:N+3].*imageMask,title="t=$(@sprintf("%.2f", sol.t[i])), Free energy=$(@sprintf("%.2f", freeEnergies[i]))",aspect_ratio=:equal,border=:none,show=false,color=:cork)
    end
    gif(anim3,"output/$folderName/anim_del4u.gif",fps=5)

    plot(sol.t,freeEnergies)
    savefig("output/$folderName/freeEnergyVsTime.png")

    return nothing

end

# Function to import data needed for visualisation from file
# For example, run visualise(importData("ouput/folder/data.jld2")...)
function importData(path)

    data = load(path)

    sol          = data["sol"]
    N            = data["N"]
    h            = data["h"]
    folderName   = data["folderName"]
    freeEnergies = data["freeEnergies"]
    imageMask    = data["imageMask"]

    return sol,N,h,folderName,freeEnergies,imageMask

end

export visualise, importData

end
