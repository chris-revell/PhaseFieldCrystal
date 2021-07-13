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

@inline @views function visualise(sol,N,folderName,imageMask)

    anim = @animate for i=1:(size(sol.t)[1])
        heatmap(sol.u[i][4:N+3,4:N+3].*imageMask,clims=(-1,1),aspect_ratio=:equal,border=:none,show=false,color=:hawaii)
    end every 5
    gif(anim,"output/$folderName/anim.gif",fps=2)

    return nothing

end

export visualise

end
