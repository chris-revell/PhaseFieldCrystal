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

@inline @views function visualise(sol,ilow,ihigh,jlow,jhigh,N,folderName)

    anim = @animate for i=1:(size(sol.t)[1])
        sol.u[i][ilow:ihigh,jlow:jhigh] .= 0.0
        heatmap(sol.u[i][4:N+3,4:N+3],clims=(-1,1),aspect_ratio=:equal,border=:none,show=false,color=:hawaii)
    end every 5
    gui()
    gif(anim,"output/$folderName/anim.gif",fps=2)
    display(anim)

    return nothing

end

export visualise

end
