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

using Laplacian

@inline @views function visualise(sol,ilow,ihigh,jlow,jhigh,N,h,folderName)

    # Plot phase field
    anim = @animate for u in sol.u
        heatmap(u[4:N+3,4:N+3],clims=(-1,1),aspect_ratio=:equal,border=:none,show=false,color=:hawaii)
    end every 5
    gif(anim,"output/$folderName/anim_u.gif",fps=2)

    # Plot 2nd derivative of phase field
    part1 = zeros(N+6,N+6)
    anim = @animate for u in sol.u
        ∇²!(part1, u, N, h, 0)
        heatmap(part1[4:N+3,4:N+3],aspect_ratio=:equal,border=:none,show=false,color=:roma)
    end every 5
    gif(anim,"output/$folderName/anim_delsquaredu.gif",fps=2)

    # Plot 4th derivative of phase field
    part2 = zeros(N+6,N+6)
    anim = @animate for u in sol.u
        ∇²!(part1, u, N, h, 0)
        ∇²!(part2, part1, N, h, 1)
        heatmap(part2[4:N+3,4:N+3],aspect_ratio=:equal,border=:none,show=false,color=:cork)
    end every 5
    gif(anim,"output/$folderName/anim_deltothefouru.gif",fps=2)

    return nothing

end

export visualise

end
