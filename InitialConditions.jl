#
#  InitialConditions.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 13/05/2021.
#
#

module InitialConditions

@inline @views function initialConditions(N,L,α₀)

    # Gaussian initial conditions for order parameter field (u0)
    u0 = zeros(N+6,N+6)
    for i=1:N
        for j=1:N
            u0[i+3,j+3] = exp(-((L*(i-(N+1)/2)/(N-1))^2 + (L*(j-(N+1)/2)/(N-1))^2)/(2.0*5.0^2))
        end
    end

    # Set spatially varying diffusivity
    αᵢ = α₀.*ones(N+5,N+6)
    αⱼ = α₀.*ones(N+6,N+5)
    # Want to block off i = 11:20, j=101:150
    ilow = 135
    ihigh = 170
    jlow = 101
    jhigh = 150
    αᵢ[ilow-1:ihigh,jlow:jhigh] .= 0.0
    αⱼ[ilow:ihigh,jlow-1:jhigh] .= 0.0
    # Set smooth edges
    αᵢ[ilow-2,jlow-1:jhigh+1] .= 0.1
    αᵢ[ihigh+1,jlow-1:jhigh+1] .= 0.1
    αᵢ[ilow-1:ihigh+1,jlow-2] .= 0.1
    αᵢ[ilow-1:ihigh+1,jhigh+1] .= 0.1

    αᵢ[ilow-3,jlow-2:jhigh+2] .= 0.2
    αᵢ[ihigh+2,jlow-2:jhigh+2] .= 0.2
    αᵢ[ilow-2:ihigh+2,jlow-3] .= 0.2
    αᵢ[ilow-2:ihigh+2,jhigh+2] .= 0.2

    αᵢ[ilow-4,jlow-3:jhigh+3] .= 0.3
    αᵢ[ihigh+3,jlow-3:jhigh+3] .= 0.3
    αᵢ[ilow-3:ihigh+3,jlow-4] .= 0.3
    αᵢ[ilow-3:ihigh+3,jhigh+3] .= 0.3

    αᵢ[ilow-5,jlow-4:jhigh+4] .= 0.3
    αᵢ[ihigh+4,jlow-4:jhigh+4] .= 0.3
    αᵢ[ilow-4:ihigh+4,jlow-5] .= 0.3
    αᵢ[ilow-4:ihigh+4,jhigh+4] .= 0.3

    # Allocate additional arrays for later calculations
    deriv  = zeros(N+6,N+6)
    part1  = zeros(N+6,N+6)
    part2  = zeros(N+6,N+6)
    graduᵢ = zeros(N+5,N+6)
    graduⱼ = zeros(N+6,N+5)

return u0,deriv,part1,part2,αᵢ,αⱼ,graduᵢ,graduⱼ

end

export initialConditions

end
