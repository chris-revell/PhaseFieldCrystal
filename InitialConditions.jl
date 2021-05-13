#
#  InitialConditions.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 13/05/2021.
#
#

module InitialConditions

@inline @views function initialConditions(N,L)

    # Gaussian initial conditions 
    u0 = zeros(N+6,N+6)
    for i=1:N
        for j=1:N
            u0[i+3,j+3] = exp(-((L*(i-(N+1)/2)/(N-1))^2 + (L*(j-(N+1)/2)/(N-1))^2)/(2.0*5.0^2))
        end
    end

    deriv = zeros(N+6,N+6)
    part1 = zeros(N+6,N+6)
    part2 = zeros(N+6,N+6)

return u0,deriv,part1,part2

end

export initialConditions

end
