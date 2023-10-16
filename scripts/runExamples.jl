using DrWatson 
@quickactivate
using PhaseFieldCrystal

for q2 in [0.5, 1.5]
    for ϕ0_2 in [0.41, 0.43, 0.45]
        for r2 in [0.5, 0.7, 0.9]
            phaseFieldCrystal(q2=q2, ϕ0_2=ϕ0_2, r2=r2, subFolder="Examples",tMax=500.0)
        end
    end
end
