using DrWatson 
@quickactivate
using PhaseFieldCrystal

for q2 in [0.5, 0.9, 1.1, 1.5]
    for ϕ0_2 in [0.30, 0.35, 0.40]
        for r2 in [0.5, 0.7, 0.9]
            phaseFieldCrystal(q2=q2, ϕ0_2=ϕ0_2, r2=r2, subFolder="Examples")
        end
    end
end
