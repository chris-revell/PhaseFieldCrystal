using DrWatson
using FromFile
@quickactivate "PhaseFieldCrystal"
@info "Loading test parameters"
include("scripts/TestParameters.jl")
@info "Precompiling PhaseFieldCrystal"
@from "src/PhaseFieldCrystal.jl" using PhaseFieldCrystal
