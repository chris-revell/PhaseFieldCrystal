push!(LOAD_PATH,"src/")
using DrWatson
@quickactivate "PhaseFieldCrystal"

println(
"""
Currently active project is: $(projectname())
Path of active project: $(projectdir())
"""
)
@info "Loading test parameters"
include("scripts/TestParameters.jl")
# @info "Precompiling PhaseFieldCrystal"
using PhaseFieldCrystal
