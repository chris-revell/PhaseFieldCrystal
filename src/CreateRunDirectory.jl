#
#  CreateRunDirectory.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 22/03/2021.
#
#
# Function to create a directory in which to store run results and parameters, with directory name given by currentdate and time

module CreateRunDirectory

# Import Julia packages
using Dates
using Base.Filesystem
using DelimitedFiles

# Import local modules
# include("<Module>.jl"); using .Module

function createRunDirectory(L,N,h,r,ϕ₀,α,q,outInt,tMax)

    # Create directory for run data labelled with current time.
    foldername = Dates.format(Dates.now(),"yyyy-mm-dd-HH-MM-SS")
    mkpath("output/$(foldername)")

    # Store system parameters.
    open("output/$(foldername)/conditions.txt","w") do conditionsfile
        println(conditionsfile, "L,      $L     ")
        println(conditionsfile, "N,      $N     ")
        println(conditionsfile, "h,      $h     ")
        println(conditionsfile, "r,      $r     ")
        println(conditionsfile, "ϕ₀,     $ϕ₀    ")
        println(conditionsfile, "α,      $α     ")
        println(conditionsfile, "q,      $q     ")
        println(conditionsfile, "tMax,   $tMax  ")
        println(conditionsfile, "outInt, $outInt")
    end

    return foldername

end

export createRunDirectory

end
