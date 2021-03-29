#
#  CreateRunDirectory.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 22/03/2021.
#
#
# Function to create a directory in which to store run results and parameters, with directory name given by currentdate and time

module CreateRunDirectory

# Julia packages
using Dates
using Base.Filesystem
using DelimitedFiles

function createRunDirectory(L,h,α,γ,outInt,tMax)

    # Create directory for run data labelled with current time.
    foldername = Dates.format(Dates.now(),"yyyy-mm-dd-HH-MM-SS")
    mkpath("output/$(foldername)")

    # Store system parameters.
    open("output/$(foldername)/conditions.txt","w") do conditionsfile
        println(conditionsfile, "L,      $L     ")
        println(conditionsfile, "h,     $h    ")
        println(conditionsfile, "α,      $α     ")
        println(conditionsfile, "γ,      $γ     ")
        println(conditionsfile, "tMax,   $tMax  ")
        println(conditionsfile, "outInt, $outInt")
    end

    return foldername

end

export createRunDirectory

end
