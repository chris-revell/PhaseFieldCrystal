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

function createRunDirectory(nGrid,lSpace,r,ϕ₀,a,δt,outInt,tMax,integrator)

    # Create directory for run data labelled with current time.
    folderDate = Dates.format(Dates.now(),"yyyy-mm-dd-HH-MM-SS")
    folderName = "$(pwd())/output/$(folderDate)"
    mkpath(folderName)

    # Store system parameters.
    open("$folderName/conditions.txt","w") do conditionsfile
        println(conditionsfile, "lSpace,     $lSpace    ")
        println(conditionsfile, "nGrid,      $nGrid     ")
        println(conditionsfile, "h,          $h         ")
        println(conditionsfile, "r,          $r         ")
        println(conditionsfile, "ϕ₀,         $ϕ₀        ")
        println(conditionsfile, "a,          $a         ")
        println(conditionsfile, "δt,         $δt        ")
        println(conditionsfile, "tMax,       $tMax      ")
        println(conditionsfile, "outInt,     $outInt    ")
        println(conditionsfile, "integrator, $integrator")
    end

    return folderName

end

export createRunDirectory

end
