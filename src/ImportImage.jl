#
#  ImportImage.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 07/07/2021.
#
#
# Import electron micrograph or test image as a matrix of 1s and 0s to map inter-cellular spaces

# module ImportImage

using Images

@inline function importImage(imagePath)

    grayImage = Gray.(load(imagePath))
    imageMask = Float64.(grayImage .> 0.5)
    nY = size(imageMask)[1]
    nX = size(imageMask)[2]

    return imageMask, nX, nY

end

# export importImage

# end
