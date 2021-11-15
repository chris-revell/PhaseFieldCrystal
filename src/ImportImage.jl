#
#  ImportImage.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 07/07/2021.
#
#
# In place Laplacian function

module ImportImage

using Images

@inline function importImage(imagePath)

    grayImage = Gray.(load(imagePath))
    imageMask = Float64.(grayImage .> 0.5)
    nY = size(imageMask)[2]
    nX = size(imageMask)[1]

    return imageMask, nY, nX

end

export importImage

end
