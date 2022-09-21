#
#  ColourFunctions.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 14/09/2022.
#
#
# 

module ColourFunctions

using Random
using Images
using Colors
using ImageSegmentation

# Function to set random colour for each segment
function get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
end

# Function to set white for the largest segment and black for others
function maskColour(i,seg)
    if seg.segment_pixel_count[i]==maximum(values(seg.segment_pixel_count))
        return Gray{}(1)
    else
        return Gray{}(0)
    end
end

export get_random_color, maskColour

end