using Images
using ImageBinarization
using FileIO
using ImageSegmentation
using ImageTransformations
using ImageView
using ImageSegmentation
using Random
using Base.Filesystem

function get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
end

function maskColour(i,seg)
    if seg.segment_pixel_count[i]==maximum(values(seg.segment_pixel_count))
        return RGB{}(0,0,0)
    else
        return RGB{}(1,1,1)
    end
end


fileDirectory = "/Users/christopher/Dropbox (The University of Manchester)/EM-images/mmp13ko3view"


binarizeFuns = [Intermodes(),
    Otsu(),
    Polysegment(),
]

binarizeFunNames = ["Intermodes",
    "Otsu",
    "Polysegment",
]


distance = 2.0

files = [f for f in readdir(fileDirectory) if f[end-2:end]=="png"]

for f in files
    mkpath("$fileDirectory/$(f[1:end-4])")

    image = load("$fileDirectory/$f")
    filteredImage  = imfilter(Gray.(image),Kernel.gaussian(distance))
    binarizedImage = copy(filteredImage)
    segmentedImage = copy(image)

    for (n,fun) in enumerate(binarizeFuns)

        binarizedImage .= binarize(filteredImage,fun)
        saveName = "$fileDirectory/$(f[1:end-4])/$(binarizeFunNames[n]).png"
        save(saveName,binarizedImage)

        seg = fast_scanning(binarizedImage, 0.1)
        segmentedImage .= map(i->get_random_color(i), labels_map(seg))
        saveName = "$fileDirectory/$(f[1:end-4])/$(binarizeFunNames[n])FastScanning.png"  # threshold = 0.1
        save(saveName,segmentedImage)
        seg2 = prune_segments(seg, i->(segment_pixel_count(seg,i)<30000), (i,j)->(-segment_pixel_count(seg,j)))
        # segmentedImage .= map(i->get_random_color(i), labels_map(seg2))
        # saveName = "$fileDirectory/$(f[1:end-4])/$(binarizeFunNames[n])FastScanningPruned.png"
        # save(saveName,segmentedImage)
        segmentedImage .= map(i->maskColour(i,seg2), labels_map(seg2))
        saveName = "$fileDirectory/$(f[1:end-4])/$(binarizeFunNames[n])FastScanningPrunedMask.png"
        save(saveName,segmentedImage)
    end
end
