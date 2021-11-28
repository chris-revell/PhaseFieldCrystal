using Images
using ImageBinarization
using FileIO
using ImageSegmentation
using ImageTransformations
using ImageView
using ImageSegmentation
using Random

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

#imagePath = "$(homedir())/Dropbox (The University of Manchester)/EM-images/mmp13ko3view/mmp13ko-3wiew_4800X_hui_0002.tif";
imagePath = "$(homedir())/Postdoc/Code/PhaseFieldCrystal/data/exp_raw/mmp13ko-3wiew_4800X_hui_0002.png";

image = load(imagePath)

binarizeFuns = [Intermodes(),
    Moments(),
    Otsu(),
    Polysegment(),
    UnimodalRosin(),
    #Sauvola(),
    #MinimumError(),
    #Niblack(),
    #Yen(),
    #Balanced(),
    #AdaptiveThreshold(image;window_size=50,percentage=10)
]

binarizeFunNames = ["Intermodes",
    "Moments",
    "Otsu",
    "Polysegment",
    "UnimodalRosin",
    #"Sauvola",
    #"MinimumError",
    #"Niblack()",
    #"Yen()",
    #"Balanced()",
    #"AdaptiveThreshold(image;window_size=50,percentage=10"
]

segmentFuns = [
    #"unseeded_region_growing",
    "felzenszwalb",
    #"meanshift",
    "fast_scanning"
]


distance = 10.0
filteredImage  = Gray.(image) #imfilter(Gray.(image),Kernel.gaussian(distance))
binarizedImage = copy(filteredImage)
segmentedImage = copy(image)

for (n,fun) in enumerate(binarizeFuns)

    binarizedImage .= binarize(filteredImage,fun)
    saveName = "$(homedir())/Postdoc/Code/PhaseFieldCrystal/data/exp_pro/binarizeFuns/$(binarizeFunNames[n]).png"
    save(saveName,binarizedImage)

    seg = felzenszwalb(binarizedImage, 300, 100)
    segmentedImage .= map(i->get_random_color(i), labels_map(seg))
    saveName = "$(homedir())/Postdoc/Code/PhaseFieldCrystal/data/exp_pro/binarizeFuns/$(binarizeFunNames[n])Felzenszwalb.png"
    save(saveName,segmentedImage)
    seg2 = prune_segments(seg, i->(segment_pixel_count(seg,i)<30000), (i,j)->(-segment_pixel_count(seg,j)))
    segmentedImage .= map(i->get_random_color(i), labels_map(seg2))
    saveName = "$(homedir())/Postdoc/Code/PhaseFieldCrystal/data/exp_pro/binarizeFuns/$(binarizeFunNames[n])FelzenszwalbPruned.png"
    save(saveName,segmentedImage)
    segmentedImage .= map(i->maskColour(i,seg2), labels_map(seg2))
    saveName = "$(homedir())/Postdoc/Code/PhaseFieldCrystal/data/exp_pro/binarizeFuns/$(binarizeFunNames[n])FelzenszwalbPrunedMask.png"
    save(saveName,segmentedImage)

    seg = fast_scanning(binarizedImage, 0.1)
    segmentedImage .= map(i->get_random_color(i), labels_map(seg))
    saveName = "$(homedir())/Postdoc/Code/PhaseFieldCrystal/data/exp_pro/binarizeFuns/$(binarizeFunNames[n])FastScanning.png"  # threshold = 0.1
    save(saveName,segmentedImage)
    seg2 = prune_segments(seg, i->(segment_pixel_count(seg,i)<30000), (i,j)->(-segment_pixel_count(seg,j)))
    segmentedImage .= map(i->get_random_color(i), labels_map(seg2))
    saveName = "$(homedir())/Postdoc/Code/PhaseFieldCrystal/data/exp_pro/binarizeFuns/$(binarizeFunNames[n])FastScanningPruned.png"
    save(saveName,segmentedImage)
    segmentedImage .= map(i->maskColour(i,seg2), labels_map(seg2))
    saveName = "$(homedir())/Postdoc/Code/PhaseFieldCrystal/data/exp_pro/binarizeFuns/$(binarizeFunNames[n])FastScanningPrunedMask.png"
    save(saveName,segmentedImage)
end
