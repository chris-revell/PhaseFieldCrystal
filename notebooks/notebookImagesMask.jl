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
function final(i)
    if i==2
        return RGB{}(0,0,0)
    else
        return RGB{}(1,1,1)
    end
end

#imagePath = "$(homedir())/Dropbox (The University of Manchester)/EM-images/mmp13ko3view/mmp13ko-3wiew_4800X_hui_0002.tif";
imagePath = "$(homedir())/Postdoc/Code/PhaseFieldCrystal/data/exp_raw/mmp13ko-3wiew_4800X_hui_0002.png";

image = load(imagePath)

binarizeFun = Intermodes()
#binarizeFun = AdaptiveThreshold(image;window_size=50,percentage=10);

image2 = copy(image)

image3 = binarize(Gray.(image2),binarizeFun)
# saveName = "$(homedir())/Postdoc/Code/PhaseFieldCrystal/data/exp_pro/MaskDist=$distance$(split(split(imagePath,"/")[end],".")[1]).png"
# save(saveName,image2)


#seg = unseeded_region_growing(image3, 0.01)
#seg = felzenszwalb(image3, 300, 100)
#seg = meanshift(image3, 16, 8/255)
seg = fast_scanning(image3, 0.1)  # threshold = 0.1

segmentedImage = map(i->get_random_color(i), labels_map(seg))
# saveName = "$(homedir())/Postdoc/Code/PhaseFieldCrystal/data/exp_pro/Segmented$(split(split(imagePath,"/")[end],".")[1]).png"
# save(saveName,segmentedImage)

seg2 = prune_segments(seg, i->(segment_pixel_count(seg,i)<10000), (i,j)->(-segment_pixel_count(seg,j)))

# segmentedImage2 = map(i->get_random_color(i), labels_map(seg2))
# saveName = "$(homedir())/Postdoc/Code/PhaseFieldCrystal/data/exp_pro/SegmentedPruned$(split(split(imagePath,"/")[end],".")[1]).png"
# save(saveName,segmentedImage2)

segmentedImage3 = map(i->final(i), labels_map(seg2))
saveName = "$(homedir())/Postdoc/Code/PhaseFieldCrystal/data/exp_pro/Final$(split(split(imagePath,"/")[end],".")[1]).png"
save(saveName,segmentedImage3)

display(segmentedImage3)
