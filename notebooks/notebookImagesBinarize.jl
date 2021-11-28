using Images
using ImageBinarization
using FileIO
using ImageSegmentation
using ImageTransformations
using ImageView
using ImageSegmentation
using Random

#imagePath = "$(homedir())/Dropbox (The University of Manchester)/EM-images/mmp13ko3view/mmp13ko-3wiew_4800X_hui_0002.tif";
imagePath = "$(homedir())/Postdoc/Code/PhaseFieldCrystal/data/exp_raw/mmp13ko-3wiew_4800X_hui_0002.png";

image = load(imagePath)

binarizeFuns = [Intermodes(),
    MinimumError(),
    Moments(),
    Otsu(),
    Polysegment(),
    UnimodalRosin(),
    Sauvola(),
    Niblack(),
    Yen(),
    Balanced(),
    AdaptiveThreshold(image;window_size=50,percentage=10)
]

binarizeFunNames = ["Intermodes()",
    "MinimumError()",
    "Moments()",
    "Otsu()",
    "Polysegment()",
    "UnimodalRosin()",
    "Sauvola()",
    "Niblack()",
    "Yen()",
    "Balanced()",
    "AdaptiveThreshold(image;window_size=50,percentage=10"
]


#distance = 10.0
#image .= imfilter(image2,Kernel.gaussian(distance))


for (n,fun) in enumerate(binarizeFuns)
    image2 = binarize(Gray.(image),fun)
    saveName = "$(homedir())/Postdoc/Code/PhaseFieldCrystal/data/exp_pro/binarizeFuns/$(binarizeFunNames[n]).png"
    save(saveName,image2)
end
