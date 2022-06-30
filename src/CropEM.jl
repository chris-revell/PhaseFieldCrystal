module CropEM

using Images
using ImageView
using Gtk.ShortNames
using DrWatson
using IntervalSets
using ImageBinarization
using ImageSegmentation
using Random

function maskColour(i,seg)
    if seg.segment_pixel_count[i]==maximum(values(seg.segment_pixel_count))
        return Gray{}(1)
    else
        return Gray{}(0)
    end
end

function get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
end

function cropEM(fileName)

    img = load(fileName)
    guidict = imshow(img)

    ranges = [1..2,1..2]
    # Create a condition object
    c = Condition()

    # Get the window
    win = guidict["gui"]["window"]

    # Notify the condition object when the window closes
    signal_connect(win, :destroy) do widget
        ranges[1] = (guidict["roi"]["zoomregion"][]).currentview.x
        ranges[2] = (guidict["roi"]["zoomregion"][]).currentview.y
        notify(c)
    end

    # Wait for the notification before proceeding ...
    wait(c)

    savename = datadir("exp_raw","cropped","cropped_$(splitpath(fileName)[end])")
    counter = 0
    nameNotFound = 1
    while nameNotFound==1
        savename = datadir("exp_raw","cropped","cropped_$(splitpath(fileName)[end][1:end-4])_$counter.png")
        if isfile(savename)
            counter+=1
        else
            nameNotFound=0
        end
    end

    croppedImage = img[ranges[2].left:ranges[2].right,ranges[1].left:ranges[1].right]

    save(savename,croppedImage)

    grayImage = Gray.(croppedImage)
    distance = 1.0
    filteredImage  = imfilter(grayImage,Kernel.gaussian(distance))
    binarizedImage = binarize(filteredImage,Intermodes())

    dilate!(binarizedImage)
    dilate!(binarizedImage)
    dilate!(binarizedImage)

    erode!(binarizedImage)
    erode!(binarizedImage)
    erode!(binarizedImage)

    seg1 = fast_scanning(binarizedImage, 0.01)

    seg2 = prune_segments(seg1, i->(segment_pixel_count(seg1,i)<1000), (i,j)->(segment_pixel_count(seg1,j)))
    vals = [seg2.segment_pixel_count[i] for i in keys(seg2.segment_pixel_count)]
    segmentLabelsOrderedBySize = [i for i in keys(seg2.segment_pixel_count)]
    segmentLabelsOrderedBySize .= segmentLabelsOrderedBySize[sortperm(vals)]
    seg3 = prune_segments(seg2, i->(segment_pixel_count(seg2,i)<segment_pixel_count(seg2,segmentLabelsOrderedBySize[end-1])), (i,j)->(-segment_pixel_count(seg2,j)))
    prunedImage = map(i->maskColour(i,seg3), labels_map(seg3))

    savename = datadir("exp_pro","mask","mask_$(splitpath(fileName)[end])")
    counter = 0
    nameNotFound = 1
    while nameNotFound==1
        savename = datadir("exp_pro","mask","mask_$(splitpath(fileName)[end][1:end-4])_$counter.png")
        if isfile(savename)
            counter+=1
        else
            nameNotFound=0
        end
    end

    # segmentedImage = map(i->get_random_color(i), labels_map(seg1))
    # guidict = imshow(segmentedImage)
    # c = Condition()
    # win = guidict["gui"]["window"]
    # signal_connect(win, :destroy) do widget
    #     notify(c)
    # end
    #
    # # Wait for the notification before proceeding ...
    # wait(c)

    save(savename,prunedImage)


end

# cropEM("data/exp_raw/png/16tailT_4800X_HUI_0001.png")

export cropEM

end
