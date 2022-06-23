module CropEM

using Images
using ImageView
using Gtk.ShortNames
using DrWatson
using IntervalSets

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

    save(savename,img[ranges[2].left:ranges[2].right,ranges[1].left:ranges[1].right])

end

# cropEM("data/exp_raw/png/16tailT_4800X_HUI_0001.png")

export cropEM

end
