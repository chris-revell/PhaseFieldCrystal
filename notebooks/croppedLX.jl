using DataFrames
using CSV

lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","lengthMeasurements.csv")))

cropped = [f for f in readdir(datadir("exp_pro","cropped")) if f[end-3:end]==".png"]

newMeasures = DataFrame(file = cropped, lX = zeros(length(cropped)))

for (n,row) in enumerate(eachrow(newMeasures))
    originalData = filter(:File => F->F=="$(row[:file][1:end-6]).png", lengthMeasurements)
    μmPerPixel = originalData[!,:length]/originalData[!,:Pixels]
    
    croppedImage = load(datadir("exp_pro","cropped",row[:file]))
    imSize       = size(croppedImage)
    
    # display(row[:file])
    # display(imSize)

    imageWidth = μmPerPixel*imSize[2]
    
    newMeasures[n,:lX] = imageWidth[1]

end

writedlm(datadir("exp_pro","lengthMeasurements","croppedLX.csv"),eachrow(newMeasures),",")