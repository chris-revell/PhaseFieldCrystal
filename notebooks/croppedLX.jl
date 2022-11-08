using DataFrames
using CSV
using DelimitedFiles

lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","lengthMeasurements.csv")))

cropped = [f for f in readdir(datadir("exp_pro","cropped")) if occursin("endo",f)]#f[end-3:end]==".png"]

newMeasures = DataFrame(file = cropped, lX = zeros(length(cropped)))

for (n,row) in enumerate(eachrow(newMeasures))
    originalData = filter(:File => F->occursin(row[:file][1:end-6],F), lengthMeasurements)
    μmPerPixel = originalData[!,:length]/originalData[!,:Pixels]
    croppedImage = load(datadir("exp_pro","cropped",row[:file]))
    imSize       = size(croppedImage)
    imageWidth = μmPerPixel*imSize[2]
    newMeasures[n,:lX] = imageWidth[1]
end

# writedlm(datadir("exp_pro","lengthMeasurements","croppedLX.csv"),eachrow(newMeasures),",")
open(datadir("exp_pro","lengthMeasurements","croppedLX.csv"), "a") do io
    writedlm(io, eachrow(newMeasures),",")
end