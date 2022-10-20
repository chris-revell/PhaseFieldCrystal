using DataFrames
using CSV

scalingLX = 200.0/1.82314

lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","lengthMeasurements.csv")))
croppedLX = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","croppedLX.csv")))
spacingData = DataFrame(CSV.File(datadir("exp_pro","emCentroidMeasurements","spacingData.csv")))

for r in eachrow(spacingData)
    croppedLXRow  = filter(:file => f->f==r[:file], croppedLX)
    display(croppedLXRow[1,:lX]*1000.0/r[:mean])
    display(croppedLXRow[1,:lX]*scalingLX)
end