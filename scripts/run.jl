
using DrWatson 
@quickactivate
using Base.Threads
using PhaseFieldCrystal
using CSV
using DataFrames
using DelimitedFiles

okMasks = Vector(readdlm(datadir("exp_pro","filesToUse.txt"))[:,1])
croppedLX = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","croppedLX.csv")))
EMspacingData = DataFrame(CSV.File(datadir("exp_pro","emCentroidMeasurements","EMspacingData.csv")))

c = Dict()

c[:imagelXpairs] = Tuple[]
for m in okMasks
    filteredCroppedLX = filter(:file => f->f==m, croppedLX)[1,:lX]
    filteredEMspacingData = filter(:file => f->f==m, EMspacingData)[1,:mean]
    # scalingLX = (fibril spacing in q units)/(fibril spacing in nm)
    scalingLX = 2π/filteredEMspacingData
    push!(c[:imagelXpairs],(m,filteredCroppedLX*1000*scalingLX))
end
c[:r]             = [0.80]
c[:m]             = [0.1]
c[:ϕ0]            = [0.40]
c[:λ]             = [10.0]
c[:a]             = [2.0]
c[:δt]            = [0.1]
c[:tMax]          = [1000.0]
c[:loggerFlag]    = [0]
c[:outCount]      = [1000]
c[:outputFlag]    = [1]
c[:visualiseFlag] = [0]
c[:freeEnergyFlag]= [1]
c[:nBlasThreads]  = [1]
c[:subFolderName] = "timeResolution"

dl = dict_list(c)

ps = dl[2]

# Warmup
phaseFieldCrystal("data/exp_pro/testMask.png",200.0,0.7,-0.41,0.1,2.0,1.0,0.5,10.0,1,0,0,0,0,1)
# Run 
phaseFieldCrystal(datadir("exp_pro","masksCompressed",(ps[:imagelXpairs][1])[1:end-4],ps[:imagelXpairs][1]),ps[:imagelXpairs][2],ps[:r],ps[:ϕ0],ps[:m],ps[:a],ps[:λ],ps[:δt],ps[:tMax],ps[:outCount],ps[:loggerFlag],ps[:outputFlag],1,ps[:freeEnergyFlag],ps[:nBlasThreads];subFolder=ps[:subFolderName])
