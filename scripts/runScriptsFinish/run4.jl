
using DrWatson
@quickactivate
using PhaseFieldCrystal

imagePath     = "data/exp_pro/croppedMask.png"
lX            = 200.0
rs            = [0.30,0.35,0.40,0.45,0.50,0.55]
ms            = [0.1,0.2,0.3,0.4,0.5]
ϕ0s           = [-0.420,-0.425,-0.430,-0.435,-0.440,-0.445]
a             = 2.0
δt            = 0.5
tMax          = 10000.0
loggerFlag    = 0
outputFlag    = 1
visualiseFlag = 0
nBlasThreads  = Threads.nthreads()


# Warmup
phaseFieldCrystal(imagePath,lX,0.5,-0.4,0.1,a,δt,10.0,0,0,0,nBlasThreads;subFolder="")
r = rs[4]
for ϕ0 in ϕ0s[end-1:end]
    phaseFieldCrystal(imagePath,lX,r,ϕ0,0.7,a,δt,tMax,loggerFlag,outputFlag,visualiseFlag,nBlasThreads;subFolder="PhaseSpace2")
end
for m in ms
    for ϕ0 in ϕ0s
        phaseFieldCrystal(imagePath,lX,r,ϕ0,m,a,δt,tMax,loggerFlag,outputFlag,visualiseFlag,nBlasThreads;subFolder="PhaseSpace2")
    end
end
