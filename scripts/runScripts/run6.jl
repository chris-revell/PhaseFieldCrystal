
using DrWatson
@quickactivate
using PhaseFieldCrystal

imagePath     = "data/exp_pro/croppedMask.png"
lX            = 200.0
rs            = [0.5,0.6,0.7,0.8,0.9,1.0]
ms            = [0.1,0.3,0.5,0.7,0.9,1.1]
ϕ0s           = [-0.36,-0.38,-0.39,-0.40,-0.42,-0.44]
a             = 2.0
δt            = 0.5
tMax          = 20000.0
loggerFlag    = 0
outputFlag    = 1
visualiseFlag = 0
nBlasThreads  = Threads.nthreads()


# Warmup
phaseFieldCrystal(imagePath,lX,0.5,-0.4,0.1,a,δt,10.0,0,0,0,nBlasThreads;subFolder="")
r = rs[6]
for m in ms
    for ϕ0 in ϕ0s
        phaseFieldCrystal(imagePath,lX,r,ϕ0,m,a,δt,tMax,loggerFlag,outputFlag,visualiseFlag,nBlasThreads;subFolder="PhaseSpace")
    end
end
