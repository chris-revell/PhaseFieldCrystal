
using DrWatson
@quickactivate
using PhaseFieldCrystal

imagePath     = "data/exp_pro/croppedMask.png"
lX            = 200.0
r             = 0.5
ϕ0            = [-0.37,-0.4,-0.45,-0.50,-0.55,-0.6]
a             = 2.0
δt            = 0.01
tMax          = 2000.0
loggerFlag    = 0
outputFlag    = 1
visualiseFlag = 0
nBlasThreads  = Threads.nthreads()

# phaseFieldCrystal(imagePath,lX,r,ϕ0[1],a,δt,tMax,loggerFlag,outputFlag,visualiseFlag,nBlasThreads)
phaseFieldCrystal(imagePath,lX,r,ϕ0[3],a,δt,tMax,loggerFlag,outputFlag,visualiseFlag,nBlasThreads)
# phaseFieldCrystal(imagePath,lX,r,ϕ0[3],a,δt,tMax,loggerFlag,outputFlag,visualiseFlag,nBlasThreads)
# phaseFieldCrystal(imagePath,lX,r,ϕ0[4],a,δt,tMax,loggerFlag,outputFlag,visualiseFlag,nBlasThreads)
# phaseFieldCrystal(imagePath,lX,r,ϕ0[5],a,δt,tMax,loggerFlag,outputFlag,visualiseFlag,nBlasThreads)
# phaseFieldCrystal(imagePath,lX,r,ϕ0[6],a,δt,tMax,loggerFlag,outputFlag,visualiseFlag,nBlasThreads)
