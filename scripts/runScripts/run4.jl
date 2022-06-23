
using DrWatson
@quickactivate
using PhaseFieldCrystal

imagePath     = "data/exp_pro/croppedMask.png"
lX            = 200.0
r             = 0.5
ϕ0            = [-0.41,-0.42,-0.43,-0.44,-0.45,-0.46,-0.37]
a             = 2.0
δt            = 0.5
tMax          = 4000.0
loggerFlag    = 0
outputFlag    = 1
visualiseFlag = 0
nBlasThreads  = Threads.nthreads()

phaseFieldCrystal(imagePath,lX,r,ϕ0[1],a,δt,δt*5,0,0,0,nBlasThreads)
# phaseFieldCrystal(imagePath,lX,r,ϕ0[1],a,δt,tMax,loggerFlag,outputFlag,visualiseFlag,nBlasThreads)
# phaseFieldCrystal(imagePath,lX,r,ϕ0[2],a,δt,tMax,loggerFlag,outputFlag,visualiseFlag,nBlasThreads)
# phaseFieldCrystal(imagePath,lX,r,ϕ0[3],a,δt,tMax,loggerFlag,outputFlag,visualiseFlag,nBlasThreads)
phaseFieldCrystal(imagePath,lX,r,ϕ0[4],a,δt,tMax,loggerFlag,outputFlag,visualiseFlag,nBlasThreads)
# phaseFieldCrystal(imagePath,lX,r,ϕ0[5],a,δt,tMax,loggerFlag,outputFlag,visualiseFlag,nBlasThreads)
# phaseFieldCrystal(imagePath,lX,r,ϕ0[6],a,δt,tMax,loggerFlag,outputFlag,visualiseFlag,nBlasThreads)
