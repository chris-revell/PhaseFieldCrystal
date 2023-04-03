#
#  ArrayLoop.jl
#  PhaseFieldCrystal
#
#  Created by Christopher Revell on 15/06/2022.
#
#
#

module ArrayLoop

arrayLoop(a,nGrid) = mod(nGrid+a-1,nGrid)+1

export arrayLoop

end
