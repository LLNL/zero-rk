


import numpy
from .surrogate_mixture import TargetBase, NegativeDifferenceErrorFn

class SparsityTarget(TargetBase):
    def __init__(self, weight, target_count=5, mode="max", cutoff=1.0e-5):
        assert(mode in ["max", "sum"])
        name = "Sparsity"
        parent_target = 0.0
        self.mode = mode
        self.cutoff = cutoff
        self.target_count = target_count 
        TargetBase.__init__(self, name, parent_target, weight, NegativeDifferenceErrorFn())
    def __call__(self, mixture):
        def limiter(x):
            invx = 1/x
            invc = 1/self.cutoff
            return 1/(invx*invc/(invx+invc))
        vfs_active = [limiter(vf) for vf in mixture.getVolumeFractions() if vf > 0.0]
        if(self.mode == "max"):
            max_vf = numpy.sort(vfs_active)[-(self.target_count+1)]
            val = 1.0/max_vf
        elif(self.mode == "sum"):
            sum_vf = numpy.sum(numpy.sort(vfs_active)[:-self.target_count])
            val = sum_vf
        error = TargetBase.__call__(self,val)
        return error

