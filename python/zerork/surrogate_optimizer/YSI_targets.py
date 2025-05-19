
from .surrogate_mixture import TargetBase, RelativeErrorFn

class YSITarget(TargetBase):
    def __init__(self, target, weight, target_fn=RelativeErrorFn(), mode="mass"):
        assert(mode == "mass" or mode == "mole" or mode == "CT")
        self.mode = mode
        TargetBase.__init__(self,"YSI ({})".format(mode),target,weight,target_fn)
    def __call__(self, mixture):
        YSI = mixture.calcYSI(self.mode)
        error = TargetBase.__call__(self,YSI)
        return error
    def get_idt_sweeps(self):
        return None

