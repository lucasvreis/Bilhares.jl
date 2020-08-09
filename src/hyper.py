from matplotlib.transforms import *
import numpy as np

class HyperTransf(Transform):
    def __init__(self, **kwargs):
        Transform.__init__(self, **kwargs)
    
    input_dims = 2
    output_dims = 2
    is_separable = False
    has_inverse = False

    def transform(self,values):
        ssq = np.diag(np.inner(values,values))
        k = 1 / (1 + np.sqrt(1 - ssq))
        return k * values

    def transform_path(self, path):
        return path.interpolated(10) # ??