#
#   Calculate H.O. approximations hw and b
#   Fox
#

import sys
import numpy as np


def oscb(A):
    hbarc = 197.327
    mp = 938.272
    mn = 939.565
    mass = 0.5*(mp+mn)
    #hw = 45.0*A**(-0.3) - 25.0*A**(-0.6)
    hw = 45.0*A**(-1/3) - 25.0*A**(-2/3)
    #hw = 41.0*A**(-1/3)

    b = hbarc / np.sqrt(hw*mass)
    return b


