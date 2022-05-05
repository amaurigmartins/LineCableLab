#  -*- coding: utf-8 -*-
"""

Author: Amauri Martins-Britto
Date: 09/12/2021

"""

import sys

from mpmath import *


mp.dps = 25


def foo(zin_re, zin_im, param_re, param_im, param2_re, param2_im):
    Zin_tmp1 = complex(float(zin_re), float(zin_im))
    param = complex(float(param_re), float(param_im))
    param2 = complex(float(param2_re), float(param2_im))
    Zin_tmp2 = besseli(0, param) * besselk(1, param2) + besselk(0, param) * besseli(1, param2)
    D = besseli(1,param)*besselk(1,param2)-besselk(1,param)*besseli(1,param2)
    zin=Zin_tmp1/D*Zin_tmp2
    return zin


if __name__ == '__main__':

    zin_re, zin_im, param_re, param_im, param2_re, param2_im = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]

    print(foo(zin_re, zin_im, param_re, param_im, param2_re, param2_im))
