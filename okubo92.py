import sys
import numpy as np

e = 1
n = 2
depth = 3
strike = 4
dip = 5
L = 6
W = 7
rake = 8
slip = 9.5
U3 = 10
rho = 11
rhop = 12.2


def Parameters(*args):
    if len(sys.argv) < 11 or len(sys.argv) > 14:
        print('wrong number of input arguments')
    # if len(sys.argv) == 11:  # input argument number
    #     sys.argv.append[11] = rho  # too see later
    #     print('test1')
    # elif len(sys.argv) == 12:
    #     sys.argv[11] = rhop
    #     print('test')
    for arg in args:
        print('another arg :', arg)


Parameters(e, n, depth, strike, dip, L, W,
           rake, slip, U3, rho, rhop)

#  Newton's gravitational constant
G = 6.67384e-11
# (m^3/kg/s^2)

#  Default values for optional input arguments
nu = 0.25
# isotropic Poisson's ratio
beta = 0.309e-5
# free-air gravity gradient (m.s^-2/m)

# args = float(*args)
# e = args[1]
# n = args[2]
# depth = args[3]
# strike = args[4]
# dip = args[5]
# L = args[6]
# W = args[7]
# rake = args[8]
# slip = args[9]
# U3 = args[10]
# rho = args[11]
