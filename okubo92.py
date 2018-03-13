import sys
import numpy as np
from math import *


G = 6.67384e-11
nu = 0.25
beta = 0.309e-5


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


def forward(x, y, xoff=0, yoff=0,
            depth=5e3, length=1e3, width=1e3,
            slip=0.0, opening=10.0,
            strike=0.0, dip=0.0, rake=0.0,
            nu=0.25):
    e = x
    n = y
    x = x - xoff
    y = y - yoff
    strike = np.deg2rad(strike)
    dip = np.deg2rad(dip)
    rake = np.de2rad(rake)
    L = lenght
    W = width
    U1 = cos(rake)*slip
    U2 = sin(rake)*slip
    U3 = opening

    d = depth + sin(dip)*W/2
    ec = e + cos(strike)*cos(dip)*W/2
    nc = n - sin(strike)*cos(dip)*W/2
    x = cos(strike)*nc + sin(strike)*ec + L/2
    y = sin(strike)*nc - cos(strike)*ec + cos(dip)*W

    p = y*cos(dip) + d*sin(dip)
    q = y*sin(dip) - d*cos(dip)

    def chinnery(f, x, p, L, W, q, dip, nu):
        # ''' % Chinnery's notation [equation (24) p. 1143]'''
        u = (f(x, p, q, dip, nu) -
             f(x, p - W, q, dip, nu) -
             f(x - L, p, q, dip, nu) +
             f(x - L, p - W, q, dip, nu))
        return u

    # strike-slip displacement subfunction [equation (58) p. 7139]
    def Sh(xi, eta, q, dip, nu):
        R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
        db = eta*sin(dip) - q*cos(dip)
        u = - db*q/(R*(R + eta)) - \
            q*sin(dip)/(R + eta) - \
            I4(db, eta, q, dip, nu, R)*sin(dip)
        return u

    # dip-slip displacement subfunction [equation (59) p. 7139]
    def Dh(xi, eta, q, dip, nu):
        R = sqrt(xi**2 + eta**2 + q**2)
        db = eta*sin(dip) - q*cos(dip)
        u = - db*q/(R*(R + xi)) - sin(dip)*np.arctan(xi*eta/(q*R)) + \
            I5(xi, eta, q, dip, nu, R, db)*sin(dip)*cos(dip)
        return u

    # tensile fault displacement subfunction [equation (60) p. 7139]
    def Th(xi, eta, q, dip, nu):
        R = sqrt(xi**2 + eta**2 + q**2)
        db = eta*sin(dip) - q*cos(dip)
        u = (eta*cos(dip) + q*sin(dip))*q/(R*(R + xi)) + \
            cos(dip)*(xi*q/(R*(R + eta)) -
                      np.arctan(xi*eta/(q*R))) - \
            I5(xi, eta, q, dip, nu, R, db)*sin(dip)**2
        return u

    #  I4 displacement subfunction [equations (61) and (63) p. 7139-7140]
    def I4(db, eta, q, dip, nu, R):
        if cos(dip) == 0:
            I = -(1 - 2*nu) * q/(R + db)
        else:
            I = (1 - 2*nu) * 1/cos(dip) * (log(R + db) - sin(dip)*log(R + eta))
        return I

    # I5 displacement subfunction [equations (62) and (64) p. 7140]
    def I5(xi, eta, q, dip, nu, R, db):
        if cos(dip) == 0:
            I = -(1 - 2*nu) * xi*sin(dip)/(R + db)
        else:
            I = (1 - 2*nu) * 2/cos(dip) * \
                np.arctan((-q*cos(dip) + (1 + sin(dip))*(R + eta))/(xi*cos(dip)))
        return I

    # Gravity subfunctions
    #  strike-slip gravity subfunction [equation (52) p. 7139]
    def Sg(xi, eta, q, dip, nu):
        R = sqrt(xi**2 + eta**2 + q**2)
        u = - q*sin(dip)/R + (q**2 * cos(dip))/(R*(R + eta))
        return u

    # dip-slip gravity subfunction [equation (53) p. 7139]
    def Dg(xi, eta, q, dip, nu):
        R = sqrt(xi**2 + eta**2 + q**2)
        db = eta*sin(dip) - q*cos(dip)
        u = 2*I2(xi, eta, q, R)*sin(dip) - q*db/(R*(R + xi))
        return u

    # tensile fault gravity subfunction [equation (54) p. 7139]
    def Tg(xi, eta, q, dip, nu):
        R = sqrt(xi**2 + eta**2 + q**2)
        yb = eta*cos(dip) + q*sin(dip)
        u = 2*I2(xi, eta, q, R)*cos(dip) + q*yb/(R*(R+xi)) + \
            q*xi*cos(dip)/(R*(R+eta))
        return u

    # cavity-filling gravity subfunction [equation (55) p. 7139]
    def Cg(xi, eta, q, dip, nu):
        R = sqrt(xi**2 + eta**2 + q**2)
        u = 2*I2(xi, eta, q, R)*cos(dip) - sin(dip)*log(R + xi)
        return u

    # I2 subfunction [equation (32) p. 7139]
    def I2(xi, eta, q, R):
        I = np.arctan((R + xi + eta)/q)
        return I

    # Elevation changes dH (must be computed first) [equation (57) p. 7139]
    dH = 1/(2*pi)*(U1*chinnery(Sh, x, p, L, W, q, dip, nu) +
                   U2*chinnery(Dh, x, p, L, W, q, dip, nu) +
                   U3*chinnery(Th, x, p, L, W, q, dip, nu))

    # Total gravity changes dG [equation (49) p. 7139]
    dG = rho*G*(U1*chinnery(Sg, x, p, L, W, q, dip, nu) +
                U2*chinnery(Dg, x, p, L, W, q, dip, nu) +
                U3*chinnery(Tg, x, p, L, W, q, dip, nu)) + \
        (rhop-rho)*G*U3*chinnery(Cg, x, p, L, W, q, dip, nu) - beta*dH


print(dG, dH)
