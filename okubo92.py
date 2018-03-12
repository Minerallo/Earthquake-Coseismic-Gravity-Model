import sys
import numpy as np
from math import *


G = 6.67384e-11
nu = 0.25
beta = 0.309e-5
eps = 1e-14

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


ux = - U1 / (2 * np.pi) * chinnery(ux_ss, x, p, L, W, q, dip, nu) - \
    U2 / (2 * np.pi) * chinnery(ux_ds, x, p, L, W, q, dip, nu) + \
    U3 / (2 * np.pi) * chinnery(ux_tf, x, p, L, W, q, dip, nu)
uy = - U1 / (2 * np.pi) * chinnery(uy_ss, x, p, L, W, q, dip, nu) - \
    U2 / (2 * np.pi) * chinnery(uy_ds, x, p, L, W, q, dip, nu) + \
    U3 / (2 * np.pi) * chinnery(uy_tf, x, p, L, W, q, dip, nu)
uz = - U1 / (2 * np.pi) * chinnery(uz_ss, x, p, L, W, q, dip, nu) - \
    U2 / (2 * np.pi) * chinnery(uz_ds, x, p, L, W, q, dip, nu) + \
    U3 / (2 * np.pi) * chinnery(uz_tf, x, p, L, W, q, dip, nu)
ue = np.sin(strike) * ux - np.cos(strike) * uy
un = np.cos(strike) * ux + np.sin(strike) * uy
return ue, un, uz


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

# def ux_ss(xi, eta, q, dip, nu):
#     R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
#     u = xi * q / (R * (R + eta)) + \
#         I1(xi, eta, q, dip, nu, R) * np.sin(dip)
#     k = (q != 0)
#     # u[k] = u[k] + np.arctan2( xi[k] * (eta[k]) , (q[k] * (R[k])))
#     u[k] = u[k] + np.arctan((xi[k] * eta[k]) / (q[k] * R[k]))
#     return u


# dip-slip displacement subfunction [equation (59) p. 7139]
def Dh(xi, eta, q, dip, nu):
    R = sqrt(xi**2 + eta**2 + q**2)
    db = eta*sind(dip) - q*cosd(dip)
    u = - db.*q./(R.*(R + xi)) ...
        - sind(dip).*atan(xi.*eta./(q.*R)) ...
        + I5(xi, eta, q, dip, nu, R, db).*sind(dip).*cosd(dip)
    return u


def uy_ss(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    u = (eta * np.cos(dip) + q * np.sin(dip)) * q / (R * (R + eta)) + \
        q * np.cos(dip) / (R + eta) + \
        I2(eta, q, dip, nu, R) * np.sin(dip)
    return u


def uz_ss(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    db = eta * np.sin(dip) - q * np.cos(dip)
    u = (eta * np.sin(dip) - q * np.cos(dip)) * q / (R * (R + eta)) + \
        q * np.sin(dip) / (R + eta) + \
        I4(db, eta, q, dip, nu, R) * np.sin(dip)
    return u


def ux_ds(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    u = q / R - \
        I3(eta, q, dip, nu, R) * np.sin(dip) * np.cos(dip)
    return u


def uy_ds(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    u = ((eta * np.cos(dip) + q * np.sin(dip)) * q / (R * (R + xi)) -
         I1(xi, eta, q, dip, nu, R) * np.sin(dip) * np.cos(dip))
    k = (q != 0)
    u[k] = u[k] + np.cos(dip) * np.arctan((xi[k] * eta[k]) / (q[k] * R[k]))
    return u


def uz_ds(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    db = eta * np.sin(dip) - q * np.cos(dip)
    u = (db * q / (R * (R + xi)) -
         I5(xi, eta, q, dip, nu, R, db) * np.sin(dip) * np.cos(dip))
    k = (q != 0)
    # u[k] = u[k] + np.sin(dip) * np.arctan2(xi[k] * eta[k] , q[k] * R[k])
    u[k] = u[k] + np.sin(dip) * np.arctan((xi[k] * eta[k]) / (q[k] * R[k]))
    return u


def ux_tf(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    u = q ** 2 / (R * (R + eta)) - \
        I3(eta, q, dip, nu, R) * (np.sin(dip) ** 2)
    return u


def uy_tf(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    u = - (eta * np.sin(dip) - q * np.cos(dip)) * q / (R * (R + xi)) - \
        np.sin(dip) * xi * q / (R * (R + eta)) - \
        I1(xi, eta, q, dip, nu, R) * (np.sin(dip) ** 2)
    k = (q != 0)
    # u[k] = u[k] + np.sin(dip) * np.arctan2(xi[k] * eta[k] , q[k] * R[k])
    u[k] = u[k] + np.sin(dip) * np.arctan((xi[k] * eta[k]), (q[k] * R[k]))
    return u


def uz_tf(xi, eta, q, dip, nu):
    R = np.sqrt(xi**2 + eta**2 + q**2)
    db = eta * np.sin(dip) - q * np.cos(dip)
    u = (eta * np.cos(dip) + q * np.sin(dip)) * q / (R * (R + xi)) + \
        np.cos(dip) * xi * q / (R * (R + eta)) - \
        I5(xi, eta, q, dip, nu, R, db) * np.sin(dip)**2
    k = (q != 0)   # not at depth=0?
    u[k] = u[k] - np.cos(dip) * np.arctan((xi[k] * eta[k]) / (q[k] * R[k]))
    return u


def I1(xi, eta, q, dip, nu, R):
    db = eta * np.sin(dip) - q * np.cos(dip)
    if np.cos(dip) > eps:
        I = (1 - 2 * nu) * (- xi / (np.cos(dip) * (R + db))) - \
            np.sin(dip) / np.cos(dip) * \
            I5(xi, eta, q, dip, nu, R, db)
    else:
        I = -(1 - 2 * nu) / 2 * xi * q / (R + db) ** 2
    return I


def I2(eta, q, dip, nu, R):
    I = (1 - 2 * nu) * (-np.log(R + eta)) - \
        I3(eta, q, dip, nu, R)
    return I


def I3(eta, q, dip, nu, R):
    yb = eta * np.cos(dip) + q * np.sin(dip)
    db = eta * np.sin(dip) - q * np.cos(dip)
    if np.cos(dip) > eps:
        I = (1 - 2 * nu) * (yb / (np.cos(dip) * (R + db)) - np.log(R + eta)) + \
            np.sin(dip) / np.cos(dip) * \
            I4(db, eta, q, dip, nu, R)
    else:
        I = (1 - 2 * nu) / 2 * (eta / (R + db) + yb * q / (R + db) ** 2 - np.log(R + eta))
    return I


def I4(db, eta, q, dip, nu, R):
    if np.cos(dip) > eps:
        I = (1 - 2 * nu) * 1.0 / np.cos(dip) * \
            (np.log(R + db) - np.sin(dip) * np.log(R + eta))
    else:
        I = - (1 - 2 * nu) * q / (R + db)
    return I


def I5(xi, eta, q, dip, nu, R, db):
    X = np.sqrt(xi**2 + q**2)
    if np.cos(dip) > eps:
        I = (1 - 2 * nu) * 2 / np.cos(dip) * \
            np.arctan((eta * (X + q*np.cos(dip)) + X*(R + X) * np.sin(dip)) /
                      (xi*(R + X) * np.cos(dip)))
        I[xi == 0] = 0
    else:
        I = -(1 - 2 * nu) * xi * np.sin(dip) / (R + db)
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

    ux = - U1 / (2 * np.pi) * chinnery(ux_ss, x, p, L, W, q, dip, nu) - \
        U2 / (2 * np.pi) * chinnery(ux_ds, x, p, L, W, q, dip, nu) + \
        U3 / (2 * np.pi) * chinnery(ux_tf, x, p, L, W, q, dip, nu)
    uy = - U1 / (2 * np.pi) * chinnery(uy_ss, x, p, L, W, q, dip, nu) - \
        U2 / (2 * np.pi) * chinnery(uy_ds, x, p, L, W, q, dip, nu) + \
        U3 / (2 * np.pi) * chinnery(uy_tf, x, p, L, W, q, dip, nu)
    uz = - U1 / (2 * np.pi) * chinnery(uz_ss, x, p, L, W, q, dip, nu) - \
        U2 / (2 * np.pi) * chinnery(uz_ds, x, p, L, W, q, dip, nu) + \
        U3 / (2 * np.pi) * chinnery(uz_tf, x, p, L, W, q, dip, nu)
