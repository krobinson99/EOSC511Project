#!/usr/bin/env python
import numpy as np

def glacier(grid, ngrid, dt, T, small=False):  # return eta
    '''recommended values ngrid=11, dt=150, T=4*3600 (4 hours)???? CHANGE FOR OUR PROJECT
    '''
  
# Define our constants and dimensions
    g = 10 
    D = 200           #depth of our domain in x direction [m]
    L = 20e3          #length of our domain in x direction [m]
    dx = L/(ngridx-1)
    dz = D/(ngridz-1)
    C0 = 10           # input concentration of methane NOT TRUE
    S0 = 0            # Input concentration of salinity 
    zz = 2            # Depth of glacier outflow

# set up temporal scale T is total run time
    ntime = np.int(T/dt)

# initialize
    u, w, rho, S, C = initial(ngridx,ngridz)

    # main loop (leap-frog)
    for k in range(ntime):
        u, v, eta = stepper[grid](ngrid, f, g, Hu, Hv, dt, rdx, u, v, eta, up, vp, etap)

        # add forcing
        t = k * dt / 8640
        eta = eta + 0.1 * (1 - np.exp(-t*t)) * spatial

        # periodic boundary conditions
        u, v, eta = periodicbc(ngrid, u, v, eta)

        # exchange values
        u, v, eta, up, vp, etap = exchange(u, v, eta, up, vp, etap)

    # plot contours
    plotit(grid, ngrid, dx, x, y, u, v, eta)

    return


def initialize(ngrid):
    '''initialize a ngrid x ngrid domain, u, v,, all zero 
     we need density salinity, ch4''' 

    u = np.zeros((ngrid, ngrid))
    v = np.zeros_like(u)

    up = np.zeros_like(u)
    vp = np.zeros_like(u)

    return u, v, eta, up, vp, etap


def stepgrid2(ngrid, f, g, Hu, Hv, dt, rdx, u, v, eta, up, vp, etap):
    '''take a step forward using grid 2 (B)'''
    n = ngrid
    nm1 = ngrid-1
    nm2 = ngrid-2
    u[1:nm1, 1:nm1] = u[1:nm1, 1:nm1] + 2*dt * (f * vp[1:nm1, 1:nm1]   #u, v center grid points same as C and S
                                                 -g * (etap[2:n, 1:nm1] - etap[1:nm1, 1:nm1]
                                                        + etap[2:n, 2:n] - etap[1:nm1, 2:n]) * 0.5 * rdx)
    v[1:nm1, 1:nm1] = v[1:nm1, 1:nm1] + 2*dt * (-f*up[1:nm1, 1:nm1]
                                                 - g * (etap[1:nm1, 2:n] - etap[1:nm1, 1:nm1]
                                                      + etap[2:n, 2:n] - etap[2:n, 1:nm1]) * 0.5 * rdx)
    eta[1:nm1, 1:nm1] = eta[1:nm1, 1:nm1] + 2*dt * (-(Hu[1:nm1, 1:nm1] * up[1:nm1, 1:nm1] #eta grid points same as u and w
                                                      - Hu[0:nm2, 1:nm1] * up[0:nm2, 1:nm1]
                                                         + Hu[1:nm1, 0:nm2] * up[1:nm1, 0:nm2]
                                                         - Hu[0:nm2, 0:nm2] * up[0:nm2, 0:nm2]
                                                         + Hv[1:nm1, 1:nm1] * vp[1:nm1, 1:nm1]
                                                         - Hv[1:nm1, 0:nm2] * vp[1:nm1, 0:nm2]
                                                         + Hv[0:nm2, 1:nm1] * vp[0:nm2, 1:nm1]
                                                         - Hv[0:nm2, 0:nm2] * vp[0:nm2, 0:nm2]) * 0.5 * rdx)
    return u, v, eta


def boundary_steady(C, S):
    '''Sets the boundary conditions for the steady state and for the source and sink stage'''
    ## open water boundary
    C[-1, :] = 4.5 ## nM
    S[-1, :] = 35 ## PSU
    ## Glacier wall
    C[0, :] = C[1,:]
    C[0, zz] = C0
    S[0, :] = S[1,:]
    S[0, zz] = 0

    return C, S
