#!/usr/bin/env python
import numpy as np

def glacier(grid, ngrid, dt, T, small=False):  # return eta
    '''recommended values ngrid=11, dt=150, T=4*3600 (4 hours)???? CHANGE FOR OUR PROJECT
    '''
  
# Define our constants and dimensions
    g = 10 
    D = 200 #depth of our domain in x direction [m]
    L = 20e3 #length of our domain in x direction [m]
    dx = L/(ngridx-1)
    dz = D/(ngridz-1)

# write out dx/dt and sqrt(g*H) for comparison
    print('dx/dt {0:.3f}\n'.format(dx/dt))
    print('sqrt(g*H0) {0:.3f}\n'.format(np.sqrt(g*H0)))
    rdx = 1./dx

# set up Gaussian for forcing (its a bump in the center of the domain)
    x, y, spatial = findforcing(L, dx, ngrid)

# set up temporal scale T is total run time
    ntime = np.int(T/dt)

# initialize
    u, v, eta, up, vp, etap = initial(ngrid)

# set-up topography
    Hu, Hv = find_depth[grid](H0, ngrid)

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

def findforcing(L, dx, ngrid):
    '''define a two dimensional Gaussian of forcing with half-width L'''
    x = np.arange(0, L + 0.1*dx, dx)
    y = x  # symmetrical case
    spatial = np.zeros((ngrid, ngrid))

    for i in range(ngrid):
        for j in range(ngrid):
            spatial[i,j] = np.exp(-((y[j] - L/2.) * (y[j] - L/2)
                                    +(x[i] - L/2.) * (x[i] - L/2))
                                    /(L*L / 48.))
    return x, y, spatial


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


def periodicbc(ngrid, u, v, eta):
    '''do periodic boundary conditions'''
    eta[0, :] = eta[-2, :]
    eta[-1, :] = eta[1, :]
    eta[:, 0] = eta[:, -1]
    eta[:, -1] = eta[:, 1]

    u[0, :] = u[-2, :]
    u[-1, :] = u[1, :]
    u[:, 0] = u[:, -1]
    u[:, -1] = u[:, 1]

    v[0, :] = v[-2, :]
    v[-1, :] = v[1, :]
    v[:, 0] = v[:, -1]
    v[:, -1] = v[:, 1]

    return u, v, eta

def exchange(u, v, eta, up, vp, etap):
    '''swap new and old values'''
    store = np.zeros_like(vp)  # make sure store is its own array, not just another name for vp
    store = vp
    vp = v
    v = store
    store = up
    up = u
    u = store
    store = etap
    etap = eta
    eta = store

    return u, v, eta, up, vp, etap

def plotit(grid, ngrid, dx, x, y, u, v, eta):
    '''Contour plots of u, v, eta'''

    shift = {1: (0, 0), 2: (0.5, 0.5), 3: (0.5, 0)}

    fig, ax = plt.subplots(2,2, figsize=(10,10))
    for i in range(2):
        ax[1,i].set_xlabel('x (km)')
        ax[i,0].set_ylabel('y (km)')
    ax[0,0].set_title('$\eta$')
    ax[0,1].set_title('velocity')
    ax[1,0].set_title('u')
    ax[1,1].set_title('v')
    ax[0,0].contour(x/1000, y/1000, eta.transpose())
    ax[1,0].contour((x + shift[grid][0] * dx)/1000,
                  (y + shift[grid][1] * dx)/1000, u.transpose())
    ax[1,1].contour((x + shift[grid][1] * dx)/1000,
                  (y + shift[grid][0] * dx)/1000, v.transpose())

    if grid == 3:
        ax[0,1].quiver(x[1:]/1000., y[1:]/1000,
                       0.5 * (u[1:,1:] + u[:-1,1:]).transpose(),
                       0.5 * (v[1:,1:] + v[1:,:-1]).transpose())
    else:
        ax[0,1].quiver(x/1000, y/1000, u.transpose(), v.transpose())

        