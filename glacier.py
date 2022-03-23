#!/usr/bin/env python
import numpy as np

<<<<<<< HEAD
def glacier(ngridx, ngridz, dt, zinput, T, motion = False):  # return eta
=======
def glacier(ngridx, ngridz, dt, T, zz, mode = 0):  # return eta
>>>>>>> 7acf0d8f0aa8e7e828f32d9ce4df70d2c56eb4fe
    '''recommended values ngrid=11, dt=150, T=4*3600 (4 hours)???? CHANGE FOR OUR PROJECT
    if motion = True motion case for BCs, initial, stepper (eventually) will be used
    '''

    g = 10 
<<<<<<< HEAD
    D = 200            # depth of our domain in x direction [m]
    L = 20e3           # length of our domain in x direction [m]
    C0 = 10            # input concentration of methane          NOT TRUE
    S0 = 0             # Input concentration of salinity 
    zz = zinput        # *** Scale so input height matches grid *** 
    dx = L/(ngridx-1)
    dz = D/(ngridz-1)
    u0 = 0             # Input velocity of FW plume              NEEDS A REAL VALUE
  

# set up temporal scale T is total run time
    ntime = np.int(T/dt)

# initialize
    u, w, rho, S, C = init0(ngridx,ngridz,ntime)
    C, S = inital_steady(C, S)
    C,S=boundaries(C, S, u, w, uo, C0, S0, zz, D, motion)
    return C,S


    # main loop (leap-frog)
    # for k in range(ntime):
    #     u, v, eta = stepper[grid](ngrid, f, g, Hu, Hv, dt, rdx, u, v, eta, up, vp, etap)

    #     # add forcing
    #     t = k * dt / 8640
    #     eta = eta + 0.1 * (1 - np.exp(-t*t)) * spatial

    #     # periodic boundary conditions
    #     u, v, eta = periodicbc(ngrid, u, v, eta)

    #     # exchange values
    #     u, v, eta, up, vp, etap = exchange(u, v, eta, up, vp, etap)


    #return


def init0(ngridx,ngridz,ntime):
    '''initialize a ngrid x ngrid domain, u, v,, all zero 
     we need density salinity, ch4''' 

    u = np.zeros((ntime,ngridx, ngridz))       ## [time, x, z]
    w = np.zeros_like(u)
    S = np.ones((ntime,ngridx-1, ngridz-1))    ## Needs to be 1 for setting steady init conditions. If problem can set this in init func   
    C = np.ones_like(S)
    rho = np.zeros_like(S)
=======
    D = 200           #depth of our domain in x direction [m]
    L = 20e3          #length of our domain in x direction [m]
    C0 = 10           # input concentration of methane NOT TRUE
    S0 = 0            # Input concentration of salinity 
    dx = L/(ngridx-1)
    dz = D/(ngridz-1)
    zz = np.int(zz/dz)

# set up temporal scale T is total run time
    ntime = np.int(T/dt)
    #S, C = init0(ngridx,ngridz,ntime,mode)
    #C[0,:,:],S[0,:,:]=boundary_steady(C[0,:,:], S[0,:,:],C0,S0,zz) #this line is not needed, just to check output
    if mode==0:
        # initialize
        S, C = init0(ngridx,ngridz,ntime,mode)
        #S, C = initial_steady(S,C)
    # main loop (Euler forward)
        for nt in range(ntime):
            C[nt,:,:], S[nt,:,:] = stepper_steady(ngrid,dx,dz,dt,C,S,Kx,Kz)
    # periodic boundary conditions
            C[nt,:,:], S[nt,:,:] = boundary_steady(C, S, C0, S0,zz)
        return C,S


def init0(ngridx,ngridz,ntime,mode):
    '''initialize a ngrid x ngrid domain, u, v,, all zero 
     we need density salinity, ch4''' 
    S = np.zeros((ntime,ngridx-1, ngridz-1))
    C = np.zeros_like(S)
    if mode==1:
        u = np.zeros((ntime,ngridx, ngridz))
        w = np.zeros_like(u)
        rho = np.zeros_like(S)
        return  u, w, rho, S, C 
    else:
        return  S, C 
>>>>>>> 7acf0d8f0aa8e7e828f32d9ce4df70d2c56eb4fe

def stepper_steady(ngrid,dx,dz,dt,C,S,Kx,Kz):




<<<<<<< HEAD
def boundaries(C, S, u, w, uo, C0, S0, zz, D, motion):
    '''Sets the boundary conditions for the steady state if motion = False, boundaries for motion case if true'''
=======
def boundary_steady(C, S, C0, S0,zz):
    '''Sets the boundary conditions for the steady state and for the source and sink stage'''
>>>>>>> 7acf0d8f0aa8e7e828f32d9ce4df70d2c56eb4fe
    ## open water boundary
    C[-1, :] = 4.5 ## nM
    S[-1, :] = 35 ## PSU 
    ## Glacier wall
    C[0, :] = C[1,:]
    C[0, zz] = C0
    S[0, :] = S[1,:]
    S[0, zz] = S0
<<<<<<< HEAD
    
    if motion:
        ''' set BCs for velocity'''
        ## Surface and depth boundary
        w[:, 0] = w[:, D] = 0
        u[:, -1] = u[:, -2]
        u[:, 0] = u[:, 1]
        ## Wall boundary
        w[0, :] = u[0, :] = 0
        u[0, zz] = u0
        ## open boundary
        w[-1, :] = w[-2, :]
        u[-1, :] = w[-1, :]
        
    return C, S, u, w
'''
def boundary_motion(C, S, u, w, uo, C0, S0, zz, D):
    
    C, S = boundary_steady(C, S, C0, S0, zz)
    ## Surface and depth boundary
    w[:, 0] = w[:, D] = 0
    u[:, -1] = u[:, -2]
    u[:, 0] = u[:, 1]
    ## Wall boundary
    w[0, :] = u[0, :] = 0
    u[0, zz] = u0
    ## open boundary
    w[-1, :] = w[-2, :]
    u[-1, :] = w[-1, :]
    
    return C, S, u, w
'''

def initial_steady(C, S):
    ''' sets the inital conditions for the steady state stages'''
    C  = 4.5*C
    S = 35*S
    return C, S

def inital_motion(C, S, u, w):
    ''' sets initial conditions for motion case '''
    C, S = initial_steady(C, S)       ## A placeholder until we have a steady field solution. maybe use the steady state stepper
    u, w = 0                          ## should be set when creating grids anyway
=======
    return C, S



# def stepgridB(ngrid, f, g, Hu, Hv, dt, rdx, u, v, eta, up, vp, etap):
#     '''take a step forward using grid 2 (B)'''
#     n = ngrid
#     nm1 = ngrid-1
#     nm2 = ngrid-2
#     u[1:nm1, 1:nm1] = u[1:nm1, 1:nm1] + 2*dt * (f * vp[1:nm1, 1:nm1]   #u, v center grid points same as C and S
#                                                  -g * (etap[2:n, 1:nm1] - etap[1:nm1, 1:nm1]
#                                                         + etap[2:n, 2:n] - etap[1:nm1, 2:n]) * 0.5 * rdx)
#     v[1:nm1, 1:nm1] = v[1:nm1, 1:nm1] + 2*dt * (-f*up[1:nm1, 1:nm1]
#                                                  - g * (etap[1:nm1, 2:n] - etap[1:nm1, 1:nm1]
#                                                       + etap[2:n, 2:n] - etap[2:n, 1:nm1]) * 0.5 * rdx)
#     eta[1:nm1, 1:nm1] = eta[1:nm1, 1:nm1] + 2*dt * (-(Hu[1:nm1, 1:nm1] * up[1:nm1, 1:nm1] #eta grid points same as u and w
#                                                       - Hu[0:nm2, 1:nm1] * up[0:nm2, 1:nm1]
#                                                          + Hu[1:nm1, 0:nm2] * up[1:nm1, 0:nm2]
#                                                          - Hu[0:nm2, 0:nm2] * up[0:nm2, 0:nm2]
#                                                          + Hv[1:nm1, 1:nm1] * vp[1:nm1, 1:nm1]
#                                                          - Hv[1:nm1, 0:nm2] * vp[1:nm1, 0:nm2]
#                                                          + Hv[0:nm2, 1:nm1] * vp[0:nm2, 1:nm1]
#                                                          - Hv[0:nm2, 0:nm2] * vp[0:nm2, 0:nm2]) * 0.5 * rdx)
#     return u, v, eta
>>>>>>> 7acf0d8f0aa8e7e828f32d9ce4df70d2c56eb4fe
