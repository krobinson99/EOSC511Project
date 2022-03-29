#!/usr/bin/env python
import numpy as np

def glacier(ngridx, ngridz, dt, zinput, T, motion = False, steady = True):  # return eta
    '''recommended values ngridx=50, ngridz = 20, dt=200, T=10*86400 (10 days)
    if motion = True motion case for BCs, initial, stepper (eventually) will be used
    '''
    g = 10 
    D = 200            # depth of our domain in x direction [m]
    L = 20e3           # length of our domain in x direction [m]
    C0 = 10           # input concentration of methane          NOT TRUE
    S0 = 0             # Input concentration of salinity 
    dx = L/(ngridx-1)
    dz = D/(ngridz-1)
    u0 = 0             # Input velocity of FW plume              NEEDS A REAL VALUE
    zz = int(zinput/dz)        # *** Scale so input height matches grid *** 
    Kx= 5 * L/dx
    Kz= 1e-4 * D/dz
    alpha =  0.4/86400    # Oxidation rate constant (0.4 day^-1, converted to seconds).
    mu = 0.00183       # Dynamic viscosity of water at 1 C [m2/s]
    Kd = 0.0087e-5      # molecular diffusivity of methane at 4C in seawater (couldn't find for 1C) [m2/s]
    #Sc = mu/(Kd*rho)    # Schmidt number for water at 1 C.
    zn = np.linspace(0,D,ngridz)
    Sop = 33 + np.log(1e-3+zn/D)
    Cop = 3.7

# set up temporal scale T is total run time
    ntime = int(T/dt)
    if motion==False:
    # initialize
        C, S = init0(ngridx,ngridz,ntime,motion)
        C[0,:,:], S[0,:,:] = initial_steady(C[0,:,:],S[0,:,:],Sop)
        C[0,:,:], S[0,:,:] = boundary_steady(C[0,:,:], S[0,:,:], C0, S0,zz,Sop)
    # main loop (Euler forward)
        for nt in range(1,ntime):
            if steady:
                C[nt,:,:], S[nt,:,:] = stepper_steady(dx,dz,dt,C[nt-1,:,:],S[nt-1,:,:],Kx,Kz)
            else:
                C[nt,:,:], S[nt,:,:] = stepper_sink(dx,dz,dt,C[nt-1,:,:],S[nt-1,:,:],Kx,Kz,alpha,Kd,ngridz)
    # periodic boundary conditions
            C[nt,:,:], S[nt,:,:] = boundary_steady(C[nt,:,:], S[nt,:,:], C0, S0,zz,Sop)
        C = C + Cop
        return C,S


def init0(ngridx,ngridz,ntime,motion):
    '''initialize a ngrid x ngrid domain, u, v,, all zero 
     we need density salinity, ch4''' 
    S = np.ones((ntime,ngridx, ngridz))
    C = np.ones_like(S)
    if motion:
        u = np.zeros_like(S)
        w = np.zeros_like(S)
        rho = np.zeros_like(S)
        return  u, w, rho, S, C 
    else:
        return  C, S 

def diffx(C,dx):
    Cdx=np.zeros_like(C)
    for i in range(C.shape[0]):
        for j in range(C.shape[1]):
            if i==0:
                Cdx[i,j]= (dx**2) * (2*C[i,j] - 5*C[i+1,j] + 4*C[i+2,j] - C[i+3,j])/(dx**3)
            elif i==C.shape[0]-1:
                Cdx[i,j]= (dx**2) * (2*C[i,j] - 5*C[i-1,j] + 4*C[i-2,j] - C[i-3,j])/(dx**3)
            else:
                Cdx[i,j]=C[i+1,j]-2*C[i,j]+C[i-1,j]
    return Cdx

def diffz(C,dz):
    Cdz=np.zeros_like(C)
    for i in range(C.shape[0]):
        for j in range(C.shape[1]):
            if j==0:
                Cdz[i,j]= (dz**2) * (2*C[i,j] - 5*C[i,j+1] + 4*C[i,j+2] - C[i,j+3])/(dz**3)
            elif j==C.shape[1]-1:
                Cdz[i,j]= (dz**2) * (2*C[i,j] - 5*C[i,j-1] + 4*C[i,j-2] - C[i,j-3])/(dz**3)
            else:
                Cdz[i,j]=C[i,j+1]-2*C[i,j]+C[i,j-1]    
    return Cdz

def boundary_steady(C, S, C0, S0, zz ,Sop):
    '''Sets the boundary conditions for the steady state if motion = False, boundaries for motion case if true'''
    ## open water boundary
    C[-1, :] = C[-2,:] ## nM
    S[-1, :] = Sop ## PSU a function for S with depth 

    ## Glacier wall
    C[0, :] = C[1,:]
    C[0, zz] = C0
    S[0, :] = S[1,:]
    S[0, zz] = S0
    return C,S

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

def initial_steady(C, S,Sop):
    ''' sets the inital conditions for the steady state stages'''
    C  = 0*C ## 3.7*C # Changed from 4.5 to 3.7 because solubility of methane at 33 PSU and 0.5 C (closest to our conditions) is 3.7 nM  
    for i in range(S.shape[0]):
        S[i,:] = Sop # Changed to 33 from 35, closest to actual environmental conditions
    return C, S

def inital_motion(C, S, u, w):
    ''' sets initial conditions for motion case '''
    C, S = initial_steady(C, S)       ## A placeholder until we have a steady field solution. maybe use the steady state stepper
    u, w = 0                          ## should be set when creating grids anyway
    return C, S
    
def stepper_steady(dx,dz,dt,C,S,Kx,Kz):
    Cn = C + dt*(diffx(C,dx)*Kx/(dx**2)+Kz*diffz(C,dz)/(dz**2))
    Sn = S + dt*(diffx(S,dx)*Kx/(dx**2)+Kz*diffz(S,dz)/(dz**2))
    return Cn,Sn

def stepper_sink(dx,dz,dt,C,S,Kx,Kz,alpha,Kd,ngridz):
    ML = 10
    Dp = int(ML/(ngridz-1))   ## space step closest to mixing layer for MLayer = 20m
    Cn = C + dt*(diffx(C,dx)*Kx/(dx**2)+Kz*diffz(C,dz)/(dz**2)-alpha*C)
    Cn[:,0:Dp] = Cn[:, 0:Dp] -dt*(Kd/(ML*0.000025)*(C[:,0:Dp]))
    #Cn[:,0] = C[:,0] + dt*(diffx(C[:,0])*Kx/(dx**2)+Kz*diffz(C[:,0])/(dz**2)-Ro*C[:,0]-Kd/(dz*0.000025)*(C[:,0] - 3.7))
    Sn = S + dt*(diffx(S,dx)*Kx/(dx**2)+Kz*diffz(S,dz)/(dz**2))
    return Cn,Sn

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

