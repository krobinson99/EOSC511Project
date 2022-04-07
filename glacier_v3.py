#!/usr/bin/env python
from curses import A_HORIZONTAL
import numpy as np
from scipy.sparse import spdiags

def glacier(ngridx, ngridz, dt, zinput, T, ML, alpha, motion = False, steady = True):  # return eta
    '''recommended values ngridx=50, ngridz = 20, dt=60, zinput=0 (unstable if motion=True at other depths), T=5*86400 (5 days),
    ML at least dz, alpha=0.4/86400.
    examples: 
    C,S,u,w = glacier(ngridx, ngridz, dt, zinput, T, ML, alpha, motion = True,steady=false)
    C,S= glacier(ngridx, ngridz, dt, zinput, T, ML, alpha,steady=False)
    '''
    g = 10              # gravity
    rhoc = 1026         # standard density sw
    gr = g/rhoc        
    D = 200             # depth of our domain [m]
    L = 20e3            # length of our domain in x direction [m]
    C0 = 10             # input concentration of methane          NOT TRUE
    S0 = 0              # Input concentration of salinity 
    dx = L/(ngridx)     # spatial resolution x direction
    dz = D/(ngridz)     # spatial resolution z direction
    Q = 5e-4            # inflow at glacier wall m/s
    u0 = Q * (D/(2*dz)) # velocity inflow                          NOT WORKING
    zz = int(zinput/dz) # *** Scale so input height matches grid *** 
    Kx= 5 * L/dx        # diffusivity coeff x direction
    Kz= 1e-4 * D/dz     # diffusivity coeff z direction
    Ah= 2 * L/dx        # eddy horizontal diff coeff 
    Av= 1e-3 * D/dz     # eddy vertical diff coeff 
    #alpha =  0.4/86400 # Oxidation rate constant (0.4 day^-1, converted to seconds).
    mu = 0.00183       # Dynamic viscosity of water at 1 C [m2/s]
    Kd = 0.0087e-5      # molecular diffusivity of methane at 4C in seawater (couldn't find for 1C) [m2/s]
    #Sc = mu/(Kd*rho)    # Schmidt number for water at 1 C.
    zn = np.linspace(0,D,ngridz+1)
    Sop = 33 #+ np.log(1e-3+zn/D) #initial condition of salinity
    Cop = 3.7                     #[ch4] at equilibrium with atmosphere
# set up temporal scale T is total run time
    ntime = int(T/dt)
    if motion==False:
    # initialize
        C, S = init0(ngridx,ngridz,ntime,motion)
        S[0,:,:] = initial(S[0,:,:],Sop,dz)
        C[0,:,:], S[0,:,:] = boundary_steady(C[0,:,:], S[0,:,:], C0, S0,zz,Sop)
    # main loop (Euler forward)
        for nt in range(1,ntime):
            if steady:
                C[nt,:,:], S[nt,:,:] = stepper_sink(dx,dz,dt,C[nt-1,:,:],S[nt-1,:,:],Kx,Kz,alpha,Kd,ngridz,ML)
            else:
                C[nt,:,:], S[nt,:,:] = stepper_sink(dx,dz,dt,C[nt-1,:,:],S[nt-1,:,:],Kx,Kz,alpha,Kd,ngridz,ML,steady)
    # periodic boundary conditions
            C[nt,:,:], S[nt,:,:] = boundary_steady(C[nt,:,:], S[nt,:,:], C0, S0,zz,Sop)
        C = C + Cop
        return C,S
    if motion:
    # initialize
        u, w, S, C = init0(ngridx,ngridz,ntime,motion)
        S[0,:,:] = initial(S[0,:,:],Sop,dz)
        C[0,:,:], S[0,:,:], u[0,:,:], w[0,:,:] = boundary_motion(C[0,:,:], S[0,:,:],u[0,:,:], w[0,:,:], u0, C0, S0, zz, D, Sop)
    # main loop (Euler forward)
        for nt in range(1,ntime):
            C[nt,:,:], S[nt,:,:],u[nt,:,:],w[nt,:,:]= stepper_motion(gr,dx,dz,dt,C[nt-1,:,:],S[nt-1,:,:],Kx,Kz,u[nt-1,:,:],w[nt-1,:,:],D,Q,Ah,Av,alpha,ML,Kd,steady)
    # periodic boundary conditions
            C[nt,:,:], S[nt,:,:], u[nt,:,:], w[nt,:,:] = boundary_motion(C[nt,:,:], S[nt,:,:],u[nt,:,:], w[nt,:,:], u0, C0, S0, zz, D, Sop)
        C = C + Cop
        return C,S,u,w


def init0(ngridx,ngridz,ntime,motion):
    '''initialize a ngrid x ngrid domain, u, v, all zero 
     we need density salinity, ch4''' 
    S = np.ones((ntime,ngridx+1, ngridz+1))  #Arakawa B grid with S,C in the gridpoints.
    C = np.zeros_like(S)
    if motion:
        u = np.zeros((ntime,ngridx, ngridz))
        w = np.zeros_like(u)
        return  u, w, S, C
    else:
        return  C, S 
    
def initial(S,Sop,dz,motion=False):
    ''' sets the inital conditions for the steady state stages'''
    #C  = 0*C ## 3.7*C # Changed from 4.5 to 3.7 because solubility of methane at 33 PSU and 0.5 C (closest to our conditions) is 3.7 nM  
    for i in range(S.shape[0]):
        S[i,:] = Sop # Changed to 33 from 35, closest to actual environmental conditions
    return S

def dens(S,dz):
    a0 = 1.665e-1
    b0 = 7.6554e-1
    lam1 = 5.952e-2
    lam2 = 5.4914e-4
    cab = 2.4341e-3
    nu1 = 1.497e-4 
    nu2 = 1.109e-5
    rhoc=1026
    Ta = -9
    Sa = S - 35
    z = np.arange(0,S.shape[1]*dz,dz)
    Z=np.tile(z, (S.shape[0],1)) 
    return rhoc*(1+(-a0*(1 + 0.5*lam1*Ta + nu1*Z)*Ta + b0*(1 - 0.5*lam2*Sa - nu2*Z)*Sa - cab*Ta*Sa)/rhoc)

def boundary_steady(C, S, C0, S0, zz ,Sop):
    '''Sets the boundary conditions for the steady state if motion = False, boundaries for motion case if true'''
    ## open water boundary
    C[-1, :] = C[-2,:] ## nM
    S[-1, :] = S[-2, :] ## PSU a function for S with depth 
    #surface
    C[:, 0] = C[:,1] 
    S[:, 0] = S[:,1] 
    #Bottom
    C[:, -1] = C[:,-2] 
    S[:, -1] = S[:,-2] 
    ## Glacier wall
    C[0, :] = C[1,:]
    C[0, zz] = C0
    C[1, zz] = (C0 + C[2, zz])/2
    S[0, :] = S[1,:]
    S[0, zz] = S0
    S[1, zz] = (S0 + S[2, zz])/2
    return C,S

def boundary_motion(C, S, u, w, u0, C0, S0, zz, D,Sop):
    C, S = boundary_steady(C, S, C0, S0, zz,Sop)
    ## Surface and bottom boundary
    w[:, 0] = w[:, -1] = 0
    u[:, -1] = 0 #u[:, -2] #With friction we get no flow boundary
    u[:, 0] = u[:, 1]  #overwrites velocities surface!!!!
    ## Wall boundary
    u[0, :] = 0
    w[0, :] = 0#w[1, :]
    #if zz == 0:
    #     u[0, zz] = 2*u0 
    #else:
    #    u[0, zz-1:zz+1] = u0     #inflow over two points (allows us to calculate flow at gridpoint C,S)
    ## open boundary
    w[-1, :] = w[-2, :]
    u[-1, :] = u[-2, :]  
    return C, S, u, w
    
def stepper_sink(dx,dz,dt,C,S,Kx,Kz,alpha,Kd,ngridz,ML,steady=True):
    if steady:
        Cn = C + dt*(diffx(C,dx)*Kx+Kz*diffz(C,dz))
    else:
        Dp = int(ML/(ngridz-1))   ## space step closest to mixing layer for MLayer = 20m
        Cn = C + dt*(diffx(C,dx)*Kx+Kz*diffz(C,dz)-alpha*C) #oxidation
        Cn[:,0:Dp] = Cn[:, 0:Dp] -dt*(Kd/(ML*0.000025)*(C[:,0:Dp])) #evasion
    Sn = S + dt*(diffx(S,dx)*Kx+Kz*diffz(S,dz))
    return Cn,Sn

def stepper_motion(gr,dx,dz,dt,C,S,Kx,Kz,u,w,D,Q,Ah,Av,alpha,ML,Kd,steady):
    Dp = int(ML/(C.shape[1]))
    Uc,Wc,rhox = grid_mid(u,w,S,dz,dx)
    if steady:
        Cn = C + dt*(diffx(C,dx)*Kx + Kz*diffz(C,dz) - upstr(C,Uc,0)/dx  - upstr(C,Wc,1)/dz)
    else:
        Cn = C + dt*(diffx(C,dx)*Kx + Kz*diffz(C,dz) - upstr(C,Uc,0)/dx  - upstr(C,Wc,1)/dz -alpha*C)
        Cn[:,0:Dp] = Cn[:, 0:Dp] -dt*(Kd/(ML*0.000025)*(C[:,0:Dp]))
    Sn = S + dt*(diffx(S,dx)*Kx + Kz*diffz(S,dz) - upstr(S,Uc,0)/dx - upstr(S,Wc,1)/dz)
    z = np.arange(0,u.shape[1]*dz,dz)
    Z=np.tile(z, (u.shape[0],1)) 
    un =  u + dt*(diffx(u,dx)*Ah + diffz(u,dz)*Av -gr*p_grad(rhox,dz) - upstr(u,u,0)/dx  - upstr(u,w,1)/dz) 
    q = Q - (np.sum(un*dz,axis=1))/D
    un = un + np.tile(q, (u.shape[1],1)).T 
    wn = vert_vel(un,dz,dx)
    return Cn,Sn,un,wn

def vert_vel(u,dz,dx):
    w = np.zeros_like(u)
    for i in range(u.shape[0]-3,0,-1):
        for j in range(u.shape[1]-2,0,-1):
            w[i,j] = w[i,j+1] + (dz/dx)*(u[i+1,j]-u[i,j])
    return w

def p_grad(rhox,dz):
    pgrad = np.zeros_like(rhox)
    for j in range(rhox.shape[1]):
        pgrad[:,j] = 0.5*dz*rhox[:,j] + np.sum(dz*rhox[:,0:j],axis=1)
    return pgrad

def grid_mid(u,w,S,dz,dx):
    rho = dens(S,dz)
    rhox = np.zeros_like(u)
    uc = np.zeros_like(rho)
    wc = np.zeros_like(rho)
    for j in range(rho.shape[1]):
        if j==0:
            uc[:,j] = matprod(u,0.5,j)
            wc[:,j] = matprod(w,0.5,j)
        if j==u.shape[1]:
            uc[:,j] = matprod(u,0.5,j-1)
        else:
            uc[:,j] = matprod(u,0.25,j) + matprod(u,0.25,j-1)
            wc[:,j] = matprod(w,0.25,j) + matprod(w,0.25,j-1)
        for i in range(1,rho.shape[0]):
            if i < rho.shape[0]-1 and j < rho.shape[1]-1:
                rhox[i,j] = ((rho[i+1,j]+rho[i+1,j+1])/2 - (rho[i,j]+rho[i,j+1])/2)/dx
                
    return uc,wc,rhox

def matprod(U,val,ind):
    n = U.shape[0]
    e = np.ones([1,n])
    dat = np.vstack((val*e,val*e)) 
    diags = np.array([-1,0])
    A = spdiags(dat, diags, n, n + 1).toarray()
    A[0,0:2] = np.array([2*val, 0])
    A[-1,-2:] = np.array([0, 2*val])
    return U[:,ind]@A
        
def diffx(C, dx): 
    n   = C.shape[0];
    e   = np.ones([1,n])
    dat = np.vstack((e,-2*e,e))/(dx**2)# Centered difference approximation, second derivative.
    diags = np.array([-1,0,1])
    A = spdiags(dat, diags, n, n).toarray()
    
    # Sets simple forward and backward difference scheme at end points, change as needed w/ BCs
    A[0,0:4] = np.array([2, -5, 4, -1])*(dx**2)/(dx**3) # Forward difference approximation, second derivative
    A[-1,-4:] = np.array([-1, 4, -5, 2])*(dx**2)/(dx**3) # Backward difference approximation, second derivative.

    Cdx = A@C
    return Cdx

def diffz(C,dz): 
    n   = C.shape[1];
    e   = np.ones([1,n])
    dat = np.vstack((e,-2*e,e))/(dz**2) # Centered difference, second derivative
    diags = np.array([-1,0,1])
    A   = spdiags(dat, diags, n, n).toarray()
    
    # Sets simple forward and backward difference scheme at end points, change as needed w/ BCs
    A[0,0:4] = np.array([2, -5, 4, -1])/(dz**3) # Forward difference approximation, second derivative
    A[-1,-4:] = np.array([-1, 4, -5, 2])/(dz**3) # Backward difference approximation, second derivative.

    Cdz = A@C.transpose()
    Cdz = Cdz.transpose()
    return Cdz

def upstr(C,U,ax):
    adv=np.zeros_like(C)
    for i in range(0,C.shape[0]-1):
        for j in range(0,C.shape[1]-1):
            if ax == 0:
                if U[i,j]>=0:
                    adv[i,j] = U[i,j]*(C[i,j] - C[i-1,j])
                else:
                    adv[i,j] = abs(U[i,j])*(C[i,j] - C[i+1,j])
            else:
                if U[i,j]>=0:
                    adv[i,j] = U[i,j]*(C[i,j] - C[i,j-1])
                else:
                    adv[i,j] = abs(U[i,j])*(C[i,j] - C[i,j+1])
    return adv
