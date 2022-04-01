##This code helps to compare three different mathematical schemes for numerical modeling of a concentration C advection (last stepspace of the domain conected to first)
###COMMENT(select and ctrl+3) A SHCEME WHEN NOT USING AND UNCOMMENT(ctrl+4) TO USE (from line 63)  
##Changing the timestep(dt) will affect on the performance of the model (line:14)
import numpy as np
import matplotlib.pyplot as plt
import math
import ffmpeg
from numpy import random
from matplotlib import animation
plt.style.use('seaborn-talk')

#Define parameters
N = 50     #number spacesteps
M = 20
J = int(N/2)     #waves number's for fourier
L = 20e3    #2500km
D = 200
dx = L/N    #25km
dz = D/M
Cou = 0.05 #Courant number
u0 = -0.1   #advection velocity km/s 
dt = int(Cou*dx/abs(u0))    #timestep seconds 
C0 = 100    #Initial concentration
tmax = dt*300  #number of timesteps in order to compute a full cycle
#Initial cond
# itions and vectors
x = np.arange(0,L,dx) 
z = np.arange(0,D,dz) 

C = np.zeros([tmax,N,M])
C[0,int(L/(5*dx)):int(L/(10*dx)),int(D/(5*dz)):int(D/(10*dz))] = C0

X,Z=np.meshgrid(x,z)
print(C)
plt.contour(X,Z,C[0,:,:].T)
plt.show()

#Define the scheme functions 
#EULER FORWARD UPSTREAM SCHEME
def upstr(c,ca):
    cf = c - dt*(abs(u0)*(c - ca)/dx)
    return cf

#LAX-WENDROFF SCHEME
def lax(ca,c,cs):
    cf = c - dt*(u0*(cs - ca)/(2*dx)) + dt*((u0**2)*(dt/2)*(cs -2*c + ca)/(dx**2))
    return cf


######################################

###############  UNCOMMENT TO USE EULER SCHEME ###############
# for n in range(0,tmax-1):
#     for j in range(0,N):
#         C[n+1,j] = upstr(C[n,j],C[n,j-1])
# titulo=str('Concentration advection computed with Euler Math.Scheme')

###############  UNCOMMENT TO USE EULER SCHEME negative ###############
for n in range(0,tmax-1):
    for j in range(0,N-1):
        C[n+1,j] = upstr(C[n,j],C[n,j+1])
titulo=str('Concentration advection computed with Euler Math.Scheme')


#Generate animation and play it in loop but don't save it
fig = plt.figure()
def animate(frames):
    if frames > tmax:
        factor=int(frames/tmax)
        frames=frames-(tmax*factor)
        
    plt.clf()
    ax = fig.add_subplot(111)
    ax.set_ylim(-20, C0+20)
    a = C[frames,:]
    oo = plt.plot(X,a)
    plt.title(titulo)
    plt.ylabel('Concentration')
    plt.xlabel('x space (km)')
    plt.grid()
    
    return oo

#ani = animation.FuncAnimation(fig, animate, frames=np.arange(0,tmax,50))
#plt.show()

