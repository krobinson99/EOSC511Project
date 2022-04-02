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
u0 = 5   #advection velocity km/s 
dt = 2   #timestep seconds 
C0 = 100    #Initial concentration
tmax = 1000  #number of timesteps in order to compute a full cycle
#Initial cond
# itions and vectors
x = np.arange(0,L,dx) 
z = np.arange(0,D,dz) 

u = u0*np.linspace(-1,1,int(D/dz))
C = np.zeros([tmax,N,M])
C[0,23:27,0:D] = C0
X,U=np.meshgrid(x,u)
X,Z=np.meshgrid(x,z)
U=U.T
# plt.contourf(X,Z,U.T)
# plt.show()
#Define the scheme functions 
#EULER FORWARD UPSTREAM SCHEME
def upstr(c,ca,uj):
    cf = c - dt*(abs(uj)*(c - ca)/dx)
    return cf

#LAX-WENDROFF SCHEME
def lax(ca,c,cs,uj):
    cf = c - dt*(uj*(cs - ca)/(2*dx)) + dt*((uj**2)*(dt/2)*(cs -2*c + ca)/(dx**2))
    return cf


######################################

###############  UNCOMMENT TO USE EULER SCHEME ###############

# for n in range(0,tmax-1):
#     for i in range(0,N-1):
#         for j in range(0,M):
#             if U[i,j]>=0:
#                 C[n+1,i,j] = upstr(C[n,i,j],C[n,i-1,j],U[i,j])
#             else:
#                 C[n+1,i,j] = upstr(C[n,i,j],C[n,i+1,j],U[i,j])
# titulo=str('Concentration advection computed with Euler Math.Scheme')

##############  UNCOMMENT TO USE LAX SCHEME  ###############

# for n in range(0,tmax-1):
#     for i in range(0,N-1):
#         for j in range(0,M):
#             C[n+1,i,j] = lax(C[n,i-1,j],C[n,i,j],C[n,i+1,j],U[i,j])
        
# titulo=str('Concentration advection computed with Lax Wendrof Math.Scheme')


#Generate animation and play it in loop but don't save it
fig = plt.figure()
def animate(frames):
    if frames > tmax:
        factor=int(frames/tmax)
        frames=frames-(tmax*factor)
        
    plt.clf()
    ax = fig.add_subplot(111)
    a = C[frames,:,:] 
    oo = plt.contourf(X,-Z,a.T,vmin=0,vmax=100)
    plt.title(titulo)
    plt.ylabel('Concentration')
    plt.xlabel('x space (km)')
    plt.grid()
    
    return oo

ani = animation.FuncAnimation(fig, animate, frames=np.arange(0,tmax,10))
plt.show()

