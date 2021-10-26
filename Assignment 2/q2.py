

import sys, os
import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import scipy as sp
import scipy.sparse as sps
import scipy.sparse.linalg as LAS
import math
from matplotlib.gridspec import GridSpec

from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

# Setup basic domain
L = 1
D = 1

N = 10    # Number of lattice points
h = 0.1  # Spatial step size

# Create the x-grid
x = np.linspace(0, L, N+1)

# Setup final solution time + time-step
tf = 1
k  = 0.001
n_steps = round(tf / k)

# Print the stability condition; Recall that this ratio must be less than 1/2
print('h = ', h)
print('D k / h2 = ', (D * k / h**2))

# Setup to create a movie
n_frames = 10
frame_step = round(n_steps / n_frames)

# Assemble the matrix
d1 = (1.0 - (2. * D * k / h**2)) * np.ones(N+1)
d2 = (D * k / h**2) * np.ones(N)

# Assemble the matrices -> See documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.diags.html
A = sps.diags((d2, d1, d2), offsets=(-1, 0, 1))

# Set initial condition
w = np.ones(11)
peak=np.ones(6)
u_t0 = w
time=['0','0.1','0.03','0.06','1.0','10']
print("t = 0: ", u_t0)
plt.scatter(x,w)
plt.xlabel('length')
plt.ylabel('concentration')
plt.title('Concentration at different time')

for j in range(n_steps):
    w[:] = A.dot(w)+k*math.sin((j+1)/1000)
    if (j == int(0.05*n_steps-1)):
        u_t005 = w
        u_t005[0]=1
        u_t005[10]=1
        peak[0]=u_t005[6]
        print("t = 0.005: ", u_t005)
        plt.scatter(x, u_t005)
                
    if (j == int(0.1*n_steps-1)):
        u_t01 = w
        u_t01[0]=1
        u_t01[10]=1
        peak[1]=u_t01[6]
        print("t = 0.01: ", u_t01)
        plt.scatter(x, u_t01)
        
    if (j == int(0.2*n_steps-1)):
        u_t02 = w
        u_t02[0]=1
        u_t02[10]=1
        peak[2]=u_t02[6]
        print("t = 0.02: ", u_t02)
        plt.scatter(x, u_t02)
        
    if (j == int(1*n_steps-1)):
        u_t1 = w
        u_t1[0]=1
        u_t1[10]=1
        peak[3]=u_t1[6]
        print("t = 1: ", u_t1)
        plt.scatter(x, u_t1)
        
    if (j == int(10*n_steps-1)):
        u_t10 = w
        u_t10[0]=1
        u_t10[10]=1
        peak[10]=u_t10[6]
        print("t = 1: ", u_t10)
        plt.scatter(x, u_t10)

        
    # enforce BC
    w[0] = 1
    w[N] = 1

fig = plt.figure(figsize=(10, 5))
gs = GridSpec(nrows=2, ncols=2)

ax0 = fig.add_subplot(gs[0, 0])
ax0.plot(time, peak)
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Oscillation of midpoint concentration')
