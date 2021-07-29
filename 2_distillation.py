import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# define model
def distill(x,t,rr,Feed,x_Feed):
    
    # Parameters
    # Distillate Flowrate (mol/min)
    D=0.5*Feed
    # Flowrate of the Liquid in the Rectification Section (mol/min)
    L=rr*D
    # Vapor Flowrate in the Column (mol/min)
    V=L+D
    # Flowrate of the Liquid in the Stripping Section (mol/min)
    FL=Feed+L
    # Relative Volatility = (yA/xA)/(yB/xB) = KA/KB = alpha(A,B)
    vol=1.6
    # Total Molar Holdup in the Condenser
    atray=0.25
    # Total Molar Holdup on each Tray
    acond=0.5
    # Total Molar Holdup in the Reboiler
    areb=1.0

    #Volatility
    y = np.empty(len(x))
    for i in range(32):
        y[i] = x[i] * vol/(1.0+(vol-1.0)*x[i])

    # Compute xdot
    xdot = np.empty(len(x))
    xdot[0] = 1/acond*V*(y[1]-x[0])
    for i in range(1,16):
        xdot[i] = 1.0/atray*(L*(x[i-1]-x[i])-V*(y[i]-y[i+1]))
    xdot[16] = 1/atray*(Feed*x_Feed+L*x[15]-FL*x[16]-V*(y[16]-y[17]))
    for i in range(17,31):
        xdot[i] = 1.0/atray*(FL*(x[i-1]-x[i])-V*(y[i]-y[i+1]))
    xdot[31] = 1/areb*(FL*x[30]-(Feed-D)*x[31]-V*y[31])
    return xdot

# Steady State Initial Conditions for the 32 tstages
x_ss =np.array([0.935,0.900,0.862,0.821,0.779,0.738,\
0.698,0.661,0.628,0.599,0.574,0.553,0.535,0.521,    \
0.510,0.501,0.494,0.485,0.474,0.459,0.441,0.419,    \
0.392,0.360,0.324,0.284,0.243,0.201,0.161,0.125,    \
0.092,0.064])
x0 = x_ss

# Steady State Initial Condition
rr_ss = 3.0

# Time Interval (min)
t = np.linspace(0,10,100)

# Store results for plotting
xd = np.ones(len(t)) * x_ss[0]
rr = np.ones(len(t)) * rr_ss
ff = np.ones(len(t))
xf = np.ones(len(t)) * 0.5
sp = np.ones(len(t)) * 0.97

# Step in reflux ratio
#rr[10:] = 4.0
#rr[40:] = 2.0
#rr[70:] = 3.0

# Feed Concentration (mol frac)
#xf[50:] = 0.42

# Feed flow rate
ff[80:] = 0.2

# Simulate
for i in range(len(t)-1):
    ts = [t[i],t[i+1]]
    y = odeint(distill,x0,ts,args=(rr[i],ff[i],xf[i]))
    xd[i+1] = y[-1][0]
    x0 = y[-1]

plt.figure()
plt.subplot(3,1,1)
plt.plot(t,rr,'b--',linewidth=3)
plt.ylabel('$RR$')
plt.legend(['Reflux Ratio'],loc='best')


plt.subplot(3,1,2)
plt.plot(t,xf,'k:',linewidth=3,label='Feed composition')
#plt.plot(t,ff,'g-',linewidth=3,label='Feed Flow (mol/min)')
plt.ylabel('Feed')
plt.ylim([0,1.1])
plt.legend(loc='best')


plt.subplot(3,1,3)
plt.plot(t,xd,'r-',linewidth=3,label='Distillate composition')
plt.plot(t,sp,'k.',linewidth=1,label='Set Point')
plt.ylabel('$x_d')
plt.xlabel('Time (min)')
plt.ylim([0.9,1])
plt.legend(loc='best')

plt.savefig('distillation.png')
plt.show()

# Construct results and save data file
# Column 1 = time
# Column 2 = reflux ratio
# Column 3 = distillate composition
data = np.vstack((t,rr,xd)) # vertical stack
data = data.T             # transpose data
np.savetxt('data.txt',data,delimiter=',')


