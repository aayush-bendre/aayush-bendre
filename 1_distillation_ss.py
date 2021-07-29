import numpy as np
from scipy.optimize import fsolve

stages = 5
n = stages+2
stagef = 3
al = 1.6
rr = 3
F = 1
D = 0.5*F
z = 0.5
Lr = rr*D
V = Lr + D
Ls = Lr + F
B = F-D
L1 = V + B
L2 = Lr+F
L3 = L2
L4 = Lr
L5 = Lr
k = F*z


def myFunction(z):
   xb = z[0]
   x1 = z[1]
   x2 = z[2]
   x3 = z[3]
   x4 = z[4]
   x5 = z[5]
   xd = z[6]
   yb = z[7]
   y1 = z[8]
   y2 = z[9]
   y3 = z[10]
   y4 = z[11]
   y5 = z[12]
   
   
   F = np.empty(13)
   F[0] = yb - (al*xb/(1+((al-1)*xb)))
   F[1] = y1 - (al*x1/(1+((al-1)*x1)))
   F[2] = y2 - (al*x2/(1+((al-1)*x2)))
   F[3] = y3 - (al*x3/(1+((al-1)*x3)))
   F[4] = y4 - (al*x4/(1+((al-1)*x4)))
   F[5] = y5 - (al*x5/(1+((al-1)*x5)))
   F[6] = (V*y5)-((Lr+D)*xd)   #xd
   F[7] = (Lr*xd)-(L5*x5)+(V*y4)-(V*y5)    #x5
   F[8] = (L5*x5)-(L4*x4)+(V*y3)-(V*y4)    #x4
   F[9] = (L4*x4)-(L3*x3)+(V*y2)-(V*y3)+k    #x3
   F[10] = (L3*x3)-(L2*x2)+(V*y1)-(V*y2)    #x2
   F[11] = (L2*x2)-(L1*x1)+(V*yb)-(V*y1)    #x1
   F[12] = (L1*x1)-(V*yb)-(B*xb)    #xb
   return F

guess = np.array([0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5])
z = fsolve(myFunction,guess)

x = np.zeros(n)

for i in range(n):
    x[i] = z[i]
print(x)
    
    