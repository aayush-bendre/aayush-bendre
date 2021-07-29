import numpy as np
import sympy as sp
from gekko import GEKKO

m = GEKKO()
m.options.SOLVER=1 #APOPT

x , y, k11, k12, k13, k21, k22, k23 = sp.symbols('x y k11 k12 k13 k21 k22 k23')

f = 2*x**2 + 3*y**2
g = (1/x)+(1/y)-4

dfdx = sp.Derivative(f,x).doit()
dfdy = sp.Derivative(f,y).doit()
dgdx = sp.Derivative(g,x).doit()
dgdy = sp.Derivative(g,y).doit()

def func (X,Y,ix,iy):

    X0 = [X[ix],Y[iy]]
    #print(X0)

    x0 = X0[0]
    y0 = X0[1]
    #print("x=",x0)
    #print("y=",y0)
    
    fx0 = f.subs(x,x0).subs(y,y0)
    #print("Z=",fx0)

    gx0 = g.subs(x,x0).subs(y,y0)
    #print("G=",gx0)

    df = np.array([dfdx.subs(x,x0).subs(y,y0),dfdy.subs(y,y0).subs(x,x0)])
    #print(df)
    dg = np.array([dgdx.subs(x,x0).subs(y,y0),dgdy.subs(y,y0).subs(x,x0)])
    #print(dg)

    dXn = np.array([X[ix-1]-X[ix],X[ix]-X[ix],X[ix+1]-X[ix]])
    #print("dXn:",dXn)
    Kx = np.array([[k11],[k12],[k13]])

    dYn = np.array([Y[iy-1]-Y[iy],Y[iy]-Y[iy],Y[iy+1]-Y[iy]])
    Ky = np.array([[k21],[k22],[k23]])
    #print("dYn:",dYn)
    
    dX1 = np.dot(dXn,Kx)
    dY1 = np.dot(dYn,Ky)

    dX = np.array([[dX1[0]],[dY1[0]]])
    #print("dx",dX)

    ff = fx0 + np.dot(df,dX)
    gg = gx0 + np.dot(dg,dX)

    F = ff[0]
    G = gg[0]
    #print(F,G)
    
    Af = np.zeros(7)
    Af[0] = F.subs(k11,0).subs(k12,0).subs(k13,0).subs(k21,0).subs(k22,0).subs(k23,0) #constant
    Af[1] = F.subs(k11,1).subs(k12,0).subs(k13,0).subs(k21,0).subs(k22,0).subs(k23,0) - Af[0] #k11
    Af[3] = F.subs(k11,0).subs(k12,0).subs(k13,1).subs(k21,0).subs(k22,0).subs(k23,0) - Af[0] #k13
    Af[4] = F.subs(k11,0).subs(k12,0).subs(k13,0).subs(k21,1).subs(k22,0).subs(k23,0) - Af[0] #k21
    Af[6] = F.subs(k11,0).subs(k12,0).subs(k13,0).subs(k21,0).subs(k22,0).subs(k23,1) - Af[0] #k23
    #print(Af)

    Ag = np.zeros(7)
    Ag[0] = G.subs(k11,0).subs(k12,0).subs(k13,0).subs(k21,0).subs(k22,0).subs(k23,0) #constant
    Ag[1] = G.subs(k11,1).subs(k12,0).subs(k13,0).subs(k21,0).subs(k22,0).subs(k23,0) - Ag[0] #k11
    Ag[3] = G.subs(k11,0).subs(k12,0).subs(k13,1).subs(k21,0).subs(k22,0).subs(k23,0) - Ag[0] #k13
    Ag[4] = G.subs(k11,0).subs(k12,0).subs(k13,0).subs(k21,1).subs(k22,0).subs(k23,0) - Ag[0] #k21
    Ag[6] = G.subs(k11,0).subs(k12,0).subs(k13,0).subs(k21,0).subs(k22,0).subs(k23,1) - Ag[0] #k23
    #print(Ag)

    x11 = m.Var(value=0,lb=0,ub=1,integer=True)
    x12 = m.Var(value=0,lb=0,ub=1,integer=True)
    x13 = m.Var(value=0,lb=0,ub=1,integer=True)
    x21 = m.Var(value=0,lb=0,ub=1,integer=True)
    x22 = m.Var(value=0,lb=0,ub=1,integer=True)
    x23 = m.Var(value=0,lb=0,ub=1,integer=True)

    X = np.array([1,x11,x12,x13,x21,x22,x23])

    Fx = np.dot(Af,X.T)
    Gx = np.dot(Ag,X.T)

    m.Equation(Gx<=0)
    m.Equation(x11 + x12 + x13 == 1)
    m.Equation(x21 + x22 + x23 == 1)
    m.Obj(Fx)

    m.solve(disp=False)

    Xf = [x11.value[0],x12.value[0],x13.value[0]]
    Yf = [x21.value[0],x22.value[0],x23.value[0]]
    
    xf = 0
    yf = 0
    
    if Xf[0]==1:
        xf = 0
    elif Xf[1]==1:
        xf = 1
    else:
        xf = 2
    
    if Yf[0]==1:
        yf = 0
    elif Yf[1]==1:
        yf = 1
    else:
        yf = 2
    
    return fx0,xf,yf
    
# Discrete Values
X = [0.3,0.7,0.8,1.2,1.5,1.8]
Y = [0.4,0.8,1.1,1.4,1.6]

iM = 3
iN = 2

i = 0
z,bx,by = func(X,Y,iM,iN)   


#print(bx)
#print(by)

zmin = np.inf
xmin = 0
ymin = 0

while z < zmin:
    
    
    zmin = z
    print("Z:",zmin)
    xmin = X[iM]
    ymin = Y[iN]
    print("x:",xmin)
    print("y:",ymin)
    
    if bx == 0:
        iM = iM - 1
    elif bx == 2:
        iM = iM + 1
    else:
        iM = iM
    
    if by == 0:
        iN = iN - 1
    elif by == 2:
        iN = iN + 1
    else:
        iN = iN

    z,bx,by = func(X,Y,iM,iN)


print("Optimized solution:",zmin)
print("X =", xmin)
print("Y =", ymin)
