import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import math

v0 =1.5
f = 3732 #4500
Y_in = [0.5, 0.1]
C = 62.9e-09
L = 32.9e-03
R = 700
Ga = -0.00115
Gb = 0.005
Bp = 0.62
xfmax = 4000
dxf = (xfmax - f)/500 #500
dxfn = dxf*R*C
xfnmax = xfmax*R*C
a = R*Ga
bl = R*Gb
bg = (C*(R**2))/L
B = v0*bg/Bp
w = 2*np.pi*f*R*C
xfn = f*R*C
NN = 2
tp = 1/xfn
NIT = 500  #1000
n1 = 100
n2 = 250
NTR = 10000
tstep = tp/NIT
xmax = xfnmax
xmin = xfn
ymax = 3
ymin = -3
xstep = (xmax-xmin)/10
ystep = (ymax-ymin)/10

YNEW1 = []
YNEW2 = []
rkk4list = []

def yprim(Y, yprim_k_helper, t):
    if Y[0] >= 1:
        gx = bl*Y[0]+a-bl

    if Y[0] <= -1:
        gx = bl*Y[0]-a+bl

    if Y[0] > -1 and Y[0] < 1:
        gx = a*Y[0]

    if yprim_k_helper == 0:
        YPRIM = Y[1]-gx

    if yprim_k_helper == 1:
        YPRIM = -bg*Y[1]-bg*Y[0]+B*math.cos(2*np.pi*xfn*t)

    return YPRIM

def RKK4(t):
    y1list = []
    YY1list = []
    y2list = []
    YY2list = []
    y3list = []
    YY3list=[]
    y4list = []
    YNEWlist = []
    for k in range(NN):      
        y1 = tstep*yprim(Y_in, k, t)
        y1list.append(y1)
    for k in range(NN):       
        YY1 = Y_in[k] + y1list[k]/2
        YY1list.append(YY1)
    for k in range(NN):        
        y2 = tstep*yprim(YY1list, k, t)
        y2list.append(y2)
    for k in range(NN):       
        YY2 = Y_in[k]+y2list[k]/2
        YY2list.append(YY2)
    for k in range(NN):       
        y3 = tstep*yprim(YY2list, k, t)
        y3list.append(y3)
    for k in range(NN):       
        YY3 = Y_in[k]+y3list[k]
        YY3list.append(YY3)
    for k in range(NN):        
        y4 = tstep*yprim(YY3list, k, t)
        y4list.append(y4)
    for k in range(NN):      
        YNEW = Y_in[k] + (y1list[k] + 2*y2list[k] + 2*y3list[k] + y4list[k])/6
        YNEWlist.append(YNEW)       
    return YNEWlist

xfnlist=[]
while xfn < xfnmax:
    for K1 in range(n2):
        t = 0
        for KK in range(NIT):
            rkk4list = RKK4(t)
            YNEW1.append(rkk4list[0])
            YNEW2.append(rkk4list[1])
            t = t + tstep
            Y_in = rkk4list
            xfnlist.append(xfn)

    xfn = xfn + dxfn
    if xfn >= xfnmax:
        break

final = [YNEW1, YNEW2]
plt.scatter(xfnlist, YNEW1, s=1, marker= ',')

axis = [xmin, xmax, ymin, ymax]
plt.axis(axis)
plt.show()


num=7
for i in range(2,num):
    if num%i==0:
        print('oxi')
        break
else:
    print('nai')
    


