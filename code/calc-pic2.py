import math
from pylab import mpl
import matplotlib.pyplot as plt
import numpy as np


EPS = 1e-15
PI = math.pi

h = 0.001
y = []
st = 10
L = 1.

def sgn(x):
    if(x>1e-20): return 1
    elif(x<-1e-20): return -1
    else: return 0

def get_y(x):   
    if(sgn(x-L)==0):return y[int(1./h)]
    id = int(x/h)
    v = y[id] * (x/h-id) + y[id+1] * (id+1-x/h)
    return v**2

def Calc(L,R):
    mid = (L + R) / 2
    return (R-L)*(get_y(L)+get_y(R)+get_y(mid)*4)/6

def Simpson(L,R):
    mid = (L + R) / 2
    ans = Calc(L,R)
    Lans = Calc(L,mid)
    Rans = Calc(mid,R)
    if(math.fabs(Lans + Rans - ans) < EPS): return ans
    return Simpson(L,mid)+Simpson(mid,R)

alpha = 10 # m*omega/(\overline{h})
def V_1(x):
    if(x<0 or x>L):return 1e200
    else: return (alpha*(x-L/2))**2

def V_0(x):
    return 0

def V(x):
    return V_0(x)

def f(E_n,id):
    return E_n - V(id*h)

def insert_pic():
    S = Simpson(0,L)
    x = []
    now = 0
    while(now <=  L + 1e-12):
        x.append(now)
        now += h
    # y2 = []
    for i in range(0,len(y)):
        y[i] /= math.sqrt(S)
        # y2.append(math.sqrt(2/L) * math.sin(5*PI/L * x[i]))
        # y[i] = y[i]**2
    xpoints = np.array(x)
    ypoints = np.array(y)
    x2points = np.linspace(0,1,1000000,endpoint = True)
    y2points = math.sqrt(2/L) * np.sin(5*PI/L * x2points)
    plt.plot(xpoints,ypoints,label="实验曲线")
    plt.plot(x2points,y2points,ls="--",label="理论曲线")

def Numerov(E_n):
    y.clear()
    y.append(0)
    y.append(st * h)
    i = h
    n = 1
    while(i + h <= L + 1e-12):
        a = (2-5*h*h*f(E_n,n)/6)*y[n] - (1+h*h*f(E_n,n-1)/12)*y[n-1]
        b = 1 + h*h*f(E_n,n+1) / 12
        y.append(a/b)
        n += 1
        i += h
    insert_pic()

mpl.rcParams['font.sans-serif'] = 'SimHei'

plt.figure(figsize = (12,6))
plt.xlim((0,1))

Numerov(246.74011)

plt.xlabel("位置$x$")
plt.ylabel("波函数$\psi$")
plt.title("基于Numerov方法和打靶法的数值曲线与理论曲线比较")
plt.legend()
plt.savefig("tmp.png")