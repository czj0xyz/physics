import math
from pylab import mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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

def get_Area():
    ret = 0
    for i in range(0,len(y)-1):
        ret += h * ((y[i] + y[i+1])/2)**2
    return ret

def V_2(x):
    if(x<0 or x>L): return 1e200
    elif(x > L*0.2 or x < L*0.8): return 50
    else: return 0

alpha = 10 # m*omega/(\overline{h})
def V_1(x):
    if(x<0 or x>L):return 1e200
    else: return (alpha*(x-L/2))**2

def V_0(x):
    return 0

def V(x):
    return V_1(x)

def f(E_n,id):
    return E_n - V(id*h)

def insert_pic(E):
    S = get_Area()
    x = []
    now = 0
    while(now <=  L + 1e-12):
        x.append(now)
        now += h
    for i in range(0,len(y)):
        y[i] /= math.sqrt(S)
        # y[i] = y[i]**2
    
    xpoints = np.array(x)
    ypoints = np.array(y)
    plt.cla()
    plt.figure(figsize = (12,6))
    plt.xlim((0,L))
    plt.plot(xpoints,ypoints,label="$n=$"+str(len(ans)-1))
    plt.xlabel("位置$x$")
    plt.ylabel("波函数$\psi$")
    plt.title("基于Numerov方法和打靶法的一维谐振子的数值解")
    plt.legend()
    plt.savefig("pic4-" + str(len(ans)-1) + ".png")

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
    return y[n]

ans = []

def shooting(st_len):
    ans.clear()
    length=st_len
    for i in range(0,10000):
        lx = Numerov(length*i/10000)
        rx = Numerov(length*(i+1)/10000)
        flg = 0

        if(sgn(lx) == 0 and i!=0): 
            flg = 1
            ans.append(lx)
        
        if(sgn(rx) == 0):
            flg = 1
            ans.append(rx)
        
        if(sgn(lx)*sgn(rx)<0 and ~flg):
            flg = 1
            l = length/10000*i
            r = length/10000*(i+1)
            for _ in range(0,100):
                mid = (l+r)/2
                if(sgn(Numerov(mid))*sgn(rx) <= 0): l = mid
                else: r = mid
            if(sgn(l) == 1):
                ans.append(l)

        if(flg):
            insert_pic(l)
            print("OK")
            print(ans[len(ans)-1])


mpl.rcParams['font.sans-serif'] = 'SimHei'

L = 9
alpha = 3
st = 10
h = 0.01
st_len = 30.

shooting(st_len)
print(len(ans))