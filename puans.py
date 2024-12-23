import matplotlib.pyplot as plt
from lmfit.models import QuadraticModel
import numpy as np
import math


def rotate_and_fit(a_start, b_start, data0):
    data = np.copy(data0)
    
    x = data[:, -2]
    y = data[:, -1]
    
 
    def chi_red(data_chi, angle_deg_chi):
        x = data[:, -2]*math.cos(angle_deg_chi*math.pi/180) - data[:, -1]*math.sin(angle_deg_chi*math.pi/180) #поворот по x
        y = data[:, -2]*math.sin(angle_deg_chi*math.pi/180) + data[:, -1]*math.cos(angle_deg_chi*math.pi/180) #поворот по y
        model = QuadraticModel()
        pars = model.guess(y, x=x) 
        out = model.fit(y, pars, x=x) 
        return out.result.redchi

    phi = 0.5 * (1.0 + math.sqrt(5.0))
    eps = 0.01
    a = a_start
    b = b_start
    while (True):
        x1 = b - (b-a)/phi
        x2 = a + (b-a)/phi
        if chi_red(data1, x1) >= chi_red(data1, x2):
            a = x1
        else:
            b = x2
        if abs(b-a) < eps:
            break
        o = (a+b)/2.0
        
    x = data[:, -2]*math.cos(o*math.pi/180) - data[:, -1]*math.sin(o*math.pi/180)
    y = data[:, -2]*math.sin(o*math.pi/180) + data[:, -1]*math.cos(o*math.pi/180)

    #для снимков после FIB, k = 1.236
    k = 1
    #для ветвей вниз:
    #y = k*y - min(y)*k 
    #x = x - x[y.argmax()]

    #для ветвей вверх:
    y = -k*y + max(y)*k
    x = x - x[y.argmin()]
    
    model = QuadraticModel()
    pars = model.guess(y, x=x) 
    out = model.fit(y, pars, x=x)
    a_p = out.params['a'].value
    
    plt.figure(dpi=300)
    plt.rcParams["figure.figsize"] = 7,7
    ax1 = plt.subplot2grid((3,1), (0,0), rowspan=2)
    #ax1.grid(which='major', color = 'gray', linestyle =":")
    plt.plot(x, y, 'bo',color='#1f77b4ff', markersize=6)
    plt.plot(x, out.best_fit, color='#ff7f0eff', linewidth=3)
    plt.legend (('Data', 'Quadratic fit'), loc='lower left')
    plt.ylabel('y, μm', fontsize=12)
    plt.tick_params(labelsize=12)
    ax1.grid()
    ax2 = plt.subplot2grid((3,1), (2,0))
    ax2.grid(which='major', color = 'gray', linestyle =":")
    out.plot_residuals()
    plt.title(None)
    plt.ylabel(None)
    #plt.savefig('profile')
    
    rad = abs(1/(2*a_p))
    return "Radius:", round(rad, 2), "Standard error: ", round(math.sqrt(a_p), 3), "Angle of rotation:", round(o, 7), "Reduced chi square:", round(a_p, 7), "a:", out.params['a'].value, "b:",out.params['b'].value, "c:", out.params['c'].value

#границы поиска углов разворота
a = -30 
b = 30

data1 = np.loadtxt('')
print(rotate_and_fit(a,b,data1))