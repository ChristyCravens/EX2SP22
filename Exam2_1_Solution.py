from scipy.integrate import odeint
from math import sin
import numpy as np
from matplotlib import pyplot as plt

def odeSystem(X, t, *args):
    """
    this is the odeSystem callback I'm using for odeint().  The *args contains a callback for the input voltage.
    :param X: the current values of the state variables
    :param t: the current time from odeint
    :param args: fn(t), L, R, C
    :return: list of derivatives of state variables
    """
    #unpack X into convenient variables
    fn, L,R,C=args
    #assign friendly variable names
    i1=X[0]
    i2=X[1]
    #calculate the current input voltage
    vt=fn(t)
    #calculate derivatives for the state variables
    i1dot=(vt-(i1-i2)*R)/L
    # i1dot=(-vt+(i1-i2)*R)/L
    i2dot=i1dot-i2/(R*C)
    return [i1dot, i2dot]

def main():
    """
    For solving problem 1 on exam.  Note:  I'm passing a callback through odeint to model in input voltage
    :return:
    """
    vin=lambda t:  20*sin(20*t)  # a callback to be passed as an argument to odeSystem
    L=20; R=10; C=0.05
    myargs=(vin, L, R, C)  # callback, L, R, C
    x0=[0,0]  # initial conditions I[0]=0, I[1]=0
    tList=np.linspace(0,10,500)  # time vector
    I=odeint(odeSystem,x0,tList,myargs)  #solve using odeint

    #the plotting part
    plt.plot(tList,I[:,0], linestyle='solid', color='k', label=r'$i_1(t)$')
    plt.plot(tList,I[:,1], linestyle='dashed', color='k', label=r'$i_2(t)$')
    ax=plt.gca()
    ax.set_xlim(0,10)
    ax.set_ylim(-0.1,0.1)
    ax.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=12)  # format tick marks
    ax.tick_params(axis='both', grid_linewidth=1, grid_linestyle='solid', grid_alpha=0.5)
    ax.tick_params(axis='both', which='minor')
    ax.set_xticks([x for x in range(11)])
    ax.grid(True)
    ax.set_xlabel('t (s)')
    ax.set_ylabel(r'$i_1, i_2 (A)$')
    ax1=ax.twinx()
    ax1.plot(tList,R*(I[:,1]-I[:,0]), linestyle='dotted', color='k', label=r'$v_c(t)$')
    ax1.set_ylim(-0.45,0.1)
    ax1.tick_params(axis='y', which='both', direction='in', top=True, right=True, labelsize=12)  # format tick marks
    ax.legend()
    ax1.legend(loc='lower right')
    ax1.set_ylabel(r'$V_c(t) (V)$')
    plt.show()
    pass

main()