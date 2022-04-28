import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

# These are our given values for the RLC circuit system
R = 10  # ohms
L = 20  # H
C = 0.05  # F

# Although t has not been identified yet, so adding this here would give an undefined value,
#    I still wanted to be sure to show the initial v(t) equation given from the problem.
# vt = 20*np.sin(20*t)

# Using KVLs from the loops in the circuit given, I end up with:
# Li1' + R(i1-i2)=v(t)
# R(i1'-i2') + (i2/C) = 0

# I then use these equations to solve for i1' and i2' to use in the code below.


# Defining the function needed to solve
def f(i,t):
    """
    This creates the ODEs and adds the values we need to solve for using i as a tuple, where t is
    representative of our time in seconds. It then returns the equations for the derivatives of i1 and i2
    to solve more easily using odeint within the main function.
    :param i: tuple used to solve for missing parts
    :param t: time in seconds
    :return: tuple with the equations for i1dot and i2dot
    """
    # Designate the values for i1 and i2 to call back later
    i1 = i[0]
    i2 = i[1]
    # Equations for i1dot and i2dot, solved using KVLs from the circuit
    i1dot = (((20*np.sin(20*t)) - R*(i1-i2))/L)
    i2dot = i1dot - (i2/(R*C))
    # Return the ODEs from above as a tuple to solve for i1 and i2 using odeint later
    return [i1dot, i2dot]

def main():
    """
    This function takes the above equations with the missing components and solves them, giving the graph
    with the plots for i1, i2, and Vc(t) all on the same graph, including proper legends and identifiers.
    :return:
    """
    tf = 10.0 # final value for t in seconds
    t = np.arange(0, tf, 0.02) # 500 points for each graph across the time axis

    # Initial conditions
    i0 = [0, 0]
    # Using odeint, we can determine the solution from the above equations of f
    i = odeint(f, i0, t)
    i1, i2 = i[:, 0], i[:, 1]
    # Equation for Vc(t)
    vct = R*(i2-i1)
    # Plot the values of i1, i2, and Vc(t)
    fig, ax1 = plt.subplots()
    # Plotting i1 with a solid black line
    ax1.plot(t, i1, '-', color='black')
    # Plotting i2 with a dashed black line
    ax1.plot(t, i2, '--', color='black')
    # Creating the legend for both currents
    ax1.legend(['I1(t)', 'I2(t)'], loc='lower left')
    # Plotting Vc(t) with a dotted black line
    ax2 = ax1.twinx()
    ax2.plot(t, vct, linestyle='dotted', color='black')
    # Label the axes as time, i1/i2, and Vc(t)
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('I1 and I2 (A)', color='black')
    ax2.set_ylabel('Vc(t) (V)', color='black')
    # Adding a second legend for Vc(t)
    ax2.legend(['Vc(t)'])
    # Title the plot for the RLC circuit system given
    plt.title('RLC Circuit Plot')
    # Plot grid lines
    ax1.grid(alpha=0.3)
    ax2.grid(alpha=0.3)
    plt.show()

    # For the record: It hurt my soul that you didn't let us make a pretty plot this time! Black and white
    # plots are too boring! I need to add my own pizazz!

if __name__ == "__main__":
    main()