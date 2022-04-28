from math import *
from Rankine_Solution import rankine

def main():
    '''
    A test program for rankine power cycles.
    R1 is a rankine cycle object that is instantiated for turbine inlet of saturated vapor.
    R2 is a rankine cycle object that is instantiated for turbine inlet of superheated vapor.
    :return: none
    '''
    R1=rankine(p_high=8000, p_low=8, eff_turbine=0.95, name='Rankine Cycle - Saturated at Turbine Inlet')
    R1.calc_efficiency()
    R1.plot_cycle_TS()

    R2=rankine(p_high=8000, p_low=8, t_high=500, eff_turbine=0.95, name='Rankine Cycle - Superheated at Turbine Inlet')
    R2.calc_efficiency()
    R2.plot_cycle_TS()

main()