import quantum_functions as fn
import numpy as np

def potential(x_axis):
    
    a_begin = 10
    a_end = 15
    V_0 = 5
    
    return np.where((x_axis>a_begin)&(x_axis<a_end),V_0,0)


x_axis, deltax = fn.grid_init()

V = potential(x_axis)


initialisation = True
while initialisation == True:
    print('This is a program to simulate quantum tunnelling of a 1-D wave function')
    print('ATTENTION: in order to insert/modify the potential you have to change it in the script (have a look at READme file)')
    
    print('Input the initial position of the wave packet (inside the range of the space grid) \n')
    get_initial_position = float(input())
    print('Input the momentum of the wave packet (suggested values: 0.5-5.0) \n')
    get_momentum = float(input())
    print('Input the variance of the wave packet (suggested values: 1.5-2.5) \n')
    get_variance = float(input())
    
    psi = fn.wave_pacekt(pos = get_initial_position, mom = get_momentum, sigma = get_variance)
    
    print('Input the point after which you want to measure the trasmitted wave \n')
    get_measuring_point = float(input())
    
    fn.quantum_tunnelling(psi, V, get_measuring_point)
    
    break
