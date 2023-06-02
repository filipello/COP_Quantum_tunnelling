import quantum_functions as fn
import numpy as np

def potential(x_axis):
    '''In this function the user can insert any potential they want. Since this script was intentionally created to study the eigenstates of a constrained particle, please insert a potential which reaches the value infinity at the extremes of space grid. As an example, a quadratic potential well is already written.'''
    a = 0.0002
    
    return a*np.square(x_axis)



x_axis, deltax = fn.grid_init()

V = potential(x_axis)


initialisation = True

while initialisation == True:
    print('This is a program to simulate the first 3 eigenstates of a 1-D wave function')
    print('ATTENTION: in order to insert/modify the potential you have to change it in the script (have a look at READme file)')
    
    print('Input the initial position of the wave packet (inside the range of the space grid) \n')
    get_initial_position = float(input())
    print('Input the momentum of the wave packet (suggested values: 0.5-5.0) \n')
    get_momentum = float(input())
    print('Input the variance of the wave packet (suggested values: 1.5-2.5) \n')
    get_variance = float(input())
    
    psi = fn.wave_packet(pos = get_initial_position, mom = get_momentum, sigma = get_variance)
    
    print('The code is running...')
    
    fn.eigenstates(psi, V)
    break
    