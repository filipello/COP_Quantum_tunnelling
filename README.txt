Code written by: Annet Konings and Lorenzo Filipello
Computational Physics 2022/2023 course, third project: "Quantum tunnelling and eigenstates."
The purpose of the code is to simulate the quantum tunnelling effect of various potential barriers and find the eigenstates numerically, 
starting from a Guassian wave packet and evolving the wave function through a fourth order Runge-Kutta algorithm.
There are three files. The first one is quantum_functions.py, containing all the functions.
The file quantum_tunnelling.py is suited for simulating the the tunnelling effect. The file 
quantum_eigenstates.py is suited for simulating the eigenstates of a constrainded wave function.
In order to change the potential the user have to manually edit the method potential() in the last two mentioned scripts, in order to give 
more freedom. Guidelines are given in the docstring of the method potential() while examples can be found at the bottom of this file.
Finally, by running the files quantum_tunnelling.py and quantum_eigenstates.py, the user interface will
ask for initial parameters (suggested values will be given).
The scripts will calculate the phsyical observable (transmission/reflection coefficient or eigenfunction coefficient) and save
a GIF of the evolving simulation.
If the user want to save the GIF of the forming eigenstates, please uncomment the commented lines of code in the methods
ground_state(), first_state() and second_state(). Please note that this will substantially increase the running time. 



Here some examples of potentials, if the user wants to use them, just copy paste the indented lines inside the function potential()

Look-a-like Dirac's delta potential:
	
	a = x_axis[2700]
	V_0 = 100
	return np.where((x_axis == a), V_0, 0)

Gaussian potential well:
		
	a = 0.0005
        b = 10
        return np.exp(a*np.square(x_axis))- b*np.exp(-a*np.square(x_axis)) + b - 1

Mexican-hat potential well:

	a = 0.000005
	b = 0.015
	c = 10
	return 0.000005*np.square(np.square(x_axis)) - 0.015*np.square(x_axis) + c