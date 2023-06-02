import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import matplotlib.animation as animation


x_axis = np.linspace(-200,200,5000) #space grid: the only global variable of the code
deltax = x_axis[1]-x_axis[0]



def grid_init():
    '''This function returns the above initialized space grid, so that the array in the library will automatically match with the array in the user-interface scripts'''
    deltax = x_axis[1]-x_axis[0]
    return x_axis, deltax


def norm(phi):
    '''This method takes a wave function and normalize it'''
    norm = np.sum(np.square(np.abs(phi)))*deltax
    return phi/np.sqrt(norm)


def complex_plot(x_value, y_value, title_plot = '', save_option = False, name_file = ''):
    '''This method takes a complex array and plots all its component (real part, imaginary part, absolute value) as a function of a real valued array'''
    real = np.real(y_value)
    imag = np.imag(y_value)
    real_plot,*_ = plt.plot(x_value,real, ':', label='Re')
    im_plot,*_ = plt.plot(x_value,imag, ':', label='Im')
    norm_plot,*_ = plt.plot(x_value ,np.abs(y_value),label='$\sqrt{P}$')
    plt.title(title_plot)
    if save_option == True:
        plt.savefig(fname = name_file, dpi = 300)
    return real_plot, im_plot ,norm_plot
   
    
def wave_packet(pos=0, mom=0, sigma=2):
    '''This method takes the two parameters of a Gaussian distribution and returns the wave function of a particle with a certain momentum moduled as a Gaussian wave packet'''
    return (2 * np.pi * np.square(sigma))**(-1/4) * np.exp(1j*mom*x_axis)*np.exp(-np.square((x_axis-pos)/(np.square(sigma))),dtype=complex)


def double_x_derivative(wave_func):
    '''This method takes the wave function and returns its discretized second order derivative with respect to space coordinate'''
    second_derivative = -2*wave_func
    second_derivative[:-1] += wave_func[1:]
    second_derivative[1:] += wave_func[:-1]
    return second_derivative/np.square(deltax)


def hamiltonian(wave_func, h=1, m=1, potential=0):
    '''This method takes a wave function and a potential and returns the corresponding Hamiltonian. Physical parameters:
    -h:  reduced planck constant (or hbar), set to unity by default (float)
    -m:  mass of the particle, set to unity by default (float)'''
    return 1j*h/2/m * double_x_derivative(wave_func) - 1j*potential*wave_func/h


def runge_kutta(wave_func, dt, **kwargs):
    '''This method implements time evolution of a wave function according to Schroedinger equation for time value dt. In order to do so a fourth-order Runge-Kutta algorithm is used'''
    k1 = hamiltonian(wave_func, **kwargs)
    k2 = hamiltonian(wave_func+dt/2*k1, **kwargs)
    k3 = hamiltonian(wave_func+dt/2*k2, **kwargs)
    k4 = hamiltonian(wave_func+dt*k3, **kwargs)
    
    return wave_func + dt/6*(k1+2*k2+2*k3+k4)


def simulation(wave_func, potential=0, steps=100000, dt=0.003, orthogonal_to = 0, save_every=1000):
    '''This method simulates a wave function evolving in time according to Schroedinger equation. Condition for normalization and orthogonality check to some certain eigenstates are implemented. It returns an array with the values of wave function at each timestep and at each spacial point. Input variables:
    -steps: number of iterations of time integration (int)
    -dt: value of time differential. If the potential is too high, the wave function might collapse at zero everywhere. In that case, a smaller value for dt would help (float)
    -orthogonal_to: all eigenstates the input wave function will be forced to be orthogonal to at each timestep (python list)
    -save_every: the number of iterations after which the value of the wave function will be saved and stored (int)'''
    evolving_wave = np.zeros((int(steps/save_every), 5000), dtype = complex)
    
    if np.sum(orthogonal_to) != 0:    
        l = 0
        for i in range(steps):
            wave_func = runge_kutta(wave_func, dt, potential=potential)
            for k in range(len(orthogonal_to)):
                wave_func = wave_func - np.sum(np.conjugate(orthogonal_to[k])*wave_func)*deltax*orthogonal_to[k]
            wave_func = norm(wave_func)
            if i % save_every == 0:
                evolving_wave[l,:] = wave_func
                l += 1
    else:
        l = 0
        for i in range(steps):
            wave_func = runge_kutta(wave_func, dt, potential=potential)
            wave_func = norm(wave_func)
            if i % save_every == 0:
                evolving_wave[l,:] = wave_func
                l += 1
                  
    return evolving_wave



def pot_init(potential):
    '''This method is used for visualization of the potential in the animation'''
    plt.fill_between(x_axis, potential, -10,color='orange',alpha=0.2)

def animate(evolving_wave, potential, init_func = None, title = ''):
    '''This method takes the array of the wave function evolving in time and returns an animation of all components (real/imaginary values, probability) out of it'''
    fig = plt.figure()
    real_wave, imag_wave, probability = complex_plot(x_axis,evolving_wave[0], title_plot = title)
    plt.xlim(-50,50)
    plt.ylim(-2,2)
    plt.xlabel('Position')
    plt.ylabel(r'$\psi$')
    plt.legend()
    
    def animation(frame):
        '''This method is required to display different Line2D objects at each frame'''
        probability.set_data((x_axis, np.abs(evolving_wave[frame])))
        real_wave.set_data((x_axis, np.real(evolving_wave[frame])))
        imag_wave.set_data((x_axis, np.imag(evolving_wave[frame])))
        return probability, real_wave, imag_wave

    anim = FuncAnimation(fig, animation, frames=int(len(evolving_wave)), init_func = pot_init(potential),  interval = 50)
    
   
    #plt.show() 

    return anim


def transmission_coeff(wave_func, measuring_point):
    '''This method takes the wave function at a certain time and a point of space grid. It returns the probability of finding, at that time the particle before or after the desired point'''
    prob_before_barrier = np.sum(np.square(np.abs(wave_func)), where=(x_axis<measuring_point))*deltax
    prob_after_barrier = np.sum(np.square(np.abs(wave_func)), where=(x_axis>measuring_point))*deltax
    print('probability to find the wave after the barrier: ', prob_after_barrier)
    print('probability to find the wave before the barrier: ', prob_before_barrier)
    return prob_before_barrier, prob_after_barrier


def quantum_tunnelling(wave_func, V, measuring_point):
    '''This function takes an initial wave function, a potential (which should be suited for the tunnelling effect) and a point of the space grid. It simulates the wave function interacting with the potential, makes an animation of it and returns the probability of finding the particle before and after the input point'''
    phi = simulation(wave_func, potential = V, steps = 50000)
        
    anim = animate(phi, V, title = 'Quantum Tunnelling')
    gif_name = r"tunnelling_animation.gif" 
    writergif = animation.PillowWriter(fps=30, metadata=dict(artist='Me')) 
    anim.save(gif_name, writer=writergif)

    
    reflected_wave, transmitted_wave = transmission_coeff(phi[-1,:], measuring_point)
    
    return reflected_wave, transmitted_wave


def wave_coefficient(eigenstate, wave_func):
    '''This method evaluates the coefficient of a certain eigenstate'''
    return np.sum(np.conjugate(eigenstate)*wave_func)*deltax


def ground_state(wave_func, V):
    '''This method takes a wave function, a potential and returns the ground state eigenfunction with its corresponding coefficient. It also makes an animation of the eigenstate evolving in time (showing that only the phase will change).'''
    ground_eigenstate = simulation(wave_func, dt=-5e-3j, potential = V, steps = 100000, save_every = 5000)
    c_0 = wave_coefficient(ground_eigenstate[-1,:],wave_func)
    
#     GIF0 = animate(ground_eigenstate, V, title = 'Ground State formation')
#     gif_name = r"ground_state_formation.gif" 
#     writergif = animation.PillowWriter(fps=30) 
#     GIF.save(gif_name, writer=writergif)
    
    ground_eigenstate_evolving = simulation(ground_eigenstate[-1,:], potential = V, dt = 0.007, steps = 120000, save_every = 400)
    
    GIF = animate(ground_eigenstate_evolving, V, title = 'Ground State')
    gif_name = r"ground_state_animation.gif" 
    writergif = animation.PillowWriter(fps=30) 
    GIF.save(gif_name, writer=writergif)
    
    return ground_eigenstate[-1,:], c_0
    
def first_state(wave_func, ground_func, V):
    '''This method takes a wave function, a potential and returns the first excited eigenfunction with its corresponding coefficient. It also makes an animation of the eigenstate evolving in time (showing that only the phase will change).'''
    phi_1 = wave_func - np.sum(np.conjugate(ground_func)*wave_func)*deltax*ground_func
    first_eigenstate = simulation(phi_1, dt=-5e-3j, orthogonal_to = [ground_func], potential = V, steps=100000,save_every=5000)
    c_1 = wave_coefficient(first_eigenstate[-1,:], wave_func)
    
#     GIF0 = animate(first_eigenstate, V, title = 'First excited State formation')
#     gif_name = r"first_state_formation.gif" 
#     writergif = animation.PillowWriter(fps=30) 
#     GIF0.save(gif_name, writer=writergif)
    
    first_eigenstate_evolving = simulation(first_eigenstate[-1,:], potential = V, dt = 0.007, steps = 120000, save_every = 400)

    GIF = animate(first_eigenstate_evolving, V, title = 'First excited State')
    gif_name = r"first_state_animation.gif" 
    writergif = animation.PillowWriter(fps=30) 
    GIF.save(gif_name, writer=writergif)
    
    return first_eigenstate[-1,:], c_1

                           
def second_state(wave_func, ground_func, first_func, V):
    '''This method takes a wave function, a potential and returns the second excited eigenfunction with its corresponding coefficient. It also makes an animation of the eigenstate evolving in time (showing that only the phase will change).'''
    phi_0 = np.sum(np.conjugate(ground_func)*wave_func)*deltax*ground_func
    
    phi_1 = np.sum(np.conjugate(first_func)*wave_func)*deltax*first_func
    
    phi_2 = wave_func - phi_0 - phi_1 
    
    second_eigenstate = simulation(phi_2, dt=-5e-3j, orthogonal_to = [ground_func, first_func], potential = V, steps=80000,save_every=5000)
    c_2 = wave_coefficient(second_eigenstate[-1,:], wave_func)
    
#     GIF0 = animate(second_eigenstate, V, title = 'Second excited State formation')
#     gif_name = r"second_state_formation.gif" 
#     writergif = animation.PillowWriter(fps=30) 
#     GIF0.save(gif_name, writer=writergif)

    second_eigenstate_evolving = simulation(second_eigenstate[-1,:], potential = V, dt = 0.007, steps = 120000, save_every = 400)

    GIF = animate(second_eigenstate_evolving, V, title = 'Second excited State')
    gif_name = r"second_state_animation.gif" 
    writergif = animation.PillowWriter(fps=30) 
    GIF.save(gif_name, writer=writergif)
    
    return second_eigenstate[-1,:], c_2


def eigenstates(wave_func, V):
    '''this method takes a wave function and a potential. It calculates the first three eigenstates and print theire corrresponding coefficients. It also makes animations of the eigenstates evolving in time (showing that only the phase will change).'''
    ground_eigenstate, c_0 = ground_state(wave_func, V)
    print('Ground state calculated. The code is still running')
    first_eigenstate, c_1 = first_state(wave_func, ground_eigenstate, V)
    print('First excited state calculated. The code is still running')

    second_eigenstate, c_2 = second_state(wave_func, ground_eigenstate, first_eigenstate, V)
    
    print('c_0 =' + str(np.abs(c_0)))
    print('c_1 =' + str(np.abs(c_1)))
    print('c_2 =' + str(np.abs(c_2)))

                           
                           
                           
                           
    
    
    