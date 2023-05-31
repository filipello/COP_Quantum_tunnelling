import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import matplotlib.animation as animation


x_axis = np.linspace(-200,200,5000)
deltax = x_axis[1]-x_axis[0]



def grid_init():
    deltax = x_axis[1]-x_axis[0]
    return x_axis, deltax


def norm(phi):
    norm = np.sum(np.square(np.abs(phi)))*deltax
    return phi/np.sqrt(norm)


def complex_plot(x_value, y_value, title_plot = '', save_option = False, name_file = ''):
    real = np.real(y_value)
    imag = np.imag(y_value)
    real_plot,*_ = plt.plot(x_value,real, ':', label='Re')
    im_plot,*_ = plt.plot(x_value,imag, ':', label='Im')
#     plt.xlim(-2,2)
    norm_plot,*_ = plt.plot(x_value ,np.abs(y_value),label='$\sqrt{P}$')
    plt.title(title_plot)
    if save_option == True:
        plt.savefig(fname = name_file, dpi = 300)
    return real_plot, im_plot ,norm_plot
   
    
def wave_packet(pos=0, mom=0, sigma=2):

    return (2 * np.pi * np.square(sigma))**(-1/4) * np.exp(1j*mom*x_axis)*np.exp(-np.square((x_axis-pos)/(np.square(sigma))),dtype=complex)


def double_x_derivative(wave_func):
    second_derivative = -2*wave_func
    second_derivative[:-1] += wave_func[1:]
    second_derivative[1:] += wave_func[:-1]
    return second_derivative/np.square(deltax)


def hamiltonian(wave_func, h=1, m=1, potential=0):
    return 1j*h/2/m * double_x_derivative(wave_func) - 1j*potential*wave_func/h


def runge_kutta(wave_func, dt, **kwargs):
    
    k1 = hamiltonian(wave_func, **kwargs)
    k2 = hamiltonian(wave_func+dt/2*k1, **kwargs)
    k3 = hamiltonian(wave_func+dt/2*k2, **kwargs)
    k4 = hamiltonian(wave_func+dt*k3, **kwargs)
    
    return wave_func + dt/6*(k1+2*k2+2*k3+k4)


def simulation(wave_func, potential=0, steps=100000, dt=0.003, orthogonal_to = 0, save_every=1000):
    
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
#     plt.axvline(x=0, ymin=0.5, ymax=1, linewidth=2, color='r')
    plt.fill_between(x_axis, potential, -10,color='orange',alpha=0.2)

def animate(evolving_wave, potential, init_func = None, title = ''):
    fig = plt.figure()
    real_wave, imag_wave, probability = complex_plot(x_axis,evolving_wave[0], title_plot = title)
    plt.xlim(-50,50)
    plt.ylim(-2,2)
    plt.xlabel('Position')
    plt.ylabel(r'$\psi$')
    if init_func:
        init_func()
    
    
    def animation(frame):
        probability.set_data((x_axis, np.abs(evolving_wave[frame])))
        real_wave.set_data((x_axis, np.real(evolving_wave[frame])))
        imag_wave.set_data((x_axis, np.imag(evolving_wave[frame])))
        return probability, real_wave, imag_wave

    anim = FuncAnimation(fig, animation, frames=int(len(evolving_wave)), init_func = pot_init(potential),  interval = 50)
    plt.show()
    plt.legend()
    

    return anim


def transmission_coeff(wave_func, measuring_point):
    prob_before_barrier = np.sum(np.square(np.abs(wave_func)), where=(x_axis<measuring_point))*deltax
    prob_after_barrier = np.sum(np.square(np.abs(wave_func)), where=(x_axis>measuring_point))*deltax
    print('probability to find the wave after the barrier: ', prob_after_barrier)
    print('probability to find the wave before the barrier: ', prob_before_barrier)
    return prob_before_barrier, prob_after_barrier


def quantum_tunnelling(wave_func, V, measuring_point):
    phi = simulation(wave_func, potential = V, steps = 50000)
    
    reflected_wave, transmitted_wave = transmission_coeff(phi[-1,:], measuring_point)
    
    GIF = animate(phi, V)
    gif_name = r"tunnelling_animation.gif" 
    writergif = animation.PillowWriter(fps=30) 
    GIF.save(gif_name, writer=writergif)
    
    return reflected_wave, transmitted_wave


def wave_coefficient(eigenstate, wave_func):
    return np.sum(np.conjugate(eigenstate)*wave_func)*deltax


def ground_state(wave_func, V):
    
    ground_eigenstate = simulation(wave_func, dt=-5e-3j, potential = V, steps = 100000, save_every = 5000)
    c_0 = wave_coefficient(ground_eigenstate[-1,:],wave_func)
    
    ground_eigenstate_evolving = simulation(ground_eigenstate[-1,:], potential = V, dt = 0.007, steps = 160000, save_every = 400)
    
    GIF = animate(ground_eigenstate_evolving, V, title = 'Ground State')
    gif_name = r"ground_state_animation.gif" 
    writergif = animation.PillowWriter(fps=30) 
    GIF.save(gif_name, writer=writergif)
    
    return ground_eigenstate[-1,:], c_0
    
def first_state(wave_func, ground_func, V):
    
    phi_1 = wave_func - np.sum(np.conjugate(ground_func)*wave_func)*deltax*ground_func
    first_eigenstate = simulation(phi_1, dt=-5e-3j, orthogonal_to = [ground_func], potential = V, steps=100000,save_every=5000)
    c_1 = wave_coefficient(first_eigenstate[-1,:], wave_func)
    
    first_eigenstate_evolving = simulation(first_eigenstate[-1,:], potential = V, dt = 0.007, steps = 160000, save_every = 400)

    GIF = animate(first_eigenstate_evolving, V, title = 'First excited State')
    gif_name = r"first_state_animation.gif" 
    writergif = animation.PillowWriter(fps=30) 
    GIF.save(gif_name, writer=writergif)
    
    return first_eigenstate[-1,:], c_1

                           
def second_state(wave_func, ground_func, first_func, V):
    
    phi_0 = np.sum(np.conjugate(ground_func)*wave_func)*deltax*ground_func
    
    phi_1 = np.sum(np.conjugate(first_func)*wave_func)*deltax*first_func
    
    phi_2 = wave_func - phi_0 - phi_1 
    
    second_eigenstate = simulation(phi_2, dt=-5e-3j, orthogonal_to = [ground_func, first_func], potential = V, steps=80000,save_every=5000)
    c_2 = wave_coefficient(second_eigenstate[-1,:], wave_func)

    second_eigenstate_evolving = simulation(second_eigenstate[-1,:], potential = V, dt = 0.007, steps = 160000, save_every = 400)

    GIF = animate(second_eigenstate_evolving, V, title = 'Second excited State')
    gif_name = r"second_state_animation.gif" 
    writergif = animation.PillowWriter(fps=30) 
    GIF.save(gif_name, writer=writergif)
    
    return second_eigenstate[-1,:], c_2


def eigenstates(wave_func, V):
    
    ground_eigenstate, c_0 = ground_state(wave_func, V)
    
    first_eigenstate, c_1 = first_state(wave_func, ground_eigenstate, V)
    
    second_eigenstate, c_2 = second_state(wave_func, ground_eigenstate, first_eigenstate, V)
    
    print('c_0 =' + str(np.abs(c_0)))
    print('c_1 =' + str(np.abs(c_1)))
    print('c_2 =' + str(np.abs(c_2)))

                           
                           
                           
                           
    
    
    