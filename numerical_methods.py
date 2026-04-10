# System info is the list of modified state vectors for each body in the system at a point in time. 
# The first entry is the mass of the body, the next 6 follow the form of a standard 3d state vector.

import numpy as np
from typing import Callable


[[[iostate1],[iostate2],[iostate3]],[ganymedestate1],[ganymedestate2]]

[[io], [gany], [next]]

def galilean_forces(initial_states:list, masses:list):
    # Extracting the state vectors for the system bodies from system_info
    system_states = [entry[-6:] for entry in system_info]
    body_masses = [entry[0] for entry in system_info]
    G = 6.67430e-20 # km^3 kg^-1 s^-2
    state_derivative:np.ndarray = np.array([
        state_vec[7], # v1
        state_vec[8], # v2
        state_vec[9], # v3
        0,            # a1
        0,            # a2
        0             # a3
    ])
    for state_vecO in system_states:
        if state_vec != state_vecO:
            r_vec:np.ndarray = np.array([state_vec[0:3] - state_vecO[0:3]])
            r_mag:np.ndarray = np.linalg.norm(r_vec) # magnitude of r vector
            force_mag = G * body_masses[system_states.index(state_vec)] * body_masses[system_states.index(state_vecO)] / r_mag**2
            force_vec = [force_mag*r_vec[0] / r_mag, force_mag*r_vec[1] / r_mag, force_mag * r_vec[2] / r_mag]
            state_derivative[7] += force_vec[0]
            state_derivative[8] += force_vec[1]
            state_derivative[9] += force_vec[2]
    
    return state_derivative

def states_time_hist(states0:list, masses:list, f:Callable, t0:int, tf:int, h:int):
    """
    This function takes one of the integrators and  masses of the moons and 
      loops through the time history, applying the times to the integrator
    :param states0: The initial states of the moons
    :param masses: An np array of all the moons' masses
    :param f: 
    :param t0: time initial
    :param tf: time final
    :param h: Step size
    :returns: A list of all the states between the initial and final time
    """
    ...
    state_ephemeris = [states0];
    states = states0

    # Set up time step loop
    i = t0
    while i < tf:
        i += h
        
    

def RK4(x:np.ndarray, h:int, t:int, f:Callable):
    """
    This function takes an initial state, timestep, inital time, and an 
    acceleration function
    :param x: Initial state in format [r1,r2,r3,v1,v2,v3]
    :param h: Time step in seconds and km
    :param t: Initial time
    :param f: Accelaration function
    :returns: Returns the next state
    """
    k1:np.ndarray = f(x,t)
    k2:np.ndarray = f(x+(h/2)*k1, t+(h/2))
    k3:np.ndarray = f(x+(h/2)*k2, t+(h/2))
    k4:np.ndarray = f(x+h*k3, t+h)
    x_next:np.ndarray = x + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    return x_next

# def RK4(f:Callable, state_vec0:np.ndarray, t0:int, tf:int, h:int):
#     t:int = t0
#     state_vec:np.ndarray = state_vec0
#     step_count = (tf - t0) / h 
#     states:list = [state_vec0]
#     for i in range(int(step_count)):
#         k1 = f(state_vec, t)
#         k2 = f(state_vec + 0.5 * h * k1, t + 0.5 * h)
#         k3 = f(state_vec + 0.5 * h * k2, t + 0.5 * h)
#         k4 = f(state_vec + h * k3, t + h)
#         state_vec += (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
#         t += h
#         states.append(state_vec)
#     return states

def stormer_verlet(state_vec:np.ndarray, h:int, t:int, f:Callable):
    x:list[float] = state_vec[0:3]
    v:list[float] = state_vec[3:6]
    v_half:float = v + 0.5 * h * f(x, t)
    x:float = x + h * v_half
    v:float = v_half + 0.5 * h * f(x, t + h)
    state_vec = (np.array([x[0], x[1], x[2], v[0], v[1], v[2]]))

    return state_vec