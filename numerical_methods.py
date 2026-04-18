# System info is the list of modified state vectors for each body in the system at a point in time. 
# The first entry is the mass of the body, the next 6 follow the form of a standard 3d state vector.

import numpy as np
from typing import Callable
from math import ceil

def galilean_forces(
        state_vec:list[float],
        system_states:list[list[float]], 
        moon_mus:list[float],
        Jupiter_mu:float,
        moon_index:int,
)->list[float,float,float,float,float,float]:
    """
    Calculates the state derivative for the specified moon
    :param state_vec: The state vector of a specific moon
    :param system_states: The state vector of all the moons at a specific time
    :param masses: Masses of all the moons
    :returns: Returns the acclerationn state vector for the specified moon
    """

    state_derivative:np.ndarray = np.array([
        state_vec[3], # v1
        state_vec[4], # v2
        state_vec[5], # v3
        0,            # a1
        0,            # a2
        0             # a3
    ])

    state_derivative[3:6] += -Jupiter_mu * (state_vec[0:3] / (np.linalg.norm(state_vec[0:3])**3)) # Jupiter's contribution to the acceleration

    for i, state_vecO in enumerate(system_states):
        if i == moon_index:
            continue

        r_vec:np.ndarray = state_vec[0:3] - state_vecO[0:3]
        r_mag:np.ndarray = np.linalg.norm(r_vec) # magnitude of r vector
        acceleration_vec = -r_vec * moon_mus[i] / r_mag**3
        state_derivative[3:6] += acceleration_vec[0:3]
    
    return state_derivative

def states_time_hist(system_states0:list, moon_mus:list, Jupiter_mu:float, galilean_forces:Callable, num_method:Callable, t0:int, tf:int, h:int):
    """
    This function takes one of the integrators, the initial states of the moons, 
    the gravitational parameters of the moons, the initial and final times,
    and returns an ephemeris for every moon.
    :param states0: The initial states of the moons
    :param moon_mus: An np array of all the moons' gravitational parameters
    :param Jupiter_mu: The gravitational parameter of Jupiter
    :param f: The integrator function to use
    :param t0: time initial
    :param tf: time final
    :param h: Step size
    :returns: A list of all the states between the initial and final time
    """
    step_count_nom = (tf-t0)/h
    step_count = ceil(step_count_nom)

    state_ephemeris:list[list[list[float]]] = [system_states0];
    system_states = system_states0
    time_vec = [t0]
    t = t0
    for k in range(step_count):
        next_states = system_states.copy()
        for j, state_vec in enumerate(system_states):
            next_states[j] = num_method(state_vec, j, system_states, moon_mus, Jupiter_mu, h, galilean_forces)
        state_ephemeris.append(next_states)
        system_states = next_states
        time_vec.append(t+h)
        t = t + h

    return (state_ephemeris, time_vec)
        
    

def RK4(state_vec:np.ndarray, moon_index:int, system_states:list[list[float]], moon_mus:list[float], Jupiter_mu:float, h:int, f:Callable):
    """
    This function takes an initial state, timestep, inital time, and an 
    acceleration function
    :param x: Initial state in format [r1,r2,r3,v1,v2,v3]
    :param h: Time step in seconds and km
    :param t: Initial time
    :param f: Accelaration function
    :returns: Returns the next state
    """
    # Important note, this system is not explicitly time-dependent which is why time isn't being passed in as an input to the state's derivative.
    k1:np.ndarray = f(state_vec, system_states, moon_mus, Jupiter_mu, moon_index)
    k2:np.ndarray = f(state_vec+(h/2)*k1, system_states, moon_mus, Jupiter_mu, moon_index)
    k3:np.ndarray = f(state_vec+(h/2)*k2, system_states, moon_mus, Jupiter_mu, moon_index)
    k4:np.ndarray = f(state_vec+h*k3, system_states, moon_mus, Jupiter_mu, moon_index)
    state_next:np.ndarray = state_vec + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

    return state_next

def stormer_verlet(state_vec:np.ndarray, moon_index:int, system_states:list[list[float]], moon_mus:list[float], Jupiter_mu:float, h:int, f:Callable):
    x:list[float] = state_vec[0:3]
    v:list[float] = state_vec[3:6]
    accel1 = f(state_vec, system_states, moon_mus, Jupiter_mu, moon_index)[3:6]

    v_half:float = v + 0.5 * h * accel1
    x_next:float = x + h * v_half
    temp_state = np.array([x_next[0], x_next[1], x_next[2], v_half[0], v_half[1], v_half[2]])

    accel2 = f(temp_state, system_states, moon_mus, Jupiter_mu, moon_index)[3:6]

    v_next:float = v_half + 0.5 * h * accel2
    state_vec = np.array([x_next[0], x_next[1], x_next[2], v_next[0], v_next[1], v_next[2]])

    return state_vec
