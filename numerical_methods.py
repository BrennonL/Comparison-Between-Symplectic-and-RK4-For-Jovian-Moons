# System info is the list of modified state vectors for each body in the system at a point in time. 
# The first entry is the mass of the body, the next 6 follow the form of a standard 3d state vector.

import numpy as np
from typing import Callable
from constants import (
    MU_JUPITER_km
)
from math import ceil
from enum import Enum
# [[[iostate1],[iostate2],[iostate3]],[ganymedestate1],[ganymedestate2]]

# [[io], [gany], [next]]

class Moon(Enum):
    IO=1
    EUROPA=2

def galilean_forces(
        state_vec:list[float],
        system_states:list[list[float]], 
        moon_mus:list[float],
        Jupiter_mu:float,
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
        if state_vec != state_vecO:
            r_vec:np.ndarray = state_vec[0:3] - state_vecO[0:3]
            r_mag:np.ndarray = np.linalg.norm(r_vec) # magnitude of r vector
            acceleration_vec = -r_vec * moon_mus[i] / r_mag**3
            state_derivative[3:6] += acceleration_vec[0:3]
    
    return state_derivative

def states_time_hist(system_state0:list, moon_mus:list, Jupiter_mu:float, galilean_forces:Callable, num_method:Callable, t0:int, tf:int, h:int):
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

    state_ephemeris:list[list[list[float]]] = [[system_state0]];
    system_state = system_state0
    time_vec = [t0]
    t = t0
    for k in range(step_count):
        next_state = system_state.copy()
        for j, state_vec in enumerate(system_state):
            next_state[j] = num_method(state_vec, h, t, galilean_forces(state_vec, system_state, moon_mus, Jupiter_mu))
        state_ephemeris.append(next_state)
        system_state = next_state
        time_vec.append(t+h)

    return state_ephemeris
        
    

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