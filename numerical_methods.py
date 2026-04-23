# System info is the list of modified state vectors for each body in the system at a point in time. 
# The first entry is the mass of the body, the next 6 follow the form of a standard 3d state vector.

import numpy as np
from typing import Callable
from math import ceil

def galilean_derivative(
        aug_state:np.ndarray,
        moon_mus:list[float],
        Jupiter_mu:float,
)->np.ndarray:
    """
    Calculates the full coupled state derivative for the Galilean moon system.
    :param aug_state: Flat augmented state vector in the form
        [x1, y1, z1, vx1, vy1, vz1, x2, ...]
    :param moon_mus: Gravitational parameters of the moons
    :param Jupiter_mu: Gravitational parameter of Jupiter
    :returns: Flat augmented derivative vector in the form
        [vx1, vy1, vz1, ax1, ay1, az1, vx2, ...]
    """
    system_states = aug_state.reshape(-1, 6)
    state_derivative = np.zeros_like(system_states)

    for i, state_vec in enumerate(system_states):
        r_vec = state_vec[0:3]
        v_vec = state_vec[3:6]

        state_derivative[i, 0:3] = v_vec

        acceleration = -Jupiter_mu * (r_vec / (np.linalg.norm(r_vec) ** 3))

        for j, other_state in enumerate(system_states):
            if i == j:
                continue

            relative_position = r_vec - other_state[0:3]
            relative_distance = np.linalg.norm(relative_position)
            acceleration += -moon_mus[j] * relative_position / (relative_distance ** 3)

        state_derivative[i, 3:6] = acceleration

    return state_derivative.reshape(-1)

def augmented_state_packer(system_states: list[list[float]]) -> np.ndarray:
    """This function creates a big augmented state matrix."""
    return np.array(system_states, dtype=float).reshape(-1)

def augmented_state_unpack(state_vec: np.ndarray) -> list[np.ndarray]:
    """This function turns the augmented state matrix back into the standard ephemeris form of list[list[float]]"""
    return [body_state.copy() for body_state in state_vec.reshape(-1, 6)]

def states_time_hist(system_states0:list, moon_mus:list, Jupiter_mu:float, galilean_derivative:Callable, num_method:Callable, t0:int, tf:int, h:int):
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

    augmented_state0 = augmented_state_packer(system_states0)
    state_ephemeris = [augmented_state_unpack(augmented_state0)]
    augmented_state = augmented_state0
    time_vec = [t0]
    t = t0
    for k in range(step_count):
        augmented_state = num_method(augmented_state, moon_mus, Jupiter_mu, h, galilean_derivative)
        state_ephemeris.append(augmented_state_unpack(augmented_state))
        t += h
        time_vec.append(t)

    return (state_ephemeris, time_vec)
        
    

def RK4(augmented_state:np.ndarray, moon_mus:list[float], Jupiter_mu:float, h:int, f:Callable):
    """
    This function takes an initial state, timestep, inital time, and an 
    acceleration function
    :param x: Initial state in format [r1,r2,r3,v1,v2,v3]
    :param h: Time step in seconds and km
    :param t: Initial time
    :param f: Accelaration function
    :returns: Returns the next state
    """
    # The Galilean moon system is autonomous, so the derivative depends only on state.
    k1:np.ndarray = f(augmented_state, moon_mus, Jupiter_mu)
    k2:np.ndarray = f(augmented_state + (h / 2) * k1, moon_mus, Jupiter_mu)
    k3:np.ndarray = f(augmented_state + (h / 2) * k2, moon_mus, Jupiter_mu)
    k4:np.ndarray = f(augmented_state + h * k3, moon_mus, Jupiter_mu)
    state_next:np.ndarray = augmented_state + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    return state_next

def stormer_verlet(augmented_state:np.ndarray, moon_mus:list[float], Jupiter_mu:float, h:int, f:Callable):
    system_states = augmented_state.reshape(-1, 6)
    positions = system_states[:, 0:3]
    velocities = system_states[:, 3:6]

    derivative = f(augmented_state, moon_mus, Jupiter_mu).reshape(-1, 6)
    accel1 = derivative[:, 3:6]

    v_half = velocities + 0.5 * h * accel1
    x_next = positions + h * v_half

    temp_system_state = np.hstack((x_next, v_half)).reshape(-1)
    accel2 = f(temp_system_state, moon_mus, Jupiter_mu).reshape(-1, 6)[:, 3:6]

    v_next = v_half + 0.5 * h * accel2
    next_state = np.hstack((x_next, v_next)).reshape(-1)

    return next_state
