# System info is the list of modified state vectors for each body in the system at a point in time. 
# The first entry is the mass of the body, the next 6 follow the form of a standard 3d state vector.

import numpy as np
from typing import Callable

def galilean_forces(state_vec:np.ndarray, system_info, t):
    # Extracting the state vectors for the system bodies from system_info
    system_states = [entry[-6:] for entry in system_info]
    body_masses = [entry[0] for entry in system_info]
    G = 6.67430e-20 # km^3 kg^-1 s^-2

    for state_vecO in system_states:
        if state_vec != state_vecO:
            r_vec:np.ndarray = np.array([state_vec[0:3] - state_vecO[0:3]])
            # r_mag = (r_vec[0]**2 + r_vec[1]**2 + r_vec[2]**2)**0.5 #! I used np.linalg.norm() instead
            r_mag:np.ndarray = np.linalg.norm(r_vec) # magnitude of r vector
            force_mag = G * body_masses[system_states.index(state_vec)] * body_masses[system_states.index(state_vecO)] / r_mag**2
            force_vec = [force_mag*r_vec[0] / r_mag, force_mag*r_vec[1] / r_mag, force_mag * r_vec[2] / r_mag]
            state_vec[7] += force_vec[0]
            state_vec[8] += force_vec[1]
            state_vec[9] += force_vec[2]


def RK4(f:Callable, state_vec0:np.ndarray, t0:int, tf:int, h:int):
    t = t0
    state_vec:np.ndarray = state_vec0
    step_count = (tf - t0) / h 
    for i in range(int(step_count)):
        k1 = f(state_vec, t)
        k2 = f(state_vec + 0.5 * h * k1, t + 0.5 * h)
        k3 = f(state_vec + 0.5 * h * k2, t + 0.5 * h)
        k4 = f(state_vec + h * k3, t + h)
        state_vec += (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += h
        # state_vec.append(state_vec) #! Why re-append the state vector? This becomse [r1,r2,r3,v1,v2,v3,[r1,r2,r3,v1,v2,v3]]
    return state_vec

def stormer_verlet(f, state_vec0, t0, tf, h):
    t:int = t0
    state_vec:list[float] = state_vec0
    step_count:float = (tf - t0) / h 
    x:list[float] = state_vec[0:3]
    v:list[float] = state_vec[3:6]
    for i in range(int(step_count)):
        v_half = v + 0.5 * h * f(x, t)
        x = x + h * v_half
        v = v_half + 0.5 * h * f(x, t + h)
        t += h
        state_vec.append([x[0], x[1], x[2], v[0], v[1], v[2]])
    return state_vec