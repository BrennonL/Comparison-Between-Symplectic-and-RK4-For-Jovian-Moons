import numpy as np
from horizons_interface import get_reference_ephemerides
from numerical_methods import (
    RK4,
    stormer_verlet,
    galilean_forces,
)

def main():
    (vectors_io, vectors_europa, vectors_ganymede, vectors_callisto) = get_reference_ephemerides()
    vectors_io_init = vectors_io[0]
    vectors_europa_init = vectors_europa[0]
    vectors_ganymede_init = vectors_ganymede[0]
    vectors_callisto_init = vectors_callisto[0]
    print("io ==", vectors_io_init)
    print("europa == ",vectors_europa_init)
    print("ganymede == ", vectors_ganymede_init)
    print("callisto == ", vectors_callisto_init)
    
    # Run RK4 for the moons:
    io_states_rk4:list[np.ndarray] = RK4(x=vectors_io_init, h=1, t=0, f=galilean_forces)
    europa_states_rk4:list[np.ndarray] = RK4(x=vectors_europa_init, h=1, t=0, f=galilean_forces)
    ganymede_states_rk4:list[np.ndarray] = RK4(x=vectors_ganymede_init, h=1, t=0, f=galilean_forces)
    callisto_states_rk4:list[np.ndarray] = RK4(x=vectors_callisto_init, h=1, t=0, f=galilean_forces)

    # stormer_verlet for the moons:
    io_states_SV:list[np.ndarray] = stormer_verlet(x=vectors_io_init, h=1, t=0, f=galilean_forces)
    europa_states_SV:list[np.ndarray] = stormer_verlet(x=vectors_europa_init, h=1, t=0, f=galilean_forces)
    ganymede_states_SV:list[np.ndarray] = stormer_verlet(x=vectors_ganymede_init, h=1, t=0, f=galilean_forces)
    callisto_states_SV:list[np.ndarray] = stormer_verlet(x=vectors_callisto_init, h=1, t=0, f=galilean_forces)


    # Run 
    # system_states = [entry[-6:] for entry in system_info]
    # body_masses = [entry[0] for entry in system_info]
# Need to define the architecture of the program
# Before I can start really coding stuff.

if __name__ == "__main__":
    main()