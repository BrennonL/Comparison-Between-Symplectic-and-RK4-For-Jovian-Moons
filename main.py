import numpy as np
from horizons_interface import get_reference_ephemerides
# from numerical_methods import (
#     RK4,
#     stormer_verlet,
#     galilean_forces,
# )

def main():
    (vectors_io, vectors_europa, vectors_ganymede, vectors_callisto) = get_reference_ephemerides()
    print(vectors_io)
    vectors_io_init = vectors_io[0]
    vectors_europa_init = vectors_europa[0]
    vectors_ganymede_init = vectors_ganymede[0]
    vectors_callisto_init = vectors_callisto[0]
    print("io ==", vectors_io_init) #1
    print("europa == ",vectors_europa_init) #2
    print("ganymede == ", vectors_ganymede_init) #3
    print("callisto == ", vectors_callisto_init) #4
    


    # Run 
    # system_states = [entry[-6:] for entry in system_info]
    # body_masses = [entry[0] for entry in system_info]
# Need to define the architecture of the program
# Before I can start really coding stuff.

if __name__ == "__main__":
    main()