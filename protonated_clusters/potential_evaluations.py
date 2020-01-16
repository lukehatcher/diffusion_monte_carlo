import numpy as np
import subprocess as sub

def write_file(numb_oo_shifts, numb_h_shifts, numb_atoms):
    for i in range(1, numb_oo_shifts + 1):
        a = np.load("dvr_oo_geoms_shift_" + str(i) + ".npy")
        a = np.reshape(a, ((numb_oo_shifts * numb_h_shifts * numb_atoms) + 10, 3))  # +10 for the start_arrays in dvr.py
        file = open("data.dat", "w+")
    return a, file


run = write_file(10, 10, 10)
# print(run)


def get_pot():
    sub.run("./getpot, cwd="?")
    return _
