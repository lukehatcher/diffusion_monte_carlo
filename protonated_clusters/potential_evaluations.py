import numpy as np
import subprocess as sub


def write_file(filename, filename_save, atoms_list):
        oo_geoms = np.load(filename)
        # numb_geoms = len(np.vstack(oo_geoms))
        f = open(filename_save, "w")
        f.write(str(len(atoms_list)) + "\n" + str(len(oo_geoms)) + "\n")
        for i in range(len(oo_geoms)):
            one_hyd_shift = oo_geoms[i]
            element = 0
            for line in one_hyd_shift:
                f.write(f"{line[0]} {line[1]} {line[2]} {atoms_list[j]} \n")
                element += 1
        f.close()
        return None


atoms_list = ["H", "H", "H", "O", "H", "H", "O", "H", "H", "O"]
run = write_file2("dvr_oo_geoms_shift_1.npy", "coord.dat", atoms_list)


def get_dat_files(numb_oo_shifts, filename): #filename must end in _#
    for i in range(1, numb_oo_shifts + 1):
        a = np.load(filenames + str(i) + ".npy")
        b = write_file2()
    return None




def get_pot():
    sub.run("./getpot, cwd="?")
    np.loadtxt(hlgjf)[:,0]
    return
