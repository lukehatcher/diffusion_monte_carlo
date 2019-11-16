import numpy as np
import numpy.linalg as linalg
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess as sub
import sys
import os

if not os.path.exists("h2o_analysis_data"):
    os.makedirs("h2o_analysis_data")

wn = 4.55634e-6  # cm^-1 to amu
au = 1822.89  # amu to au
m_hyd = (1.00784*au)
m_ox = (15.999*au)
dtau = 5.0
alpha = 1/(2*dtau)
hartree_conv = 219474.6  # hartree to cm^-1 (== 1/wn)

term_walks = sys.argv[1]
term_run = sys.argv[2]

descendant_time = 50
timeSteps = 10000

convergence_time = .25
equilibrium_time = int(timeSteps*convergence_time)  # flexible, used for calculating zpe
number_of_wfns = 10  # generally want: descendant_time * number_of_wfns = 0.5 timesteps, flexible
initial_walkers = term_walks

angst = 0.529177
num_atoms = 3  # H2O
xyz = 3
Vref_array = []
coordinates = np.zeros((initial_walkers, num_atoms, xyz))  # 1000 pairs of 3x(xyz)


def equilibrium_cds(cds):
    """
    ('''return)
    :param cds: zero array which is overrode to get starting eq. geometries
    :type cds: np.ndarray
    :return: cds: h2o geometries
    :rtype: np.ndarray
    """
    h2o_eq_coords = np.array([[0.9578400, 0.0000000, 0.0000000],
                              [-0.2399535, 0.9272970, 0.0000000],
                              [0.0000000, 0.0000000, 0.0000000]]) / angst * 1.01
    for i in range(len(cds)):  # 1000
        cds[i] = h2o_eq_coords
    return cds


def random_displacement(cds_array):  # arbitrary argument here,
    """
    randomly displaces walkers  #option + Enter > Insert documentation string stub or ('''enter)
    :param cds_array:
    :type np.ndarray
    :return: cds_array
    :rtype: np.ndarray
     """

    dispH1 = np.random.normal(loc=0.0, scale=(dtau/m_hyd)**0.5, size=(len(cds_array), 3))
    dispH2 = np.random.normal(loc=0.0, scale=(dtau/m_hyd)**0.5, size=(len(cds_array), 3))
    dispO = np.random.normal(loc=0.0, scale=(dtau / m_ox) ** 0.5, size=(len(cds_array), 3))  # could reshape(1000, 1, 3) and then hstack

    disps = np.stack((dispH1, dispH2, dispO), axis=1)
    cds_array = cds_array + disps

    return cds_array


def get_potential(cds_array):
    """

    :param cds_array:
    :type cds_array: np.ndarray
    :return: cds_array
    :rtype: np.ndarray
    """
    the_length = len(cds_array)
    cds_array = np.reshape(cds_array, (len(cds_array)*num_atoms, 3))
    np.savetxt("PES_water/hoh_coord.dat", cds_array, header=str(the_length), comments="")
    sub.run("./calc_h2o_pot", cwd="PES_water")
    vAr = np.loadtxt("PES_water/hoh_pot.dat")

    return vAr


def vref_stuff(vAr):
    VR = np.average(vAr) - alpha*((len(vAr) - initial_walkers) / initial_walkers)

    return VR


def birth_or_death(vAr, VR, cds, arb_who_from):
    birth_list = []
    death_list = []

    for i in range(len(vAr)):
        if vAr[i] < VR:
            Pb = np.exp(-1 * (vAr[i] - VR) * dtau) - 1
            if np.random.random() < Pb:
                birth_list.append(i)
        elif vAr[i] > VR:
            Pd = 1 - np.exp(-1 * (vAr[i] - VR) * dtau)
            if np.random.random() < Pd:
                death_list.append(i)

    cds = np.concatenate((cds, cds[birth_list]), axis=0)  # axis?
    cds = np.delete(cds, death_list, axis=0)
    vAr = np.concatenate((vAr, vAr[birth_list]))
    vAr = np.delete(vAr, death_list, axis=0)
    arb_who_from = np.concatenate((arb_who_from, arb_who_from[birth_list]))
    arb_who_from = np.delete(arb_who_from, death_list)

    return vAr, cds, birth_list, death_list, arb_who_from


""" call variables  """

coords_to_save = []
who_from = np.arange(len(coordinates))
coordinates = equilibrium_cds(coordinates)
# y = np.arange(1000 + 50, 10000 + 51, 1000)
weights_spots2 = np.arange((timeSteps / number_of_wfns), timeSteps + 1, timeSteps/number_of_wfns)
# x = np.arange(1000, 10000, 1000)
cds_spots2 = np.arange((timeSteps / number_of_wfns) - descendant_time, timeSteps, timeSteps/number_of_wfns)
wfn_ct = 1
wt_ct = 1


""" call """

for i in range(timeSteps + 1):
    coordinates = random_displacement(coordinates)
    V_array = get_potential(coordinates)
    if i == 0:
        Vref = vref_stuff(V_array)

    if i in cds_spots2:
        coords_to_save = (np.copy(coordinates)) * angst
        np.save("h2o_analysis_data/" + str(term_walks) + "_" + str(term_run) + "_" + "wfn" + str(wfn_ct) + ".npy",
                coords_to_save * angst)
        weights = np.zeros(len(coords_to_save))
        who_from = np.arange(len(coords_to_save))
        wfn_ct += 1
    V_array, coordinates, birth_list, death_list, who_from = birth_or_death(V_array, Vref, coordinates, who_from)

    if i in weights_spots2:
        individuals, occurrence = np.unique(who_from, return_counts=True)
        weights[individuals] = occurrence
        np.save("h2o_analysis_data/" + str(term_walks) + "_" + str(term_run) + "_" + "weights" + str(wt_ct),
                weights)
        wt_ct += 1

    Vref = vref_stuff(V_array)
    Vref_array.append(Vref)

    """ save vref """

    if i == timeSteps:
        np.save("h2o_analysis_data/" + str(term_walks) + "_" + str(term_run) + "_varray", Vref_array)

    print(len(coordinates))  # tracking walker variation
    print(i)  # tracking progress
print("gg")

