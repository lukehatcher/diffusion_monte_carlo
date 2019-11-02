import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt
import subprocess as sub

"""Harry Partridge and David W. Schwenke potential"""

wn = 4.55634e-6  # cm^-1 to amu
# omega = 4000*wn
au = 1822.89  # amu to au
m_hyd = (1.00784*au)
m_ox = (15.999*au)
dtau = 5.0
alpha = 1/(2*dtau)
hartree_conv = 219474.6  # hartree to cm^-1 (== 1/wn)
timeSteps = 1500
descendant_time = 50
initial_walkers = 100
angst = 0.529177
num_atoms = 3  # H2O
xyz = 3
Vref_array = []
coordinates = np.zeros((initial_walkers, num_atoms, xyz))  # 1000 pairs of 3x(xyz)


def equilibrium_cds(cds):
    h2o_eq_coords = np.array([[0.9578400, 0.0000000, 0.0000000],
                              [-0.2399535, 0.9272970, 0.0000000],
                              [0.0000000, 0.0000000, 0.0000000]]) / angst * 1.01
    for i in range(len(cds)):  # 1000
        cds[i] = h2o_eq_coords
    return cds


def random_displacement(cds_array):  # arbitrary argument here,
    """randomly displaces walkers  #option + Enter > Insert documentation string stub
    :param cds_array:
    :type cds_array: np.ndarray
    :return: cds_array
    :rtype:
     """

    dispH1 = np.random.normal(loc=0.0, scale=(dtau/m_hyd)**0.5, size=(len(cds_array), 3))
    dispH2 = np.random.normal(loc=0.0, scale=(dtau/m_hyd)**0.5, size=(len(cds_array), 3))
    dispO = np.random.normal(loc=0.0, scale=(dtau / m_ox) ** 0.5, size=(len(cds_array), 3))  # could reshape(1000, 1, 3) and then hstack

    disps = np.stack((dispH1, dispH2, dispO), axis=1)
    cds_array = cds_array + disps

    return cds_array


def get_potential(cds_array):
    the_length = len(cds_array)
    cds_array = np.reshape(cds_array, (len(cds_array)*num_atoms, 3))
    np.savetxt("PES_water/hoh_coord.dat", cds_array, header=str(the_length), comments="")
    sub.run("./calc_h2o_pot", cwd="PES_water")
    # sub.wait()
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
    vAr = np.delete(vAr, death_list,axis=0)
    arb_who_from = np.concatenate((arb_who_from, arb_who_from[birth_list]))
    arb_who_from = np.delete(arb_who_from, death_list)

    return vAr, cds, birth_list, death_list, arb_who_from


""" call """

coords_1450 = []
who_from = np.arange(len(coordinates))
coordinates = equilibrium_cds(coordinates)

for i in range(timeSteps + 1):
    coordinates = random_displacement(coordinates)
    V_array = get_potential(coordinates)
    if i == 0:
        Vref = vref_stuff(V_array)

    if i == timeSteps - descendant_time:
        coords_1450 = (np.copy(coordinates))*angst
        weights = np.zeros(len(coords_1450))
        who_from = np.arange(len(coords_1450))  # who_from = np.array([x for x in range(len(coordinates))])

    V_array, coordinates, birth_list, death_list, who_from = birth_or_death(V_array, Vref, coordinates, who_from)

    if i == timeSteps:
        individuals, occurrence = np.unique(who_from, return_counts=True)
        weights[individuals] = occurrence

    Vref = vref_stuff(V_array)
    Vref_array.append(Vref)

    print(len(coordinates))  # tracking walker variation
    print(i)  # tracking progress


""" get zpe """

Vref_array = Vref_array[500:]  # allow for convergence
zpe = np.average(Vref_array)
print(Vref_array)
print(zpe * hartree_conv)


""" plotting """

plt.plot(np.array(Vref_array) * hartree_conv)
plt.title("$V_{ref}$ convergence")
plt.xlabel("imaginary time")
plt.ylabel("$V_{ref}$")
plt.grid()
plt.show()


r_oh = linalg.norm(coords_1450[:, 0]-coords_1450[:, 2], axis=1)
psi_sqrd, bins = np.histogram(r_oh, bins=25, range=(.5, 2), weights=weights, density=True)
bin_centers = 0.5*(bins[:-1]+bins[1:])
plt.plot(bin_centers, psi_sqrd)
plt.title("$\\psi^2$ projection onto $r_{OH}$")
plt.xlabel("$r_{OH} (\\AA)$")
plt.ylabel("$\\psi^2$")
plt.grid()
plt.show()


""" expectation value of position <x> along r_OH """

exp_pos = np.average(r_oh, weights=weights)
print(exp_pos)
