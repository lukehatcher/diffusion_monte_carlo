import numpy as np
import matplotlib.pyplot as plt

wn = 4.55634e-6  # cm^-1 to amu
omega = 4000*wn
au = 1822.89  # amu to au
m = (1/2)*(1.00784*au)
dtau = 5.0
alpha = 1/(2*dtau)
hartree_conv = 219474.6  # hartree to cm^-1 (== 1/ wn)
timeSteps = 1500
descendant_time = 50
initial_walkers = 1000
angst = 0.529177
coordinates = np.zeros(initial_walkers)
Vref_array = []



def random_displacement(cds_array):  # arbitrary argument here,
    """randomly displaces walkers  #option + Enter > Insert documentation string stub

    :param cds_array:
    :type cds_array: np.ndarray
    :return: cds_array
    :rtype:
    """
    displacement_array = np.random.normal(0.0, (dtau/m)**0.5, len(cds_array))  # 3rd arg will change as tau propagates
    cds_array += displacement_array

    return cds_array


def get_potential(cds_array):
    vAr = 0.5*m*(omega**2)*(cds_array**2)
    return vAr


def vref_stuff(vAr):
    VR = np.average(vAr) - alpha*((len(vAr) - initial_walkers) / initial_walkers)
    return VR
#---------------------------------
#ideas for vectorization
# np.where(vAr < VR,  )
#     if np.logical_and(vAr[i] < VR, np.random.random() < Pb) = True:
#         birth_list.append(i)
#---------------------------------

def birth_or_death(vAr, VR, cds, arb_who_from):
    birth_list = []
    death_list = []

    for i in range(len(vAr)):

        if vAr[i] < VR:
            Pb = np.exp(-1 * (vAr[i] - VR) * dtau) - 1
            if np.random.random() < Pb:
                birth_list.append(i)  # used for descendant weighting as well

        elif vAr[i] > VR:
            Pd = 1 - np.exp(-1 * (vAr[i] - VR) * dtau)
            if np.random.random() < Pd:
                death_list.append(i)  # used for descendant weighting as well

    cds = np.concatenate((cds, coordinates[birth_list]))
    cds = np.delete(cds, death_list)
    vAr = np.concatenate((vAr, coordinates[birth_list]))
    vAr = np.delete(vAr, death_list)
    arb_who_from = np.concatenate((arb_who_from, arb_who_from[birth_list]))
    arb_who_from = np.delete(arb_who_from, death_list)

    return vAr, cds, birth_list, death_list, arb_who_from


""" call """


coords_1450 = []
who_from = np.arange(len(coordinates))

for i in range(timeSteps + 1):
    coordinates = random_displacement(coordinates)
    V_array = get_potential(coordinates)
    if i == 0:
        Vref = vref_stuff(V_array)

    if i == timeSteps - descendant_time:
        coords_1450 = np.copy(coordinates)
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
plt.title("Vref convergence")
plt.xlabel("time")
plt.ylabel("V ref")
plt.grid()
plt.show()

plt.scatter(coords_1450, weights)
plt.title("coords vs weights")
plt.xlabel("coordinate")
plt.ylabel("weight")
plt.grid()
plt.show()

psi, bins = np.histogram(coords_1450/angst, bins=15, range=(-2, 2), density=True, weights=np.ones(len(coords_1450)))  # normal psi weights give psi squared
bin_centers = 0.5*(bins[:-1]+bins[1:])

psi_sqrd_weighted, bins2 = np.histogram(coords_1450/angst, bins=15, range=(-2, 2), weights=weights, density=True)  # psi squared, bins2 never used
psi_sqrd_analytic = psi**2 / (np.sum((psi**2)*(bin_centers[1]-bin_centers[0])))  #normed


plt.plot(bin_centers, psi_sqrd_analytic, label="psi^2 ana")
plt.plot(bin_centers, psi_sqrd_weighted, label="psi^2 weighted")
plt.plot(bin_centers, psi, label="psi")
plt.legend()
plt.show()




