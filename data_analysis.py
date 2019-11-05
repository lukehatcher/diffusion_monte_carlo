import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as linalg

""" commented out section included in nD DMC for now """

# """ get zpe """
#
# Vref_array = Vref_array[500:]  # allow for convergence
# zpe = np.average(Vref_array)
# print(Vref_array)
# print("zpe = " + str(zpe * hartree_conv))
#
#
# """ plot Vref convergence """
#
# plt.plot(np.array(Vref_array) * hartree_conv)
# plt.title("$V_{ref}$ convergence")
# plt.xlabel("imaginary time")
# plt.ylabel("$V_{ref}$")
# plt.grid()
# plt.show()

""" project psi^2 onto OH bond """

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