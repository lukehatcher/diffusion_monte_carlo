import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
#coordinates are in bohr

# def get_vector(load_walkers):
#     v1 = load_walkers[:, 10-1] - load_walkers[:, 3-1]
#     v2 = load_walkers[:, 10-1] - load_walkers[:, 3-1]
#     return v1, v2


def get_vec(load_walkers, atom_n1, atom_n2):
    vector = load_walkers[:, atom_n1 - 1] - load_walkers[:, atom_n2 - 1]
    return vector


def get_angle(v1, v2):
    # angle = np.arccos(np.dot(v1, v2) / (la.norm(v1)*la.norm(v2)))
    angle = np.arccos(np.sum(v1*v2, axis=1) / (la.norm(v1, axis=1)*la.norm(v2, axis=1)))
    angle = np.degrees(angle)
    return angle


def hist_theta(angle, load_DW, bns=25, rng=None):
    p, bins = np.histogram(angle, bins=bns, range=rng, weights=load_DW, density=True)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    plt.plot(bin_centers, p)
    plt.show()
    plt.savefig(fname="bond_plot")
    return None


walkers = np.load("mil_prot_trimer_walks.npy")
DW = np.load("mil_prot_trimer_dw.npy")
a = get_vec(walkers, 10, 3)
b = get_vec(walkers, 9, 3)
thetas = get_angle(a, b)
r_oh = la.norm(a, axis=1) * 0.529177
# hist_theta(theta, DW, rng=[80, 100])

edges = np.linspace(70, 170, num=11)  # float


for i in range(len(edges)-1):
    bools = (thetas > edges[i]) * (thetas < edges[i + 1])
    index = np.where(bools)[0]
    print(str(len(index)))
    p, bins = np.histogram(r_oh[index], bins=25, range=[.6, 1.6], density=False, weights=DW[index])  #density needs to change later
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
#    plt.ylim(0, 20)
    plt.plot(bin_centers, ((p/700000)+i))
plt.savefig(fname="wfall_plot")
plt.show()

#test for distribution of angles
#p2, bins2 = np.histogram(thetas, bins=12, range=[70, 170])
#bin_centers2 = 0.5 * (bins2[:-1] + bins2[1:])
#plt.plot(bin_centers2, p2)
#plt.show()
