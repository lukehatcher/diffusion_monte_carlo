import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt


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


def wfall_plots():
    edges = np.linspace(70, 170, num=11)
    max_occupance = 700000

    walkers = np.load("mil_prot_trimer_walks.npy")
    DW = np.load("mil_prot_trimer_dw.npy")
    a = get_vec(walkers, 10, 3)
    b = get_vec(walkers, 9, 3)
    thetas = get_angle(a, b)
    r_oh = la.norm(a, axis=1) * 0.529177

    for i in range(len(edges)-1):
        bools = (thetas > edges[i]) * (thetas < edges[i + 1])
        index = np.where(bools)[0]
        print(str(len(index)))
        p, bins = np.histogram(r_oh[index], bins=25, range=[.6, 1.6], density=False, weights=DW[index])  #density needs to change later
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        plt.rcParams.update({'font.size': 15})
        plt.tight_layout()
        plt.plot(bin_centers, p+i)
        plt.ylabel("$N_{W}$")
        plt.yticks([0, 2, 4, 6, 8], [75, 95, 115, 135, 155])
        plt.plot(bin_centers, (p/max_occupance)+i)
        plt.xlabel("$R_{OH}(\AA)$")
        plt.ylabel("$\\theta_{HOH}(^\circ)$")
    plt.savefig(fname="wfall_plot", dpi=500)
    plt.show()
    return None
#run = wfall_plots


def check_plane(ox_1, ox_2, ox_3, load_walkers):
    oo_1 = load_walkers[:, ox_1 - 1] - load_walkers[:, ox_3 - 1]
    oo_2 = load_walkers[:, ox_2 - 1] - load_walkers[:, ox_3 - 1]
    normal_vec_list = np.cross(oo_1, oo_2)
    print(normal_vec_list)
    return ox_1, ox_2, ox_3


walkers = np.load("protonated_clusters/100_prot_trimer_walks.npy")
edges2 = np.linspace(-1, 1, num=9)
DW = np.load("protonated_clusters/100_prot_trimer_dw.npy")

# run = check_plane(1, 2, 3, walkers)


# def o2_plane(ox1, ox3, load_walkers):
#     oo_vecs = load_walkers[:, ox1 - 1] - load_walkers[:, ox3 - 1]
#     r_oh = la.norm(a, axis=1) * 0.529177
#     z_comps = oo_vecs[:, 2]
#     print(z_comps)
#     return ox1, ox3


def z_of_hydrogen(H_n, load_walkers):
    hydrogens = load_walkers[:, H_n - 1]
    z_comps = hydrogens[:, 2] * 0.529177
    # print(z_comps)
    return z_comps


def plot_z(z_comps):
    roh_vec = get_vec(walkers, 9, 3)
    r_oh2 = la.norm(roh_vec, axis=1) * 0.529177
    for i in range(len(edges2)-1):
        bools = (z_comps > edges2[i]) * (z_comps < edges2[i + 1])
        index = np.where(bools)[0]
        print(str(len(index)))
        p, bins = np.histogram(r_oh2[index], bins=25, range=[.6, 1.6], density=False, weights=DW[index])
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        plt.rcParams.update({'font.size': 15})
        plt.tight_layout()
        # plt.plot(bin_centers, p + i)

        plt.plot(bin_centers, (p/25)+i)
        plt.yticks([0, 2, 4, 6], [-.875, -.375, .125, .625])
        plt.xlabel("$R_{OH}(\AA)$")
        plt.ylabel("$Z_{H}(\AA)$")
    # plt.savefig(fname="wfall_plot", dpi=500)
    plt.show()
    return None

run = z_of_hydrogen(9, walkers)
run2 = plot_z(run)


