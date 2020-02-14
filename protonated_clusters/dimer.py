import numpy as np
import matplotlib.pyplot as plt


class MyClass:
    """
    parse txt file copy of Summary of Optimized Potential Surface Scans from gaussian
    """
    # @staticmethod
    def __init__(self,
                 filename,
                 ):

        self.filename = filename

    def method(self, n_steps):
        with open(self.filename, "r") as gauss_data:
            ar = [line.split() for line in gauss_data]
            e_list = []
            r4oh_list = []
            r2oo_list = []
            for i in ar:
                if i[0] == "Eigenvalues":
                    energy = i[2:]
                    z = energy[0]  # only one value
                    energy2 = z.split('-')
                    del energy2[0]  # removes weird first term from bad delimiter
                    # energy2 *= -1
                    e_list.append(energy2)
                elif i[0] == "r4":
                    r4oh_list.append(i[1:])
                elif i[0] == "r2":
                    r2oo_list.append(i[1:])

            """ flatten arrays """
            flat_elist = [x for list in e_list for x in list]
            flat_r4oh_list = [x for list in r4oh_list for x in list]
            flat_r2oo_list = [x for list in r2oo_list for x in list]

            """ sort flat arrays """
            oo_np = np.array([float(x) for x in flat_r2oo_list])
            oh_np = np.array([float(x) for x in flat_r4oh_list])
            eng_np = np.array([float(x) for x in flat_elist])
            eng_np *= -1.0  # fixed removed "delimiter"
            sorted_oo_pos = np.lexsort((oh_np, oo_np))

            new_oo = oo_np[sorted_oo_pos]
            new_oh = oh_np[sorted_oo_pos]
            new_eng = eng_np[sorted_oo_pos]

            """ separate by oo steps and save """
            a = n_steps - n_steps
            b = n_steps
            np.save(arr=new_oh[0:n_steps], file="dimer_r4oh")
            for i in range(1, n_steps + 1):
                np.save(arr=new_eng[a:b], file="dimer_Es_" + str(i))
                np.save(arr=new_oo[a:b], file="dimer_r2oo_" + str(i))
                a += n_steps
                b += n_steps

        return ar, e_list, r4oh_list, r2oo_list


# a = MyClass(filename="h2o_dimer_gaussian_summary")
# b = a.method(16)


class DVR:
    def __init__(self,
                 dvrgrid,
                 potentials,
                 mu,
                 ):

        """
        discrete variable representation
        :param dvrgrid: aka OH steps
        :type dvrgrid: np array
        :param potentials: for a single OO distance thus should have (nsteps) values in it
        :type potentials: np array
        :param mu: reduced mass, calc by hand prior UNITS=ATOMIC, thus *1822.89 for amu --> au
        :type mu: float
        """

        self.dvrgrid = dvrgrid
        self.potentials = potentials
        self.mu = mu
        self.deltax = self.dvrgrid[1] - self.dvrgrid[0]

    def pot_energy(self):
        v_matrix = np.diag(self.potentials)  # whats the units on the loaded self.pot?
        return v_matrix

    def kinetic_energy(self):
        t_matrix = np.zeros((len(self.dvrgrid), len(self.dvrgrid)))
        for i in range(len(self.dvrgrid)):
            for j in range(len(self.dvrgrid)):
                if i == j:
                    t_matrix[i, j] = (np.pi ** 2) / (6 * (self.deltax ** 2) * self.mu)
                else:
                    t_matrix[i, j] = ((-1) ** (i - j)) / (self.mu * ((i - j) ** 2) * (self.deltax**2))
        return t_matrix

    def run_dvr(self):
        v_matrix = self.pot_energy()
        kinetic_energy = self.kinetic_energy()
        evals, wfns_vecs = np.linalg.eigh(v_matrix + kinetic_energy)
        plt.plot(self.dvrgrid, wfns_vecs[:, 0] ** 2)
        plt.show()
        return None


""" standard OH grid for passing """
grid_list = np.load(file="dimer_r4oh.npy")
grid_arr = np.asarray(grid_list)


""" pots for a single OO distance for passing """
pots_list = np.load(file="dimer_Es_1.npy")
pots_arr = np.asarray(pots_list)


ob = DVR(grid_arr, pots_arr, 1728.3085005881399)
ob.run_dvr()



# if __name__ == '__main__':
