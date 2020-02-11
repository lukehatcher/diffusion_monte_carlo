import numpy as np


# class MyEquation:
#     @staticmethod
#     def equation1(a, b):
#         soln = a**2 + b**2
#         return soln


# a = MyEquation.equation1(2, 3)
# print(a)


# breh = np.genfromtxt(fname="h2o_dimer_gaussian_summary")

# breh2 = np.loadtxt(fname="h2o_dimer_gaussian_summary", usecols=(1, 2, 3, 4, 5), max_rows=100)

# unpack=True

# breh3 = np.loadtxt(fname="summary_test2.rtf", usecols=(1, 2, 3, 4, 5), unpack=True)


class MyClass:
    """
    parse txt file copy of Summary of Optimized Potential Surface Scans from gaussian
    """
    # @staticmethod
    def __init__(self,
                 filename,
                 ):
        # nsteps
        self.filename = filename

    def method(self, n_steps):
        # f = open("h2o_dimer_gaussian_summary", "r")
        with open(self.filename, "r") as gauss_data:
            ar = [line.split() for line in gauss_data]
            e_list = []
            r3oh_list = []
            r2oo_list = []
            for i in ar:
                if i[0] == "Eigenvalues":
                    energy = i[2:]
                    z = energy[0]  # only one value\
                    energy2 = z.split('-')
                    del energy2[0]  # removes weird first term from bad delimiter
                    e_list.append(energy2)
                elif i[0] == "r4":
                    r3oh_list.append(i[1:])
                elif i[0] == "r2":
                    r2oo_list.append(i[1:])

            """ flatten arrays """
            flat_elist = [x for list in e_list for x in list]
            flat_r3oh_list = [x for list in r3oh_list for x in list]
            flat_r2oo_list = [x for list in r2oo_list for x in list]

            """ sort flat arrays """
            oo_np = np.array([float(x) for x in flat_r2oo_list])
            oh_np = np.array([float(x) for x in flat_r3oh_list])
            eng_np = np.array([float(x) for x in flat_elist])
            sorted_oo_pos = np.lexsort((oh_np, oo_np))

            new_oo = oo_np[sorted_oo_pos]
            new_oh = oh_np[sorted_oo_pos]
            new_eng = eng_np[sorted_oo_pos]

            """ separate into oo steps and save """
            a = n_steps - n_steps
            b = n_steps
            np.save(arr=flat_r3oh_list[a:b], file="dimer_r3oh")
            for i in range(1, n_steps):
                np.save(arr=flat_elist[a:b], file="dimer_Es_" + str(i))
                np.save(arr=flat_r2oo_list[a:b], file="dimer_r2oo_" + str(i))
                a += n_steps
                b += n_steps

        return ar, e_list, r3oh_list, r2oo_list


a = MyClass(filename="h2o_dimer_gaussian_summary")
b = a.method(16)



# if __name__ == '__main__':
