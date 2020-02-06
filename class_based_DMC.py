import numpy as np

class oop_DMC:
    """
    object oriented diffusion monte carlo
    """
    def __init__(self,
                 initial_walks = 5000,
                 dtau = 5.0,
                 eq_coords = []
                 ):

# eq coords needs to be a np.ndarray

        self.my_eq_cds = my_eq_coords
        self.


        def get_eq_cds(cds):



            h2o_eq_coords = np.array([[0.9578400, 0.0000000, 0.0000000],
                                      [-0.2399535, 0.9272970, 0.0000000],
                                      [0.0000000, 0.0000000, 0.0000000]]) / angst * 1.01

            for i in range(len(cds)):  # 1000
                cds[i] = h2o_eq_coords
            return cds




# wn = 4.55634e-6  # cm^-1 to amu
# au = 1822.89  # amu to au
# m_hyd = (1.00784*au)
# m_ox = (15.999*au)
# dtau = 5.0
# alpha = 1/(2*dtau)
# hartree_conv = 219474.6  # hartree to cm^-1 (== 1/wn)
#
# convergence_time = .25
# equilibrium_time = int(param.num_tsteps * convergence_time)  # flexible, used for calculating zpe
# number_of_wfns = 10  # generally want: descendant_time * number_of_wfns = 0.5 timesteps, flexible
#
# angst = 0.529177
# num_atoms = 3  # H2O
# xyz = 3
# Vref_array = []
# coordinates = np.zeros((param.init_walks, num_atoms, xyz))  # 1000 pairs of 3x(xyz)
