import numpy as np
import matplotlib.pyplot as plt

wn = 219474.6
#kinetic energy
# h-bar = 1
# m = 1

def kinetic_energy(dvr_grid):
    t_matrix = np.zeros((len(dvr_grid), len(dvr_grid)))
    for i in range(len(dvr_grid)):
        for j in range(len(dvr_grid)):
            if i == j:
                t_matrix[i,j] = (np.pi**2) / (6*(delta_x**2))
            else:
                t_matrix[i,j] = ((-1)**(i-j)) / (((i-j)**2)*(delta_x)**2)
    return t_matrix

#harmonic osc. function / potential


def harmonic_oscillator(dvr_grid):
    w = 3000.0 /wn
    v_values = 0.5*(w**2)*(np.square(dvr_grid))
    return v_values


def potential_energy(dvr_grid):
    a = np.loadtxt("filename")
    v_matrix = np.diag(v_values)
    return v_matrix
    #returns diagonal Vmatrix


# #define dvr_grid and delta_x
# mydvr_grid = np.linspace(-100, 100, 1000)
# delta_x = (mydvr_grid[1] - mydvr_grid[0])
#
# #call
# real_v_matrix = potential_energy(mydvr_grid)
# real_kin_matrix = kinetic_energy(mydvr_grid)
#
# # plt.matshow(real_kin_matrix)
# # plt.colorbar()
# # plt.show()
# #
# # plt.matshow(real_v_matrix)
# # plt.colorbar()
# # plt.show()
#
# e, wfns = np.linalg.eigh(real_kin_matrix + real_v_matrix)
# print(e*wn,wfns)
# # plt.plot(mydvr_grid,np.diag(real_v_matrix))
# plt.plot(mydvr_grid,wfns[:,0]**2)
# plt.show()

eq_h7o3 = np.load("trimer_eq_geom_oo_fixed_2.npy")  # 2 corrected for ordering
# print(eq_h7o3)
ox_shifts = np.full(10, .25)  # bohr
ox_shifts[0] = 0

#############################

hyd_shifts = np.linspace(1.95725269-1, 1.95725269 + 1.25, 101)  # bohr


def get_oh_distance(eq_h7o3, atom_n1, atom_n2):
    oh_vector = eq_h7o3[atom_n1, :] - eq_h7o3[atom_n2, :]
    return oh_vector

# run = get_oh_distance(eq_h7o3, 3, 0)
# print(run)
# [ 1.95725269 -0.03664311 -0.10535852]

# def get_oo_distance(q_h7o3):


def get_geometries(eq_h7o3_load, ox_shifts_load, hyd_shifts_load, ox_1, ox_2, hyd_1, hyd_pull1, hyd_pull2):
    # geoms = np.zeros((100, 3))
    eq_h7o3_load[hyd_1, 1:] = 0  # y & z = 0
    eq_h7o3_load[ox_2, 0] -= (4 * .25) # shift ox back
    eq_h7o3_load[hyd_pull1, 0] -= (4 * .25)  # shift ox back
    eq_h7o3_load[hyd_pull2, 0] -= (4 * .25)  # shift ox back

    start_array = np.copy(eq_h7o3_load)
    start_array[hyd_1, 0] = hyd_shifts[0]-(hyd_shifts[1]-hyd_shifts[0])
    ct = 1
    ox_step_geoms_i = np.copy(start_array)
    for i in ox_shifts_load:
        ox_step_geoms_i[ox_2, 0] += i
        ox_step_geoms_i[hyd_pull1, 0] += i
        ox_step_geoms_i[hyd_pull2, 0] += i
        for j in hyd_shifts_load:
            h_step_geoms_i_subj = np.copy(ox_step_geoms_i)
            h_step_geoms_i_subj[hyd_1, 0] = j
            if j == hyd_shifts_load[0]:
                geoms = np.stack((ox_step_geoms_i, h_step_geoms_i_subj), axis=0)
            else:
                geoms = np.concatenate((geoms, h_step_geoms_i_subj.reshape((1, 10, 3))), axis=0)

        np.save(file="dvr_oo_geoms_shift2_" + str(ct), arr=geoms)
        ct += 1
    return geoms


run = get_geometries(eq_h7o3, ox_shifts, hyd_shifts, 3, 6, 2, 4, 5)

# print(run)
# print(np.shape(run))

# a = np.load("dvr_oo_geoms_shift_1.npy")
# b = np.load("dvr_oo_geoms_shift_2.npy")
# c = np.load("dvr_oo_geoms_shift_3.npy")
# print(a)
# print(np.shape(a))
# print(b)
# print(c)

# def writeXYZ(cds,num):
#     atms = ['H', 'H', 'H', 'O', 'H', 'H', 'O', 'H', 'H', 'O']
#     fll = open('scanPots/coord_'+str(num)+'.xyz', 'w+')
#     cdsA = cds * 0.529177
#     for i in range(len(cdsA)):
#         fll.write('%d\nasdf\n' % len(atms))
#         for k in range(len(atms)):
#             fll.write('%s %f %f %f\n' % (atms[k],cdsA[i,k,0],cdsA[i,k,1],cdsA[i,k,2]))
#         fll.write('\n')
#     fll.close()





