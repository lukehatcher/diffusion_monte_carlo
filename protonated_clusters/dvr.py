import numpy as np
import matplotlib.pyplot as plt
import argparse

wn = 219474.6
au = 1822.89


eq_h7o3 = np.load("trimer_eq_geom_oo_fixed_2.npy")  # 2 corrected for ordering
ox_shifts = np.full(10, .25)  # bohr
ox_shifts[0] = 0
hyd_shifts = np.linspace(1.95725269-1, 1.95725269 + 1.25, 101)  # bohr # len = 101


# ----------------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument(
    "--dvr_run-number",
    default=10,
    type=int,
    dest='dvr_run_numb'
)
param = parser.parse_args()
#param.dvr_run_numb
# ----------------------------------------------------------------------------------------------------------------------


def get_oh_distance(eq_h7o3, atom_n1, atom_n2):
    oh_vector = eq_h7o3[atom_n1, :] - eq_h7o3[atom_n2, :]
    return oh_vector


def get_geometries(eq_h7o3_load, ox_shifts_load, hyd_shifts_load, ox_1, ox_2, hyd_1, hyd_pull1, hyd_pull2):
    # geoms = np.zeros((100, 3))
    eq_h7o3_load[hyd_1, 1:] = 0  # y & z = 0
    eq_h7o3_load[ox_2, 0] -= (4 * .25) # shift ox back
    eq_h7o3_load[hyd_pull1, 0] -= (4 * .25)  # shift ox back
    eq_h7o3_load[hyd_pull2, 0] -= (4 * .25)  # shift ox back

    start_array = np.copy(eq_h7o3_load)
    start_array[hyd_1, 0] = hyd_shifts[0]-(hyd_shifts[1]-hyd_shifts[0])  # need to account for this in my dvr as well
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

# run = get_geometries(eq_h7o3, ox_shifts, hyd_shifts, 3, 6, 2, 4, 5)


# ----------------------------------------------------------------------------------------------------------------------
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
# ----------------------------------------------------------------------------------------------------------------------


def potential_energy(filename):
    v_vals = np.loadtxt("pots_for_dvr/" + filename)
    v_vals /= wn
    v_mat = np.diag(v_vals)
    return v_vals, v_mat


def get_grid_and_dx(v_vals, hyd_shift_array):
    mydvr_grid = np.copy(hyd_shift_array)
    # mydvr_grid = np.linspace(hyd_shift_array[0], hyd_shift_array[-1], len(hyd_shift_array))
    mydelta_x = hyd_shift_array[1] - hyd_shift_array[0]
    return mydvr_grid, mydelta_x


def kinetic_energy(mydvr_grid, mydelta_x, m):
    t_matrix = np.zeros((len(mydvr_grid), len(mydvr_grid)))
    for i in range(len(mydvr_grid)):
        for j in range(len(mydvr_grid)):
            if i == j:
                t_matrix[i,j] = (np.pi**2) / (6*(mydelta_x**2)*m)
            else:
                t_matrix[i, j] = ((-1)**(i-j)) / (m*((i-j)**2) * (mydelta_x**2))
    return t_matrix


def run_dvr(filename, m1, m2):
    m1 *= au
    m2 *= au
    m = (m1 * m2) / (m1 + m2)
    v_values, v_matrix = potential_energy(filename)
    hyd_shifts_corrected = np.insert(hyd_shifts, 0, hyd_shifts[0] - (hyd_shifts[1] - hyd_shifts[0]))  # see note above in get_geoms fwn -> start_array
    dvr_grid, delta_x = get_grid_and_dx(v_values, hyd_shifts_corrected)
    kin_matrix = kinetic_energy(dvr_grid, delta_x, m)
    myenergy, mywfns = np.linalg.eigh(v_matrix + kin_matrix)
    np.save(fname=param.dvr_run_numb, arr=mywfns)
    return myenergy, mywfns, dvr_grid, v_values



def plot_dvr_wfns(mydvr_grid, mywfns, v_vals):
    plt.plot(mydvr_grid, mywfns[:, 0]**2)
    # plt.xlim([1, 2])
    # plt.show()
    plt.plot(mydvr_grid, v_vals)
    # plt.show()
    return None


energy, wfns, dvr_grid, v_values = run_dvr("pots2_" + str(param.dvr_run_numb), 15.999, 1.00784)
run = plot_dvr_wfns(dvr_grid, wfns, v_values)
