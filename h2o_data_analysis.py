import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
# plt.rc('text', usetex=True)
# rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
# rc('text', usetex=True)
import numpy.linalg as linalg
import sys
import glob

hartree_conv = 219474.6
zpe_list = []
angst = 0.529177
run_walkers = ["1000", "3000", "5000"]
int_run_walkers = list(map(int, run_walkers))

""" energy analysis """
def zpe_convergence_plotting():
    plot_zpes = []
    plot_stds = []
    for walks in int_run_walkers:
        plt.rcParams.update({'font.size': 15})
        zpe_list = []
        for dmc_trial in range(1, 6):
            v_array_imported = np.load("h2o_analysis_data/" + str(walks) + "_" + str(dmc_trial) + "_varray.npy")
            v_array_imported = v_array_imported[1000:]  # allow for convergence
            zpe = (np.average(v_array_imported)) * hartree_conv
            zpe_list.append(zpe)

        run_std = np.std(zpe_list)
        zpe_list = np.average(zpe_list)
        plot_zpes.append(zpe_list)
        plot_stds.append(run_std)
    plt.errorbar(int_run_walkers, plot_zpes, yerr=plot_stds, solid_capstyle='projecting', capsize=5)
    plt.plot([1000, 5000], [4638.39, 4638.39], 'r--')
    plt.ylim([4620, 4650])
    plt.xlim([925, 5075])
    plt.title('ZPE convergence with walker variation')
    plt.xlabel("$N_{w}$")
    plt.ylabel('$cm^{-1}$')
    plt.grid()
    plt.savefig(fname="zpe_convergence", dpi=500, bbox_inches="tight")
    return None

zpe_convergence_plotting()

""" plot Vref convergence """

# plt.plot(np.array(Vref_array) * hartree_conv)
# plt.title("$V_{ref}$ convergence")
# plt.xlabel("imaginary time")
# plt.ylabel("$V_{ref}$")
# plt.grid()
# plt.show()
# -----------------------------------


def vref_convergence_plotting():
    plt.figure(figsize=(8, 5))
    for walks in int_run_walkers:
        plt.rcParams.update({'font.size': 15})
        # plt.figure(figsize=(8, 5))
        for dmc_trial in range(1, 2):
            v_array_imported = np.load("h2o_analysis_data/" + str(walks) + "_" + str(dmc_trial) + "_varray.npy")
            # plt.figure(figsize=(8, 5))
            plt.plot(np.array(v_array_imported) * hartree_conv, label=str(walks) + " walkers")
            plt.legend()
            plt.title('$V_{ref}$ convergence with walker variation')
            plt.xlabel("$\\tau$")
            plt.ylabel('$V_{ref} (cm^{-1})$')
            plt.grid()

    plt.savefig(fname="vref_convergence", dpi=500, bbox_inches="tight")

    return None

# vref_convergence_plotting()



""" get zpe """

# Vref_array = Vref_array[500:]  # allow for convergence
# zpe = np.average(Vref_array)
# print(Vref_array)
# print("zpe = " + str(zpe * hartree_conv))




"""walker test"""


def walker_test(high_file_num):
    for i in range(1, high_file_num + 1):
        wts_check = np.load("h2o_analysis_data/1000_1_weights" + str(i) + ".npy")
        wfn_check = np.load("h2o_analysis_data/1000_1_wfn" + str(i) + ".npy")
        print(wts_check.shape)
        print(wfn_check.shape)

"""
- combine 10 harvest wfns per trial into one overlay plot, x2 for r_oh1 and r_oh2
- results in 2 overlayplots plots x5 trials x3 walker amounts
"""


def overlayplots(walk_numbs_list):
    for num_walkers in walk_numbs_list:
        for dmc_trial in range(1, 6):
            psi_sqrd1s = []
            psi_sqrd2s = []
            for harvest_n in range(1, 11):

                wghts = np.load("h2o_analysis_data/" + str(num_walkers) + "_" + str(dmc_trial) + "_weights" + str(harvest_n) + ".npy") / angst
                coords = np.load("h2o_analysis_data/" + str(num_walkers) + "_" + str(dmc_trial) + "_wfn" + str(harvest_n) + ".npy") / angst  # multiplied twice in dmc algo

                """ project psi^2 onto OH bond """
                r_oh1 = linalg.norm(coords[:, 1] - coords[:, 2], axis=1)
                # print(r_oh1.shape)
                r_oh2 = linalg.norm(coords[:, 0] - coords[:, 2], axis=1)
                psi_sqrd1, bins1 = np.histogram(r_oh1, bins=25, range=(.5, 1.5), weights=wghts, density=True)  # bins not used
                psi_sqrd2, bins2 = np.histogram(r_oh2, bins=25, range=(.5, 1.5), weights=wghts, density=True)
                # print(r_oh1.shape)
                # print(wghts.shape)

                bin_centers1 = 0.5 * (bins1[:-1] + bins1[1:])
                bin_centers2 = 0.5 * (bins2[:-1] + bins2[1:])

                psi_sqrd1s.append(psi_sqrd1)
                psi_sqrd2s.append(psi_sqrd2)

            # -------------------------------------------------------------
            #     """ expectation value of position <x> along r_OH """
            #     expectation_pos1 = np.average(r_oh1, weights=wghts)
            #     dev_1 = (expectation_pos1)
            #     expectation_pos2 = np.average(r_oh2, weights=wghts)
            #     dev_2 = np.std(r_oh2)
            # -------------------------------------------------------------
            plt.rcParams.update({'font.size': 16})
            for j in psi_sqrd1s:
                plt.plot(bin_centers1, j)
            plt.title("$\\psi^2$ projection onto $r_{OH_1}$")
            plt.xlabel("$r_{OH_1} (\\AA)$")
            plt.ylabel("$\\psi^2$")
            plt.grid()
            plt.savefig(fname=str(num_walkers) + "_run" + str(dmc_trial) + "_wfn1plot", dpi=500, bbox_inches="tight")
            plt.close()

            plt.rcParams.update({'font.size': 16})
            for k in psi_sqrd2s:
                plt.plot(bin_centers2, k)
            plt.title("$\\psi^2$ projection onto $r_{OH_2}$")
            plt.xlabel("$r_{OH_2} (\\AA)$")
            plt.ylabel("$\\psi^2$")
            plt.grid()
            plt.savefig(fname=str(num_walkers) + "_run" + str(dmc_trial) + "_wfn2plot", dpi=500, bbox_inches="tight")
            plt.close()


# overlayplots(int_run_walkers)
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

"""
combine overlay plots into one plot for each walker amount (1000, 3000, 5000) with 5 psi^2 on each, only using r_OH1
"""


def combine_wfns_to_one(walk_numbs_list):
    for num_walkers in walk_numbs_list:
        for dmc_trial in range(1, 6):
            w = []
            c = []
            for harvest_n in range(1, 11):
                loaded_wts = np.load(
                    "h2o_analysis_data/" + str(num_walkers) + "_" + str(dmc_trial) + "_weights" + str(harvest_n) + ".npy"
                )
                w.extend(loaded_wts)
                loaded_cds = np.load(
                    "h2o_analysis_data/" + str(num_walkers) + "_" + str(dmc_trial) + "_wfn" + str(harvest_n) + ".npy"
                )
                c.extend(loaded_cds)

            w1_combo = np.array(w) / angst
            c1_combo = np.array(c) / angst

            """ project psi^2 onto OH1 bond """
            plt.rcParams.update({'font.size': 16})
            r_oh1z = linalg.norm(c1_combo[:, 1] - c1_combo[:, 2], axis=1)
            psi_sqrd1_combo, bins_hist_combo = np.histogram(r_oh1z, bins=25, range=(.5, 1.5), weights=w1_combo, density=True)
            bin_centers_combo = 0.5 * (bins_hist_combo[:-1] + bins_hist_combo[1:])
            plt.plot(bin_centers_combo, psi_sqrd1_combo)
            # plt.title("$\combining wfns for psi^2 at 1000 walkers$")
            plt.title("$\\psi^2$ combined for 1000 walker simulations")
            plt.xlabel("$r_{OH_1} (\\AA)$")
            plt.ylabel("$\\psi^2$")
            plt.grid()
        plt.savefig(fname=str(num_walkers) + "_5x_wfn1plot1000", dpi=500, bbox_inches="tight")
        plt.close()
    return None


# combine_wfns_to_one(int_run_walkers)
