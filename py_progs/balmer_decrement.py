#!/usr/bin/env python
"""
balmer_decrement.py

Runs tests of of the Balmer decrement for a one zone thin shell Python model
(found in examples/regress/matom_balmer).
Involves running py_wind on a wind_save file and reading some output files.
Compares to Osterbrock values.

Usage:
    python BalmerDecrement.py root_filename
    python BalmerDecrement.py -h for help

Requirements:
    py_wind
    numpy
    matplotlib
    py_plot_util, py_read_output from $PYTHON/py_progs in the python path
"""
import numpy as np
import py_plot_util as util
import sys
import py_read_output as rd

# some global parameters to test against
OSTERBROCK_RATIOS = np.array([2.86, 1.0, 0.470, 0.262, 0.159, 0.107, 0.0748, 0.0544])
OSTERBROCK_INTENSITY = 8.30e-26
TOLERANCE = 0.3


def balmer_test(file_root: str, plot_decrement: bool = True) -> bool:
    """
    runs tests of of the Balmer decrement for a one zone
    thin shell Python model. Involves running py_wind on
    a wind_save file and reading some output files. Compares
    to Osterbrock values.

    :param file_root: The filename to run on
    :param plot_decrement: Whether or not to create plots of the output
    :return: True if successful, False if not
    """

    print("Running Balmer Test for run {}...".format(file_root))

    # create the list of commands to run in py wind
    nlevels = 8
    cmds = ["1", "s", "n", "i", "1", "1", "1", "1", "2", "0", "M", "2"]
    for i in range(nlevels):
        cmds.append("{:d}".format(i + 3))
    cmds.append("-1")
    cmds.append("q")

    # run py wind
    util.run_py_wind(file_root, cmds=cmds)

    # these could be in principle be used to check absolute emissivity values
    # ne = rd.read_pywind("{}.ne.dat".format(root), mode="1d")[2][1]
    # vol = rd.read_pywind("{}.vol.dat".format(root), mode="1d")[2][1]
    # nh1 = rd.read_pywind("{}.ioncH1.dat".format(root), mode="1d")[2][1]
    # nh2 = rd.read_pywind("{}.ioncH2.dat".format(root), mode="1d")[2][1]
    # nprot = nh1 + nh2

    # read the emissivities
    ratios = np.zeros(nlevels)
    for i in range(nlevels):
        ratios[i] = rd.read_pywind("{}.lev{}_emiss.dat".format(file_root, i + 3), mode="1d")[2][1]

    if plot_decrement:
        # make a scatter plot of the ratios to Hβ (Balmer decrement)
        n = np.arange(nlevels) + 3
        plt.scatter(n, ratios / ratios[1], c="k")
        plt.scatter(n, OSTERBROCK_RATIOS, facecolors="None", edgecolors="k", s=80)
        plt.xlabel("$n_u$", fontsize=16)
        plt.ylabel("Balmer Decrement", fontsize=16)
        plt.savefig("BalmerDecrement_{}.png".format(file_root))

    # define a 0 or 1 pass or fail
    pass_fail = (np.fabs(ratios / ratios[1] - OSTERBROCK_RATIOS) / OSTERBROCK_RATIOS)
    print("\n----------------------------------")
    print("\n\nArray of line ratio relative errors for Balmer series:\n", pass_fail)

    return np.all(pass_fail < TOLERANCE)


if __name__ == "__main__":

    if len(sys.argv) > 1:
        if sys.argv[1] == "-h":
            print(__doc__)
        else:
            plot = True
            # only plot if matplotlib installed
            try:
                import matplotlib.pyplot as plt
            except ModuleNotFoundError:
                print("No matplotlib, so not making plot.")
                plot = False

            # Tell the user whether the test is passed.
            if not balmer_test(sys.argv[1], plot_decrement=plot):
                raise ValueError('Balmer emissivities did not match those expected')

            else:
                print('\nTest passed.')
    else:
        print(__doc__)
