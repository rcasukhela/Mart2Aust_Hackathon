"""
Created on Fri May 13, 2022
@author: StevenEgnaczyk
"""

import math
import numpy
import numpy as np
import scipy.special as sc


class Psi:
    """
    Psi Object Class.
    """

    # ============================ #
    #          Constructor         #
    # ============================ #
    def __init__(self, kappa, c, a, bandwidth):
        """
        Generates a Psi object, given Kappa, C, A, and bandwidth values

        :param kappa: (Integer)
            # Note: Add what Kappa is

        :param c: (Array)
            # Note: Add what Kappa is

        :param a: (Array)
            # Note: Add what a is

        :param bandwidth: (Float Value)
            # Note: Add what bandwidth is
        """

        self.kappa = kappa
        self.C = c
        self.A = a
        self.bandwidth = bandwidth


class ODF:
    """
    ODF Object Class.
    """

    # ============================ #
    #          Constructor         #
    # ============================ #
    def __init__(self, f_hat, weights, bandwidth, CS, SS, odfKernel, center):
        """
        :param f_hat:
            # To-Do: Add what this is

        :param weights:
            # To-Do: Add what this is

        :param bandwidth:
            # To-Do: Add what this is

        :param CS:
            # To-Do: Add what this is

        :param SS:
            # To-Do: Add what this is

        :param odfKernel:
            # To-Do: Add what this is
        """

        # Set variables
        self.f_hat = f_hat
        self.weights = weights
        self.bandwidth = bandwidth
        self.CS = CS
        self.SS = SS
        self.odfKernel = odfKernel
        self.center = center


# ============================ #
#          Methods             #
# ============================ #

def generate_kernel(halfwidth):
    """
    Creates a Kernel Object given the half-width.

    :param halfwidth: (Float)
        Add what this is

    :return:
        A Psi Object, the Kernel
    """

    # Calculate kappa and the pre-cos constant
    kappa = 0.5 * math.log(0.5) / math.log(math.cos(halfwidth / 2))
    C = sc.beta(1.5, 0.5) / sc.beta(1.5, kappa + 0.5)
    # Extract the bandwidth as an int
    L = math.floor(kappa)

    # Compute Chebyshev coefficients
    A = numpy.ones((L + 1,))  # 0 to L
    A[1] = kappa / (kappa + 2)  # Overwrite second value
    for j in range(1, L):
        A[j + 1] = ((kappa - j + 1) * A[j - 1] - (2 * j + 1) * A[j]) \
                   / (kappa + j + 2)

    for i in range(0, L):
        A[i] = (2 * i + 1) * A[i]

    A = cut(A)

    # Set the bandwidth per length of A
    bandwidth = A.size - 1

    return Psi(kappa, C, A, bandwidth)


def cut(A):
    """
    Cuts the Chebyshev coefficients off when they get too small

    :param A:
        List of Chebyshev coefficients

    :return:
        A smaller list of Chebyshev coefficients
    """

    # Epsilon value denoting the cutoff point for Chebyshev values
    epsilon = 0.0100 / 150

    returnA = A / numpy.arange(1, A.size + 1) ** 2
    ind = numpy.argwhere(returnA[1:] <= numpy.maximum(numpy.amin(
        numpy.append(returnA[1:], 10 * epsilon)), epsilon))[0]

    returnA = A[:numpy.minimum(int(ind + 1), A.size - 1) + 1]  # Added a +1 here!
    # Return the smaller list of Chebyshev coefficients
    return returnA

def eval_kernel_odf(odf, g):

    w = g-odf.center
    v = odf.psi.C * np.pow(math.cos(w/2), odf.psi.kappa)

    return v

def eval_kernel_odf(odf, g):
    w = g - odf.center
    v = odf.psi.C * np.pow(math.cos(w / 2), odf.psi.kappa)

    return v


def main():
    odfKernel = generate_kernel(0.04360)
    print(odfKernel.A)
    return odfKernel


kernel = main()
