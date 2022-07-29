import unittest
import pickle
import numpy as np
from orix.crystal_map import utilities


class TestUtilities(unittest.TestCase):
    def test_spatial_decomposition(self):
        # Open pickle
        my_pickle_jar = open(f"C:/git/Mart2Aust_Hackathon/orix/tests/crystal_map/X.pickle", "rb")
        self.X = pickle.load(my_pickle_jar)
        my_pickle_jar.close()
        my_pickle_jar = open("C:/git/Mart2Aust_Hackathon/orix/tests/crystal_map/unitcell.pickle", "rb")
        self.uc = pickle.load(my_pickle_jar)
        my_pickle_jar.close()
        my_pickle_jar = open(f"C:/git/Mart2Aust_Hackathon/orix/tests/crystal_map/V.pickle", "rb")
        self.V = pickle.load(my_pickle_jar)
        my_pickle_jar.close()
        my_pickle_jar = open(f"C:/git/Mart2Aust_Hackathon/orix/tests/crystal_map/F.pickle", "rb")
        self.F = pickle.load(my_pickle_jar)
        my_pickle_jar.close()
        my_pickle_jar = open(f"C:/git/Mart2Aust_Hackathon/orix/tests/crystal_map/I_FD.pickle", "rb")
        self.I_FD = pickle.load(my_pickle_jar)
        my_pickle_jar.close()

        self.actual_V, self.actual_F, self.actual_I_FD = utilities.spatial_decomposition(self.X, self.uc)

        # self.assertEqual(utilities.spatial_decomposition(self.X, self.uc), (self.V, self.F, self.I_FD))
        self.assertTrue(np.allclose(self.actual_V, self.V))
        self.assertTrue(np.allclose(self.actual_F, self.F))

        # Suite of tests for testing validity of output sparse matrix to due
        # inability to explicitly check equality due to shear size of matrix
        self.assertTrue((self.actual_I_FD != self.I_FD).nnz == 0)
        self.assertTrue(f"{type(self.actual_I_FD)}" == "<class 'scipy.sparse._csr.csr_matrix'>")
        self.assertTrue(np.shape(self.actual_I_FD) == np.shape(self.I_FD))


if __name__ == '__main__':
    unittest.main()
