import unittest

import numpy as np

from src.DataReader import DataReader


class DataReaderTest(unittest.TestCase):
    def test_mean_profiles(self):
        """

        :rtype: object
        """
        datareader = DataReader()

        incremental_array = np.arange(16).reshape(4, 4)
        incremental_array_expected_time_avg = [6, 7, 8, 9]
        incremental_array_actual_time_avg = datareader.__meanProfiles__(incremental_array, idx_t0=0, idx_t1=4).tolist()
        self.assertListEqual(incremental_array_expected_time_avg, incremental_array_actual_time_avg)

        one_through_four = np.array([[1, 1, 1, 1], [2, 2, 2, 2], [3, 3, 3, 3], [4, 4, 4, 4]])
        one_through_four_expected_time_avg = [2.5, 2.5, 2.5, 2.5]
        one_through_four_actual_time_avg = datareader.__meanProfiles__(one_through_four).tolist()
        self.assertListEqual(one_through_four_expected_time_avg, one_through_four_actual_time_avg)

        random_array = np.array([[72, 54, 13, 48],
                                 [943, 238, 2, 34],
                                 [245, 837, 3, 64],
                                 [42, 84, 992, 71]])

        rand_array_expected_time_avg = [325.5, 303.25, 252.5, 54.25]
        rand_array_actual_time_avg = datareader.__meanProfiles__(random_array).tolist()
        self.assertListEqual(rand_array_expected_time_avg, rand_array_actual_time_avg)

        rand_array_expected_time_avg_partial = [594., 537.5, 2.5, 49.]
        rand_array_actual_time_avg_partial = datareader.__meanProfiles__(random_array, idx_t0=1, idx_t1=3).tolist()
        self.assertListEqual(rand_array_expected_time_avg_partial, rand_array_actual_time_avg_partial)


if __name__ == '__main__':
    unittest.main()
