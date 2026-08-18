from __future__ import annotations

import unittest
import numpy as np

from subset_states.tables import table_i_array


class TableIScalingTests(unittest.TestCase):
    def test_full_range_regressions_unchanged(self) -> None:
        table = table_i_array()
        n, M, S = table[:, 0], table[:, 1], table[:, 2]
        slope_M, intercept_M = np.polyfit(n, np.log2(M), 1)
        slope_S, intercept_S = np.polyfit(n, S, 1)
        self.assertAlmostEqual(slope_M, 0.7035409513877944, places=9)
        self.assertAlmostEqual(intercept_M, -0.35773436814817483, places=9)
        self.assertAlmostEqual(slope_S, 0.5093, places=9)
        self.assertAlmostEqual(intercept_S, -0.9900909091, places=9)

    def test_low_n_robustness_comparison(self) -> None:
        table = table_i_array()
        low = table[table[:, 0] <= 20]
        slope_M_low = np.polyfit(low[:, 0], np.log2(low[:, 1]), 1)[0]
        slope_S_low = np.polyfit(low[:, 0], low[:, 2], 1)[0]
        slope_M_full = np.polyfit(table[:, 0], np.log2(table[:, 1]), 1)[0]
        slope_S_full = np.polyfit(table[:, 0], table[:, 2], 1)[0]
        self.assertLess(abs(slope_M_low - slope_M_full) / abs(slope_M_full), 0.02)
        self.assertLess(abs(slope_S_low - slope_S_full) / abs(slope_S_full), 0.02)


if __name__ == "__main__":
    unittest.main()
