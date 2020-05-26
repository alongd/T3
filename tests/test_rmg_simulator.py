"""
t3 tests test_rmg_simulator module
"""

import os
import shutil
import unittest

from tandem.simulate.rmg_simulator import RMGSimulator

################################################################################

class RMGSimulatorTest(unittest.TestCase):
    """
    Contains unit tests for the RMGSimulator module, used for performing mechanism analysis with RMG.
    """
    @classmethod
    def setUpClass(cls):
        """
        A method that is run before all unit tests in this class.
        """
        # use data from RMG's superminimal example
        cls.data_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data', 'rmg_simulator_test')
        cls.observable_list = [{'label': 'H2O2', 'smiles': 'OO'}, {'label': 'H', 'smiles': '[H]'}]
        cls.run_directory = os.path.join(cls.data_dir, 'iteration_0')
        cls.analysis_directory = os.path.join(cls.run_directory, 'rmg_analysis')
        cls.atol = 1e-6
        cls.rtol = 1e-4
        cls.sa_threshold = 0.001
        cls.sa_atol = 1e-6
        cls.sa_rtol = 1e-4
        cls.verbose = False

    def test_set_up_no_sa(self):
        """
        Runs RMG's superminimal example without SA by testing the `set_up` method in the RMGSimulator class.
        """
        RMGSimulator_adapter = RMGSimulator(run_directory=self.run_directory,
                                            atol=self.atol,
                                            rtol=self.rtol,
                                            verbose=self.verbose,
                                            )
        # check that RMG produced the concentration profiles
        concentration_profiles = os.path.isfile(os.path.join(self.analysis_directory, 'solver', 'simulation_1_12.csv'))
        self.assertTrue(concentration_profiles)

    def test_get_sa_coefficients(self):
        """
        Runs RMG's superminimal example with SA by testing the `set_up` method in the RMGSimulator class.
        Then runs the `get_sa_coefficients()` method to test that RMGSimulator correctly parses the SA csv files
        to obtain sa_dict.
        """
        RMGSimulator_adapter = RMGSimulator(observable_list=self.observable_list,
                                            run_directory=self.run_directory,
                                            atol=self.atol,
                                            rtol=self.rtol,
                                            sa_threshold=self.sa_threshold,
                                            sa_atol=self.sa_atol,
                                            sa_rtol=self.sa_rtol,
                                            verbose=self.verbose,
                                            )
        # check that RMG did run SA
        performed_sa = os.path.isdir(os.path.join(self.analysis_directory, 'sa'))
        self.assertTrue(performed_sa)

        # check that RMG produced the corresponding csv files
        files = ['sensitivity_1_SPC_3.csv', 'sensitivity_1_SPC_7.csv', 'simulation_1_12.csv']
        for file in files:
            exist = os.path.isfile(os.path.join(self.analysis_directory, 'sa', 'solver', file))
            self.assertTrue(exit)

        # check that RMG can correctly parse the csv files to product the SA dictionary
        sa_dict = RMGSimulator_adapter.get_sa_coefficients()
        # check that we have over 100 time steps as an arbitrary number to indicate that the solver completed
        self.assertGreater(len(sa_dict['time']), 100)

    @classmethod
    def tearDownClass(cls):
        """
        A method that is run after all unit tests in this class.
        Delete all project directories created during these unit tests
        """
        if os.path.isdir(cls.analysis_directory):
            shutil.rmtree(cls.analysis_directory)


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
