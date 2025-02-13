import unittest
from ler.rates import LeR
from astropy.cosmology import LambdaCDM

class TestLensingSimulation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Set up resources before any tests are run."""
        cls.ler = LeR(
            npool=4,
            z_min=0.0,
            z_max=10.0,
            event_type='BBH',
            size=100000,
            batch_size=50000,
            cosmology=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7),
            json_file_names=dict(
                ler_params="ler_params.json",
                unlensed_param="unlensed_param.json",
                unlensed_param_detectable="unlensed_param_detectable.json",
                lensed_param="lensed_param.json",
                lensed_param_detectable="lensed_param_detectable.json"
            ),
            interpolator_directory='./interpolator_pickle',
            ler_directory='./ler_data',
            verbose=True
        )
    
    def test_unlensed_event_simulation(self):
        """Test simulation of unlensed events."""
        unlensed_param = self.ler.unlensed_cbc_statistics(size=100000)
        self.assertIsNotNone(unlensed_param)
        self.assertTrue(len(unlensed_param) > 0, "Unlensed events simulation failed to generate data.")

    def test_lensed_event_simulation(self):
        """Test simulation of lensed events."""
        lensed_param = self.ler.lensed_cbc_statistics(size=100000)
        self.assertIsNotNone(lensed_param)
        self.assertTrue(len(lensed_param) > 0, "Lensed events simulation failed to generate data.")

    def test_rate_comparison(self):
        """Test rate comparison between lensed and unlensed events."""
        rate_ratio, unlensed_param_detectable, lensed_param_detectable = self.ler.rate_comparison_with_rate_calculation(
            unlensed_param=None,
            snr_threshold_unlensed=8.0,
            lensed_param=None,
            snr_threshold_lensed=[8.0, 8.0],
            num_img=[1, 1],
            nan_to_num=True,
            detectability_condition="step_function"
        )
        self.assertIsNotNone(rate_ratio, "Failed to calculate rate ratio.")
        self.assertGreater(rate_ratio, 0, "Rate ratio is not greater than 0.")

    @classmethod
    def tearDownClass(cls):
        """Clean up resources after all tests are done."""
        # Clean-up code here (if necessary)

if __name__ == '__main__':
    unittest.main()
