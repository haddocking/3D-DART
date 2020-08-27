import os
import shutil
import filecmp


class RegressionTest():
    def __init__(self):
        self.path = ''
        self.test_path = ''
        self.golden_data_path = ''

    def ini_test_path(self):
        try:
            if self.test_path:
                shutil.rmtree(self.test_path)
        except OSError:
            pass
        os.mkdir(self.test_path)

    def clean_test_path(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass


class TestRegressionNACustomBuild(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_regression_nacustombuild/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/NAcustombuild/'

    def teardown(self):
        self.clean_test_path()

    def test_workflow(self):
        os.chdir(self.test_path)
        
        assert False
