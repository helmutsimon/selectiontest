# coding=utf-8


import os
from cogent3.util.unit_test import TestCase, main
from numpy.testing import assert_allclose
from cli import selectiontestcli
from selectiontest import selectiontest
from selectiontest.selectiontest import calculate_D, sample_wf_distribution, sample_uniform_distribution
from selectiontest.selectiontest import test_neutrality, generate_sfs_array, compute_threshold
from pandas import read_csv, options
from numpy import isclose, sum, array, all, random, isfinite
from collections import Counter
from click.testing import CliRunner


__author__ = "Helmut Simon"
__copyright__ = "Copyright 2020, Helmut Simon"
__license__ = "GPL"
__version__ = "0.0.4"
__email__ = "helmut.simon@anu.edu.au"
__status__ = "Test"


abspath = os.path.abspath(__file__)
projdir = "/".join(abspath.split("/")[:-1])


from cli import selectiontestcli


def compute_sfs(variant_array):
    n = variant_array.shape[1]
    occurrences = sum(variant_array, axis=1)
    sfs = Counter(occurrences)
    sfs = [sfs[i] for i in range(1, n)]
    return array(sfs)


class TestTajima(TestCase):
    """Test with data at https://ocw.mit.edu/courses/health-sciences-and-technology/\
    hst-508-quantitative-genomics-fall-2005/study-materials/tajimad1.pdf.
    Data table is berwick_tajima_data.csv
    Results are:
    10 sequences, 41 sites
    pi: 3.888889
    Segregating sites: 16/41
    theta_hat[estimated from S]: 5.655772
    Tajima’s D: -1.446172
    a1=2.828968 a2=1.539768 b1=0.407407 b2=0.279012
    c1=0.053922 c2=0.047227 e1=0.019061 e2=0.004949  """
    options.mode.chained_assignment = None

    def __init__(self):
        # Need to put the data in github
        path = '/Users/helmutsimon/Google Drive/Genetics/Bayes SFS/Neutrality test'
        berwick = read_csv(path + '/berwick_tajima_data.csv', header=None)
        variant_array = berwick.iloc[2:]
        variant_array.drop(variant_array.columns[[-1, ]], axis=1, inplace=True)
        variant_array = variant_array.reset_index(drop=True)
        variant_array = variant_array.replace(to_replace=['A', 'T'], value=[0, 1]).T
        self.sfs = compute_sfs(variant_array)

    def test_Tajima_module(self):
        sfs = self.sfs
        tajd = calculate_D(sfs)
        assert isclose(-1.446172, tajd, atol=1e-5), 'Incorrect Taj D value (module)' + str(tajd)

    def test_Tajima_cli(self):
        sfs = self.sfs
        args = [str(x) for x in sfs]
        args.insert(0, 'calculate-d')
        runner = CliRunner()
        result = runner.invoke(selectiontestcli, args=args)
        tajd = float(result.output)
        assert result.exit_code == 0, "Exit code = %s" % result.exit_code
        assert isclose(-1.446172, tajd, atol=1e-5), 'Incorrect Taj D value (cli)' + str(tajd)


class Test_wf_distribution(TestCase):

    def __init__(self):
        self.result = [[0.63994338, 0.23721024, 0.10814627, 0.0147001]
            , [0.30087787, 0.25618537, 0.22870322, 0.21423354]
            , [0.65233883, 0.19955567, 0.09403128, 0.05407422]
            , [0.67011101, 0.12997136, 0.10227359, 0.09764404]
            , [0.61546645, 0.26632306, 0.07652315, 0.04168734]
            , [0.45996546, 0.27783446, 0.17732139, 0.0848787]
            , [0.38108608, 0.22363021, 0.19958542, 0.19569829]
            , [0.39618915, 0.24045628, 0.19744864, 0.16590593]
            , [0.60662475, 0.18794874, 0.12206768, 0.08335883]
            , [0.37366753, 0.24777493, 0.19706648, 0.18149106]]

    def test_sample_wf_distribution(self):
        random.seed(7)
        x = sample_wf_distribution(5, 10)
        assert_allclose(x, self.result, err_msg="Failed test_sample_wf_distribution")


class Test_uniform_distribution(TestCase):

    def __init__(self):
        self.result = [[0.85246025, 0.09936954, 0.02485407, 0.02331613],
           [0.60634303, 0.19092616, 0.130652  , 0.07207881],
           [0.44714758, 0.27102789, 0.16660389, 0.11522064],
           [0.49712199, 0.23361777, 0.16482629, 0.10443395],
           [0.75227604, 0.20319687, 0.0371041 , 0.00742299],
           [0.48522513, 0.26172446, 0.1626422 , 0.09040821],
           [0.81234361, 0.15649158, 0.01866754, 0.01249728],
           [0.8397336 , 0.11600754, 0.03965138, 0.00460748],
           [0.73013953, 0.13813067, 0.07873048, 0.05299932],
           [0.46740108, 0.19488643, 0.17406816, 0.16364433]]

    def test_sample_uniform_distribution(self):
        random.seed(11)
        x = sample_uniform_distribution(5, 10)
        assert_allclose(x, self.result, err_msg="Failed test_sample_uniform_distribution", atol=1e-5)


class Test_test_neutrality(TestCase):
    def __init__(self):
        self.sfs = [1, 1, 3, 0, 7, 0]
        self.result = -1.1955828464862277

    def test_test_neutrality(self):
        random.seed(3)
        rho = test_neutrality(self.sfs)
        assert isclose(rho, self.result), "Failed test of test_neutrality." + str(rho) + str(self.result)

    def test_test_neutrality_cli(self):
        args = [str(x) for x in self.sfs]
        args.insert(0, 'test-neutrality')
        random.seed(3)
        runner = CliRunner()
        result = runner.invoke(selectiontestcli, args=args)
        rho = float(result.output)
        assert result.exit_code == 0, "Exit code = %s" % result.exit_code
        assert isclose(rho, self.result), 'Failed test_neutrality (cli)' + str(rho)


class Test_compute_threshold(TestCase):
    def __init__(self):
        self.n = 5
        self.seg_sites = 11
        self.threshold = 0.45336890179930120
        self.threshold_fpr = 0.203067434354744230
        self.threshold_cli = 0.47588892627212176
        self.fpr = 0.2
        self.sfs_array = [[4, 5, 1, 1],
                          [6, 2, 2, 1],
                          [4, 1, 4, 2],
                          [4, 5, 2, 0],
                          [5, 5, 0, 1],
                          [8, 3, 0, 0],
                          [6, 5, 0, 0],
                          [6, 3, 0, 2],
                          [6, 4, 0, 1],
                          [7, 2, 1, 1]]

    def test_generate_sfs_array(self):
        random.seed = 13
        x = generate_sfs_array(self.n, self.seg_sites, reps=10)
        rowsums = x.sum(axis=1)
        assert all(rowsums == self.seg_sites), "Failed test of compute_threshold (row sums of seg sites)."
        assert all(x == self.sfs_array), "Failed test of compute_threshold (generate sfs)."

    def test_compute_threshold(self):
        random.seed = 13
        y = compute_threshold(self.n, self.seg_sites, reps=10)
        assert isclose(y, self.threshold), "Failed test of compute_threshold." + str(y) + str(self.threshold)

    def test_compute_threshold_fpr(self):
        random.seed = 5
        y = compute_threshold(self.n, self.seg_sites, reps=10, fpr=self.fpr)
        assert isclose(y, self.threshold_fpr), "Failed test of compute_threshold (fpr=0.2)." + str(y) + \
                                               str(self.threshold_fpr)

    def test_compute_threshold_cli(self):
        args = [str(self.n), str(self.seg_sites)]
        args.insert(0, 'compute-threshold')
        args.append('-r')
        args.append(str(10))
        random.seed = 3
        runner = CliRunner()
        result = runner.invoke(selectiontestcli, args=args)
        assert result.exit_code == 0, "Exit code = %s" % result.exit_code
        threshold = float(result.output)
        assert isclose(threshold, self.threshold_cli), "Failed test of compute_threshold (cli)." + str(threshold) + \
                                               str(self.threshold_cli)


def main():
    TestTajima().test_Tajima_module()
    TestTajima().test_Tajima_cli()
    Test_wf_distribution().test_sample_wf_distribution()
    Test_uniform_distribution().test_sample_uniform_distribution()
    Test_test_neutrality().test_test_neutrality()
    Test_test_neutrality().test_test_neutrality_cli()
    Test_compute_threshold().test_generate_sfs_array()
    Test_compute_threshold().test_compute_threshold()
    Test_compute_threshold().test_compute_threshold_fpr()
    Test_compute_threshold().test_compute_threshold_cli()
    print("All tests for selectiontest version:", selectiontest.__version__,  "complete.")


if __name__ == '__main__':
    main()