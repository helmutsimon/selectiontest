# coding=utf-8


import os
from cogent3.util.unit_test import TestCase, main
# from numpy.testing import assert_allclose
from cli import selectiontestcli
from selectiontest import selectiontest
from selectiontest.selectiontest import calculate_D
from pandas import read_csv
from numpy import isclose, sum, array
from collections import Counter
from click.testing import CliRunner


abspath = os.path.abspath(__file__)
projdir = "/".join(abspath.split("/")[:-1])
print('projdir=', projdir)

from cli import selectiontestcli

#sys.path.insert(0, projdir)





def compute_sfs(variant_array):
    n = variant_array.shape[1]
    occurrences = sum(variant_array, axis=1)
    sfs = Counter(occurrences)
    sfs = [sfs[i] for i in range(1, n)]
    return array(sfs)


print('selectiontest version: ', selectiontest.__version__)


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
        assert result.exit_code == 0, str(result.exit_code)
        assert isclose(-1.446172, tajd, atol=1e-5), 'Incorrect Taj D value (cli)' + str(tajd)


TestTajima().test_Tajima_module()
TestTajima().test_Tajima_cli()