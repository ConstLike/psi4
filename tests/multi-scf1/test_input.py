from addons import *

@ctest_labeler("quick;scf;multi-scf")
def test_multi_scf1():
    ctest_runner(__file__)
