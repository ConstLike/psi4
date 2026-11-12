from addons import *

@ctest_labeler("quick;reks")
def test_reks_c2h4_si_properties():
    ctest_runner(__file__, extra_infiles=["c2h4_opening.json", "c2h4_open.json"])
