from addons import *

@ctest_labeler("quick;reks")
def test_reks_h4_fon_regimes():
    ctest_runner(__file__, extra_infiles=["h4_closed.json", "h4_opening.json", "h4_open.json"])
