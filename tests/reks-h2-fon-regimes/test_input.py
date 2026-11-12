from addons import *

@ctest_labeler("quick;smoke;reks")
def test_reks_h2_fon_regimes():
    ctest_runner(__file__, extra_infiles=["h2_closed.json", "h2_opening.json", "h2_open.json"])
