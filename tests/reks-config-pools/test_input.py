from pathlib import Path

from addons import *

@ctest_labeler("quick;reks")
def test_reks_config_pools():
    ctest_runner(__file__, extra_infiles=[str(Path(__file__).parent.parent / "reks-h4-fon-regimes" / "h4_opening.json")])
