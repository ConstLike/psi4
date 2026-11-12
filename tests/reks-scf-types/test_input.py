from pathlib import Path

from addons import *

@ctest_labeler("quick;reks")
def test_reks_scf_types():
    ctest_runner(__file__, extra_infiles=[str(Path(__file__).parent.parent / "reks-h2-fon-regimes" / "h2_opening.json")])
