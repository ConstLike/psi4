from pathlib import Path

from addons import *

@uusing("dftd3")
@ctest_labeler("quick;reks;dft")
def test_reks_dft_d3():
    ctest_runner(__file__, extra_infiles=[str(Path(__file__).parent.parent / "reks-h2-fon-regimes" / "h2_opening.json")])
