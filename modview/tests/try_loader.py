import pathlib

import numpy as np
import pytest


@pytest.fixture
def rootdir():
    return pathlib.Path(__file__).parent.resolve()
# rootdir is now defined as a fixture in this module. We can now use "rootdir" to find files within the current working directory. 

# this is an example of how a pytest function is made: (extracted from G. Voet's mixsea repository)
# Read CTD example data and make it available for the following tests
# We defined rootdir as a fixture in conftest.py
# and can use it here as input now
#@pytest.fixture
#def ctd_profile(rootdir):
#    ctdtestfile = rootdir / "data/ctd_test_data.csv"
#    ctd_matrix = np.loadtxt(ctdtestfile, comments="#", delimiter=",")
#    # read variable names from comment line and divide data into a dict
#    with open(ctdtestfile) as f:
#        vars = f.readline().strip()
#    if vars[:2] == "# ":
#        v = [x.strip() for x in vars[2:].split(",")]
#        ctd = {vi: ctd_matrix[:, i] for i, vi in enumerate(v)}
#    return ctd
