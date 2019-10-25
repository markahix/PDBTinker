"""
PDB-Tinker
Converts PDBs to Tinker XYZs
"""

# Add imports here
from .pdbtinker import *
import numpy as np
import sys
import os
import PDBTinker.dictionaries as dictionaries

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

dictpath = os.path.dirname(dictionaries.__file__)
amino_acid_list = np.load(dictpath+"/amino_acid_list.npy")
nucleic_acid_list = np.load(dictpath+"/nucleic_acid_list.npy")
water_list = np.load(dictpath+"/water_list.npy")
ion_list = np.load(dictpath+"/ion_list.npy")
non_standard_list = np.load(dictpath+"/non_standard_list.npy")
cap_list = np.load(dictpath+"/cap_list.npy")
param_dict={}
