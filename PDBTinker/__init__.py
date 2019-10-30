"""
PDB-Tinker
Converts PDBs to Tinker XYZs
"""

# Add imports here
from .pdbtinker import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
