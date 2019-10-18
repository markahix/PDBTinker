"""
Unit and regression test for the PDBTinker package.
"""

# Import package, test suite, and other packages as needed
import PDBTinker
import pytest
import sys

def test_PDBTinker_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "PDBTinker" in sys.modules
