import pytest
import pandas as pd
from pathlib import Path

FIXTURES_DIR = Path(__file__).parent / "fixtures"

@pytest.fixture
def simple_region_dict():
    return {(0, 10): [["A", "B"]]}

@pytest.fixture
def ibd_segments():
    return pd.read_feather(FIXTURES_DIR / "ibd_segments.feather")

@pytest.fixture
def example_map():
    return pd.read_csv(
        FIXTURES_DIR / "chr1.map",
        delim_whitespace=True,
        header=None,
        names=["CHROM", "ID", "CM", "POS"]
    )
