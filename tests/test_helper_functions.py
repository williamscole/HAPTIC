import pytest
from helper_functions import split_regions, merge_regions, segment_coverage

class TestSplitRegions:

    def test_basic_overlap(self):
        region_dict = {(0, 10): [["A", "B"]]}
        new_region = [5, 10, ["C"]]
        result = split_regions(region_dict, new_region)
        assert result == {(0, 5): [['A', 'B']], (5, 10): [['A', 'B'], ['C']]}

    def test_no_overlap(self):
        region_dict = {(0, 10): [["A", "B"]]}
        new_region = [15, 25, ["C"]]
        result = split_regions(region_dict, new_region)
        assert result == {(0, 10): [["A", "B"]], (15, 25): [["C"]]}

    def test_complete_overlap(self):
        region_dict = {(0, 10): [["A", "B"]]}
        new_region = [0, 10, ["C"]]
        result = split_regions(region_dict, new_region)
        assert result == {(0, 10): [["A", "B"], ["C"]]}

    def test_new_region_contained_within_existing(self):
        region_dict = {(0, 20): [["A"]]}
        new_region = [5, 15, ["B"]]
        result = split_regions(region_dict, new_region)
        assert result == {(0, 5): [["A"]], (5, 15): [["A"], ["B"]], (15, 20): [["A"]]}

    def test_new_region_contains_existing(self):
        region_dict = {(5, 15): [["A"]]}
        new_region = [0, 20, ["B"]]
        result = split_regions(region_dict, new_region)
        assert result == {(0, 5): [["B"]], (5, 15): [["A"], ["B"]], (15, 20): [["B"]]}

    def test_partial_overlap_left(self):
        region_dict = {(5, 15): [["A"]]}
        new_region = [0, 10, ["B"]]
        result = split_regions(region_dict, new_region)
        assert result == {(0, 5): [["B"]], (5, 10): [["A"], ["B"]], (10, 15): [["A"]]}

    def test_partial_overlap_right(self):
        region_dict = {(0, 10): [["A"]]}
        new_region = [5, 15, ["B"]]
        result = split_regions(region_dict, new_region)
        assert result == {(0, 5): [["A"]], (5, 10): [["A"], ["B"]], (10, 15): [["B"]]}

    def test_multiple_regions_some_overlapping(self):
        region_dict = {(0, 10): [["A"]], (20, 30): [["B"]], (40, 50): [["C"]]}
        new_region = [5, 25, ["D"]]
        result = split_regions(region_dict, new_region)
        assert result == {
            (0, 5): [["A"]],
            (5, 10): [["A"], ["D"]],
            (10, 20): [["D"]],
            (20, 25): [["B"], ["D"]],
            (25, 30): [["B"]],
            (40, 50): [["C"]]
        }

    def test_empty_region_dict(self):
        region_dict = {}
        new_region = [0, 10, ["A"]]
        result = split_regions(region_dict, new_region)
        assert result == {(0, 10): [["A"]]}

class TestMergeRegions:

    def test_basic_merge(self):
        regions = [(0, 20, "A")]
        new_region = [21, 30, "A"]
        result = merge_regions(regions, new_region, gap=1)
        assert result == [(0, 30, "A")]

    def test_no_merge_different_info(self):
        regions = [(0, 20, "A")]
        new_region = [21, 30, "B"]
        result = merge_regions(regions, new_region, gap=1)
        assert (0, 20, "A") in result
        assert (21, 30, "B") in result

    def test_no_merge_gap_too_large(self):
        regions = [(0, 20, "A")]
        new_region = [25, 30, "A"]
        result = merge_regions(regions, new_region, gap=1)
        assert (0, 20, "A") in result
        assert (25, 30, "A") in result

    def test_merge_overlapping(self):
        regions = [(0, 20, "A")]
        new_region = [10, 30, "A"]
        result = merge_regions(regions, new_region, gap=0)
        assert result == [(0, 30, "A")]

    def test_merge_adjacent_no_gap(self):
        regions = [(0, 20, "A")]
        new_region = [20, 30, "A"]
        result = merge_regions(regions, new_region, gap=0)
        assert result == [(0, 30, "A")]

    def test_empty_regions(self):
        regions = []
        new_region = [0, 10, "A"]
        result = merge_regions(regions, new_region, gap=0)
        assert result == [(0, 10, "A")]

    def test_multiple_regions_one_merges(self):
        regions = [(0, 10, "A"), (20, 30, "B")]
        new_region = [11, 15, "A"]
        result = merge_regions(regions, new_region, gap=1)
        assert (0, 15, "A") in result
        assert (20, 30, "B") in result

class TestSegmentCoverage:
    # gap parameter convention:
    # positive gap: merge regions within that many cM of each other (gap=1 merges regions 1 cM apart)
    # negative gap: regions must overlap by more than abs(gap) cM to merge (gap=-3 requires 3+ cM overlap)

    def test_overlapping_regions(self):
        result = segment_coverage([[0, 10, "A"], [5, 15, "B"]])
        assert result == [[0, 15, {"A", "B"}]]

    def test_non_overlapping_regions(self):
        result = segment_coverage([[0, 10, "A"], [20, 30, "B"]])
        assert [0, 10, {"A"}] in result
        assert [20, 30, {"B"}] in result

    def test_single_region(self):
        result = segment_coverage([[0, 10, "A"]])
        assert result == [[0, 10, {"A"}]]

    def test_three_regions_all_overlapping(self):
        result = segment_coverage([[0, 10, "A"], [5, 15, "B"], [8, 20, "C"]])
        assert result == [[0, 20, {"A", "B", "C"}]]

    def test_gap_parameter_merges_close_regions(self):
        # gap=1 merges regions that are at most 1 cM apart
        result = segment_coverage([[0, 10, "A"], [11, 20, "B"]], gap=2)
        assert result == [[0, 20, {"A", "B"}]]

    def test_gap_parameter_does_not_merge_far_regions(self):
        result = segment_coverage([[0, 10, "A"], [15, 20, "B"]], gap=1)
        assert [0, 10, {"A"}] in result
        assert [15, 20, {"B"}] in result

    def test_gap_negative_requires_minimum_overlap(self):
        # gap=-3 means regions must overlap by MORE THAN 3 cM to merge
        # (0,10) and (6,20) overlap by 4 cM, so they merge
        result = segment_coverage([[0, 10, "A"], [6, 20, "B"]], gap=-3)
        assert result == [[0, 20, {"A", "B"}]]

    def test_gap_negative_no_merge_insufficient_overlap(self):
        # (0,10) and (7,20) overlap by exactly 3 cM, which is not > 3, so they don't merge
        result = segment_coverage([[0, 10, "A"], [7, 20, "B"]], gap=-3)
        assert [0, 10, {"A"}] in result
        assert [7, 20, {"B"}] in result

    def test_same_info_overlapping(self):
        result = segment_coverage([[0, 10, "A"], [5, 15, "A"]])
        assert result == [[0, 15, {"A"}]]

    def test_empty_input(self):
        result = segment_coverage([])
        assert result == []