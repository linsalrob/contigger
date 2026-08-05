from benchmarklib.reads import stable_hash
from benchmarklib.truth import sort_relationships
from benchmarklib.utilities import gzip_bytes


def test_hash():
    assert stable_hash("read1") == stable_hash("read1")


def test_gzip():
    assert gzip_bytes(b"abc") == gzip_bytes(b"abc")


def test_truth_sort():
    r = [
        {
            "case_id": "b",
            "query_id": "q",
            "target_id": "t",
            "expected_relationship": "NO_RELATIONSHIP",
        },
        {"case_id": "a", "query_id": "q", "target_id": "t", "expected_relationship": "EXACT_MATCH"},
    ]
    assert [x["case_id"] for x in sort_relationships(r)] == ["a", "b"]
