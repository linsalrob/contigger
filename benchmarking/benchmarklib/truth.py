RELATIONSHIP_ORDER = {
    k: i
    for i, k in enumerate(
        [
            "EXACT_MATCH",
            "QUERY_CONTAINED_IN_TARGET",
            "TARGET_CONTAINED_IN_QUERY",
            "QUERY_SUFFIX_TO_TARGET_PREFIX",
            "TARGET_SUFFIX_TO_QUERY_PREFIX",
            "AMBIGUOUS_OVERLAP",
            "NO_RELATIONSHIP",
        ]
    )
}


def sort_relationships(rows):
    return sorted(
        rows,
        key=lambda r: (
            r["case_id"],
            r["query_id"],
            r["target_id"],
            RELATIONSHIP_ORDER.get(r["expected_relationship"], 99),
        ),
    )


def classify(identity, overlap, q_end_gap, t_start_gap, thresholds):
    ok_id = identity >= thresholds["identity_percent"]
    ok_len = overlap >= thresholds["minimum_overlap"]
    ok_end = q_end_gap <= thresholds["end_tolerance"] and t_start_gap <= thresholds["end_tolerance"]
    return ok_id and ok_len and ok_end
