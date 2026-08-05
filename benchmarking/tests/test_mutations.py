from benchmarklib.mutations import substitute


def test_deterministic():
    assert substitute("ACGT" * 3000, 199, 7) == substitute("ACGT" * 3000, 199, 7)


def test_identity_counts():
    for n in (199, 200, 201):
        out, rows = substitute("ACGT" * 2500, n, 47291 + n)
        assert len(rows) == n
        assert sum(a != b for a, b in zip(out, "ACGT" * 2500, strict=True)) == n
