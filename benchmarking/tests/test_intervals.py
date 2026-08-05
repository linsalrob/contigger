from benchmarklib.intervals import extract, merge_suffix_prefix
def test_circular(): assert extract("ABCDEFGH",6,10,circular=True)=="GHAB"
def test_reverse(): assert extract("AACCGG",1,5,"-")=="CGGT"
def test_overlap(): assert merge_suffix_prefix("AAAACCCC","CCCCGGGG",4)=="AAAACCCCGGGG"
def test_containment(): assert "CCCG" in "AAACCCGGG"
