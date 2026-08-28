from kegg_string_mcp.evaluate.score import reported_pathways


def test_context_pathways_are_not_counted_as_claims():
    """A model annotating furA may look up katG as context. Counting katG's
    pathways as furA claims scored a correct annotation as fabrication."""
    summary = "furA has no KEGG pathways. Its neighbour katG is in mtu00360 and mtu01100."
    validation = {"citations": [{"identifier": "mtu00360", "status": "cross_target"},
                                {"identifier": "mtu01100", "status": "cross_target"}]}
    assert reported_pathways(summary, "mtu", validation) == set()
    # Without the validator's judgement the naive read counts them.
    assert reported_pathways(summary, "mtu") == {"mtu00360", "mtu01100"}

def test_own_pathways_are_still_counted():
    summary = "katG is in mtu00360."
    validation = {"citations": [{"identifier": "mtu00360", "status": "verified"}]}
    assert reported_pathways(summary, "mtu", validation) == {"mtu00360"}
