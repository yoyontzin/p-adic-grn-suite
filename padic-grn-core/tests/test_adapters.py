from padic_grn_core.adapters import from_discrete_system, from_grn_table


def test_from_grn_table_shape() -> None:
    payload = from_grn_table({0: 0, 1: 1, 2: 2, 3: 3}, n=2, p=2)
    assert payload["kind"] == "grn_table"
    assert payload["expected_states"] == 4
    assert payload["num_transitions"] == 4


def test_from_discrete_system_metadata() -> None:
    meta = {"source": "unit-test"}
    payload = from_discrete_system({0: 0, 1: 1, 2: 2, 3: 3}, n=2, p=2, metadata=meta)
    assert payload["kind"] == "discrete_system"
    assert payload["metadata"]["source"] == "unit-test"
