"""Las demos web contra fuerza bruta: niveles n = 0..N-1, B&B de mu_E y pi* de las tres redes."""

from __future__ import annotations

import importlib.util
import itertools
import json
import random
import shutil
import subprocess
from pathlib import Path

import pytest

WEB = Path(__file__).resolve().parents[1]
_BF_PATH = WEB.parent / "padic-grn-core" / "tests" / "_bruteforce.py"
_spec = importlib.util.spec_from_file_location("_bruteforce_core", _BF_PATH)
bf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bf)

DEMOS = ["padic-grn-demo.html", "padic-balls-animation.html"]
PI_STAR = {
    "athaliana": (13, [12, 3, 4, 8, 5, 1, 7, 2, 10, 11, 9, 6, 0], [26776, 53880, 25840]),
    "celegans": (8, [1, 0, 2, 3, 6, 5, 4, 7], [956, 392, 700]),
    "droso": (6, [0, 5, 3, 4, 2, 1], [64, 156, 164]),
}


def _cases() -> list[dict]:
    out = []
    for N, p, seed, n_const in bf.CASES:
        table = bf.random_network(N, p, seed, n_const)
        perms = [list(A) for A in itertools.permutations(range(N))]
        orderings = perms if len(perms) <= 24 else random.Random(seed).sample(perms, 24)
        out.append({"N": N, "p": p, "tt": [table[x] for x in range(p**N)], "table": table,
                    "orderings": orderings, "bnb": p == 2})
    return out


def _run(demo: str, cases: list[dict], tmp_path: Path) -> dict:
    path = tmp_path / "cases.json"
    path.write_text(json.dumps([{k: v for k, v in c.items() if k != "table"} for c in cases]), encoding="utf-8")
    res = subprocess.run(["node", str(WEB / "tests" / "demo_harness.js"), str(WEB / demo), str(path)],
                         capture_output=True, text=True, check=True, timeout=120)
    return json.loads(res.stdout)


@pytest.mark.skipif(shutil.which("node") is None, reason="node no disponible")
@pytest.mark.parametrize("demo", DEMOS)
def test_demo_matches_bruteforce(tmp_path: Path, demo: str) -> None:
    cases = _cases()
    out = _run(demo, cases, tmp_path)
    for c, r in zip(cases, out["cases"]):
        N, p, table = c["N"], c["p"], c["table"]
        for k, A in enumerate(c["orderings"]):
            counts = bf.level_counts(table, N, p, A)
            mu = list(bf.scores(table, N, p, A))
            assert r["fast"][k] == mu, (demo, N, p, A)
            if r["scores"]:
                s = r["scores"][k]
                assert s["mu"] == mu
                assert sorted(int(n) for n in s["E"]) == list(range(N))
                for n in range(N):
                    assert (s["E"][str(n)], s["A"][str(n)], s["I"][str(n)]) == counts[n]
            if r["levels"]:
                lv = r["levels"][k]
                assert sorted(int(n) for n in lv) == list(range(N))
                for n in range(N):
                    types = lv[str(n)]
                    assert (types.count("E"), types.count("A"), types.count("I")) == counts[n]
        if c["bnb"]:
            best_E = min(bf.scores(table, N, p, list(A))[0] for A in itertools.permutations(range(N)))
            argmin = sorted(list(A) for A in itertools.permutations(range(N)) if bf.scores(table, N, p, list(A))[0] == best_E)
            assert r["bnb"]["best"] == best_E
            assert sorted(o["A"] for o in r["bnb"]["optima"]) == argmin
            for o in r["bnb"]["optima"]:
                assert o["mu"] == list(bf.scores(table, N, p, o["A"]))


@pytest.mark.skipif(shutil.which("node") is None, reason="node no disponible")
@pytest.mark.parametrize("demo", DEMOS)
def test_demo_pi_star(tmp_path: Path, demo: str) -> None:
    nets = _run(demo, [], tmp_path)["networks"]
    assert set(nets) == set(PI_STAR)
    for key, (N, A, triple) in PI_STAR.items():
        net = nets[key]
        assert (net["N"], net["p"], net["optimal"]) == (N, 2, A)
        assert net["ttLength"] == 2**N
        assert net["fast"] == triple
        if net["expected"] is not None:
            assert net["expected"] == triple
