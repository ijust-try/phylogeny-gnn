"""Validation-only tuning of the GNN on VCF-derived SNP features, with and without
the phylogenetic graph (each mode gets its own best config -- no-graph tuned is
effectively a tuned MLP, so the comparison is fair).

Selection score (VAL only): mean val PR-AUC over drugs with >= MIN_VAL_POS val
positives -- threshold-free and not driven by ETH/LZD's 1-positive val noise.
Stage A picks the SNP budget (top-300 vs top-1000 per drug), stage B the
architecture; the winner of each mode is then rerun with 3 seeds and test is
read once. Writes results/vcf_gnn_tuned_metrics.json.
"""

import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import torch

from .data import load_gnn_data
from .run_improve import COMMON, SEEDS, run

RESULTS_PATH = Path("results/vcf_gnn_tuned_metrics.json")
MIN_VAL_POS = 20
FEATURE_SETS = {"k300": "vcf_features.npz", "k1000": "vcf_features_k1000.npz"}
BASE = dict(conv_type="sage", n_layers=1, head_hidden=128, lr=0.005)
ARCH_GRID = {
    "base": {},
    "hidden128": dict(hidden_dim=128),
    "dropout0.5": dict(dropout=0.5),
    "head256_dropout0.5": dict(head_hidden=256, dropout=0.5),
    "sage2": dict(n_layers=2),
    "lr0.002": dict(lr=0.002),
}


def val_score(r):
    pd_ = r["tuned_per_drug_threshold"]["val"]["per_drug"]
    return float(np.mean([v["pr_auc"] for v in pd_.values() if v["n_positive"] >= MIN_VAL_POS and v["pr_auc"] is not None]))


def test_f1(r):
    return r["tuned_per_drug_threshold"]["test"]["macro_f1"]


def with_features(data, path, graph):
    X = torch.tensor(sp.load_npz(f"data/processed/{path}").toarray().astype(np.float32))
    d = dict(data, X=X)
    if not graph:
        d.update(edge_index=torch.empty((2, 0), dtype=torch.long), edge_weight=torch.empty(0))
    return d


def main():
    data = load_gnn_data()
    out = {"selection_metric": f"val mean PR-AUC, drugs with >= {MIN_VAL_POS} val positives", "modes": {}}
    for mode, graph in (("with_graph", True), ("no_graph", False)):
        rec = out["modes"][mode] = {"stage_A": {}, "stage_B": {}}
        for fs, path in FEATURE_SETS.items():
            r = run(with_features(data, path, graph), {**BASE}, SEEDS[0])
            rec["stage_A"][fs] = {"val_score": val_score(r), "n_features": int(sp.load_npz(f"data/processed/{path}").shape[1])}
            print(mode, "A", fs, round(val_score(r), 4), flush=True)
        best_fs = max(rec["stage_A"], key=lambda k: rec["stage_A"][k]["val_score"])
        d = with_features(data, FEATURE_SETS[best_fs], graph)
        for name, delta in ARCH_GRID.items():
            r = run(d, {**BASE, **delta}, SEEDS[0])
            rec["stage_B"][name] = {"val_score": val_score(r), "config": {**BASE, **delta}}
            print(mode, "B", name, round(val_score(r), 4), flush=True)
            RESULTS_PATH.write_text(json.dumps(out, indent=2))
        best_arch = max(rec["stage_B"], key=lambda k: rec["stage_B"][k]["val_score"])
        runs = [run(d, rec["stage_B"][best_arch]["config"], s) for s in SEEDS]
        rec.update(best_feature_set=best_fs, best_arch=best_arch, final_runs=runs,
                   test_macro_f1_mean=float(np.mean([test_f1(r) for r in runs])),
                   test_macro_f1_std=float(np.std([test_f1(r) for r in runs])))
        print(mode, "WINNER", best_fs, best_arch, "test tuned macro-F1 %.4f +/- %.4f" % (
            rec["test_macro_f1_mean"], rec["test_macro_f1_std"]), flush=True)
        RESULTS_PATH.write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
