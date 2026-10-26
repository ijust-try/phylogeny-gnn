"""Train baselines and the GNN (with / without graph) on VCF-derived SNP features
(data/processed/vcf_features.npz, built by src/phylogeny/vcf_features.py) instead of
the curated 157-mutation table. Same fixed splits, labels, graph, metrics and the
same GNN config as run_train.py's default (not re-tuned for the new features).
Writes results/vcf_feature_metrics.json.
"""

import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import torch

from ..baselines import train_mlp, train_random_forest, train_xgboost
from ..baselines.run_all import run_model
from .data import load_gnn_data
from .run_improve import CONFIGS, SEEDS, run

RESULTS_PATH = Path("results/vcf_feature_metrics.json")
GNN_CONFIG = "sage1_mlphead"
BASELINES = {"random_forest": train_random_forest.build_model, "xgboost": train_xgboost.build_model,
             "multi_task_mlp": train_mlp.build_model}


def test_f1(r):
    return r["tuned_per_drug_threshold"]["test"]["macro_f1"]


def main():
    data = load_gnn_data()
    X = sp.load_npz("data/processed/vcf_features.npz").toarray().astype(np.float32)
    if X.shape[0] != len(data["sample_ids"]):
        raise ValueError("vcf_features.npz row count != sample_ids")
    data["X"] = torch.tensor(X)
    m, y, cols = data["masks_np"], data["y"].numpy(), data["drug_columns"]
    print("VCF feature matrix:", X.shape, flush=True)

    out = {"n_features": int(X.shape[1]), "baselines": {}}
    for name, build in BASELINES.items():
        out["baselines"][name] = run_model(build, X[m["train"]], y[m["train"]], X[m["val"]], y[m["val"]],
                                           X[m["test"]], y[m["test"]], cols)
        print(name, "test tuned macro-F1 %.3f" % test_f1(out["baselines"][name]), flush=True)
        RESULTS_PATH.write_text(json.dumps(out, indent=2))

    cfg = CONFIGS[GNN_CONFIG]
    no_graph = dict(data, edge_index=torch.empty((2, 0), dtype=torch.long), edge_weight=torch.empty(0))
    for label, d in (("gnn_with_graph", data), ("gnn_no_graph", no_graph)):
        out[label] = [run(d, cfg, s) for s in SEEDS]
        v = [test_f1(r) for r in out[label]]
        print(label, "test tuned macro-F1 %.4f +/- %.4f" % (np.mean(v), np.std(v)), [round(x, 4) for x in v], flush=True)
        RESULTS_PATH.write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
