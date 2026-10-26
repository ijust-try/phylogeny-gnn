"""No-graph ablation of the improved GNN: identical architecture, config, seeds
and splits as run_improve.py's winner, but with an empty edge set. SAGEConv then
has no neighbours to aggregate, so each isolate is processed from its own 157
features only (root-weight linear layer + MLP head) -- isolating what the
phylogenetic graph itself contributes. Writes results/gnn_no_graph_ablation.json.
"""

import json
from pathlib import Path

import numpy as np
import torch

from .data import load_gnn_data
from .run_improve import COMMON, CONFIGS, SEEDS, run

RESULTS_PATH = Path("results/gnn_no_graph_ablation.json")
CONFIG_NAME = "sage1_mlphead"


def summarize(runs):
    return {f"{s}/{m}": {"mean": float(np.mean([r[k][s][m] for r in runs])),
                         "std": float(np.std([r[k][s][m] for r in runs]))}
            for k in ("tuned_per_drug_threshold",) for s in ("val", "test")
            for m in ("macro_f1", "micro_f1", "macro_roc_auc", "macro_pr_auc")}


def main():
    data = load_gnn_data()
    cfg = CONFIGS[CONFIG_NAME]
    out = {"config": {**COMMON, **cfg}, "seeds": SEEDS}

    out["with_graph"] = [run(data, cfg, s) for s in SEEDS]
    no_graph = dict(data, edge_index=torch.empty((2, 0), dtype=torch.long), edge_weight=torch.empty(0))
    out["no_graph"] = [run(no_graph, cfg, s) for s in SEEDS]
    out["summary"] = {"with_graph": summarize(out["with_graph"]), "no_graph": summarize(out["no_graph"])}

    RESULTS_PATH.write_text(json.dumps(out, indent=2))
    for name in ("with_graph", "no_graph"):
        sm = out["summary"][name]
        print(name, "test tuned macro-F1 %.4f +/- %.4f" % (sm["test/macro_f1"]["mean"], sm["test/macro_f1"]["std"]),
              "| per-seed", [round(r["tuned_per_drug_threshold"]["test"]["macro_f1"], 4) for r in out[name]])


if __name__ == "__main__":
    main()
