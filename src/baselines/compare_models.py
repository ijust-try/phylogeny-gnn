"""Print a side-by-side test-set comparison of GNN vs. all baselines from results/*.json."""

import json
from pathlib import Path

RESULTS = Path("results")


def load_models():
    models = {}
    for name, entry in json.load(open(RESULTS / "baseline_metrics.json"))["models"].items():
        models[name] = entry
    models["cnn1d"] = json.load(open(RESULTS / "cnn_metrics.json"))
    models["gnn (phylo GCN)"] = json.load(open(RESULTS / "gnn_metrics.json"))
    return models


def main():
    print(f"{'model':18s} {'F1@0.5':>8s} {'F1 tuned':>9s} {'micro-F1':>9s} {'ROC-AUC':>8s} {'PR-AUC':>7s}   (TEST split)")
    for name, m in load_models().items():
        d, t = m["default_threshold_0.5"]["test"], m["tuned_per_drug_threshold"]["test"]
        print(f"{name:18s} {d['macro_f1']:8.3f} {t['macro_f1']:9.3f} {t['micro_f1']:9.3f} "
              f"{t['macro_roc_auc']:8.3f} {t['macro_pr_auc']:7.3f}")


if __name__ == "__main__":
    main()
