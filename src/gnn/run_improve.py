"""Improvement sweep for the multi-task GNN, with validation-only model selection.

Every config is trained on train only (early stopping on val). The winning
config is chosen by VAL macro-F1 (per-drug thresholds tuned on val); test is
never used for any choice. The winner is then re-run with extra seeds so the
reported test number is a mean +/- std, not one lucky run. Same metrics code
and fixed splits as run_train.py / baselines. Writes results/gnn_improved_metrics.json.
"""

import json
import time
from pathlib import Path

import numpy as np
import torch

from ..baselines.metrics import evaluate_multilabel, tune_per_drug_thresholds
from .data import load_gnn_data
from .train import train_model

RESULTS_PATH = Path("results/gnn_improved_metrics.json")
SEEDS = [42, 43, 44]
COMMON = dict(hidden_dim=64, weight_decay=5e-4, dropout=0.3, max_epochs=300, patience=30,
              skip_connection=True, stop_metric="val_macro_pr_auc")
CONFIGS = {
    "gcn2_linear_lr0.01_f1stop (current)": dict(conv_type="gcn", n_layers=2, head_hidden=0, lr=0.01,
                                                stop_metric="val_macro_f1", patience=20),
    "gcn2_linear_prstop": dict(conv_type="gcn", n_layers=2, head_hidden=0, lr=0.005),
    "gcn1_mlphead": dict(conv_type="gcn", n_layers=1, head_hidden=128, lr=0.005),
    "gcn2_mlphead": dict(conv_type="gcn", n_layers=2, head_hidden=128, lr=0.005),
    "sage1_mlphead": dict(conv_type="sage", n_layers=1, head_hidden=128, lr=0.005),
    "sage2_mlphead": dict(conv_type="sage", n_layers=2, head_hidden=128, lr=0.005),
}


def run(data, cfg, seed):
    kwargs = {**COMMON, **cfg}
    t0 = time.time()
    model, info = train_model(data, seed=seed, use_edge_weight=True, verbose=False, **kwargs)
    with torch.no_grad():
        proba = torch.sigmoid(model(data["X"], data["edge_index"], data["edge_weight"])).numpy()
    m, cols, y = data["masks_np"], data["drug_columns"], data["y"].numpy()
    yv, yt, pv, pt = y[m["val"]], y[m["test"]], proba[m["val"]], proba[m["test"]]
    thr = tune_per_drug_thresholds(yv, pv, cols)
    return {
        "seed": seed, "best_epoch": info["best_epoch"], "fit_seconds": round(time.time() - t0, 1),
        "default_threshold_0.5": {"val": evaluate_multilabel(yv, pv, cols, 0.5),
                                  "test": evaluate_multilabel(yt, pt, cols, 0.5)},
        "tuned_per_drug_threshold": {"thresholds": thr,
                                     "val": evaluate_multilabel(yv, pv, cols, thr),
                                     "test": evaluate_multilabel(yt, pt, cols, thr)},
    }


def val_score(r):
    return r["tuned_per_drug_threshold"]["val"]["macro_f1"]


def main():
    data = load_gnn_data()
    sweep = {}
    for name, cfg in CONFIGS.items():
        r = run(data, cfg, SEEDS[0])
        sweep[name] = r
        print(f"{name:40s} val tuned F1={val_score(r):.4f}  epoch={r['best_epoch']}  {r['fit_seconds']}s", flush=True)
        RESULTS_PATH.write_text(json.dumps({"sweep": sweep}, indent=2))

    winner = max(sweep, key=lambda n: val_score(sweep[n]))  # VAL only
    print(f"Winner by val: {winner}", flush=True)
    final_runs = [sweep[winner]] + [run(data, CONFIGS[winner], s) for s in SEEDS[1:]]
    summary = {}
    for setting in ("default_threshold_0.5", "tuned_per_drug_threshold"):
        for split in ("val", "test"):
            for metric in ("macro_f1", "micro_f1", "macro_roc_auc", "macro_pr_auc"):
                v = [r[setting][split][metric] for r in final_runs]
                summary[f"{setting}/{split}/{metric}"] = {"mean": float(np.mean(v)), "std": float(np.std(v)), "values": v}
    out = {"seeds": SEEDS, "winner_by_val": winner, "winner_config": {**COMMON, **CONFIGS[winner]},
           "sweep": sweep, "final_runs": final_runs, "summary": summary}
    RESULTS_PATH.write_text(json.dumps(out, indent=2))
    for k in ("tuned_per_drug_threshold/test/macro_f1", "default_threshold_0.5/test/macro_f1"):
        print(k, round(summary[k]["mean"], 4), "+/-", round(summary[k]["std"], 4), summary[k]["values"])


if __name__ == "__main__":
    main()
