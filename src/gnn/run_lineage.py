"""Lineage classification (single-label, 9 classes) with the GNN, with and without
the phylogenetic graph, plus a Random Forest baseline, on two feature sets:
  curated -- features.csv's 157 resistance mutations (carry little lineage info)
  vcf     -- data/processed/vcf_features.npz (932 SNPs picked from train AMR labels)
Same fixed splits as everything else. Isolates with missing lineage (lineage_id=-1)
are masked out of the loss and all metrics. Loss: class-weighted cross-entropy
(inverse train frequency). Early stopping on val macro-F1. GNN config is the
AMR default (GraphSAGE 1 layer, MLP head, lr 0.005) -- NOT re-tuned for lineage.
Writes results/lineage_metrics.json.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
import torch.nn.functional as F
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, f1_score

from .data import load_gnn_data
from .model import MultiTaskGCN

RESULTS_PATH = Path("results/lineage_metrics.json")
SEEDS = [42, 43, 44]
MAX_EPOCHS, PATIENCE, LR = 300, 30, 0.005


def load_lineage(sample_ids):
    y_df = pd.read_csv("data/processed/y_lineage.csv")
    if y_df["Name"].tolist() != sample_ids:
        raise ValueError("y_lineage.csv order differs from sample_ids.csv")
    meta = json.load(open("data/processed/y_lineage_metadata.json"))
    names = sorted(meta["classes"], key=meta["classes"].get)
    return y_df["lineage_id"].to_numpy(), names


def evaluate(y_true, y_pred, names):
    labels = list(range(len(names)))
    per_class_f1 = f1_score(y_true, y_pred, labels=labels, average=None, zero_division=0)
    return {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "macro_f1": float(per_class_f1.mean()),
        "per_class": {n: {"f1": float(per_class_f1[i]), "support": int((y_true == i).sum())} for i, n in enumerate(names)},
    }


def train_gnn(X, edge_index, edge_weight, y, masks, n_classes, seed):
    torch.manual_seed(seed)
    model = MultiTaskGCN(X.shape[1], 64, n_classes, dropout=0.3, skip_connection=True,
                         conv_type="sage", n_layers=1, head_hidden=128)
    opt = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=5e-4)
    y_train = y[masks["train"]]
    counts = torch.bincount(y_train, minlength=n_classes).float().clamp(min=1)
    weight = counts.sum() / (n_classes * counts)
    best, best_state, stale = -1.0, None, 0
    for _ in range(MAX_EPOCHS):
        model.train()
        opt.zero_grad()
        F.cross_entropy(model(X, edge_index, edge_weight)[masks["train"]], y_train, weight=weight).backward()
        opt.step()
        model.eval()
        with torch.no_grad():
            pred = model(X, edge_index, edge_weight)[masks["val"]].argmax(1).numpy()
        val_f1 = f1_score(y[masks["val"]].numpy(), pred, average="macro", zero_division=0)
        if val_f1 > best:
            best, stale = val_f1, 0
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
        else:
            stale += 1
            if stale >= PATIENCE:
                break
    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        return model(X, edge_index, edge_weight).argmax(1).numpy()


def main():
    data = load_gnn_data()
    lineage_id, names = load_lineage(data["sample_ids"])
    valid = lineage_id >= 0
    masks = {k: torch.tensor(v & valid) for k, v in data["masks_np"].items()}
    y = torch.tensor(np.where(valid, lineage_id, 0), dtype=torch.long)
    y_np = y.numpy()
    print({k: int(v.sum()) for k, v in masks.items()}, flush=True)

    feature_sets = {
        "curated": data["X"],
        "vcf": torch.tensor(sp.load_npz("data/processed/vcf_features.npz").toarray().astype(np.float32)),
    }
    empty = (torch.empty((2, 0), dtype=torch.long), torch.empty(0))
    out = {"classes": names, "results": {}}
    for fs, X in feature_sets.items():
        te = masks["test"].numpy()
        rf = RandomForestClassifier(n_estimators=300, class_weight="balanced", random_state=42, n_jobs=-1)
        rf.fit(X[masks["train"]].numpy(), y_np[masks["train"].numpy()])
        out["results"][f"{fs}/random_forest"] = {"test": evaluate(y_np[te], rf.predict(X[masks["test"]].numpy()), names)}
        for mode, (ei, ew) in (("gnn_with_graph", (data["edge_index"], data["edge_weight"])), ("gnn_no_graph", empty)):
            runs = []
            for s in SEEDS:
                pred = train_gnn(X, ei, ew, y, masks, len(names), s)
                runs.append(evaluate(y_np[te], pred[te], names))
            out["results"][f"{fs}/{mode}"] = {"test_runs": runs}
        for k in [k for k in out["results"] if k.startswith(fs)]:
            r = out["results"][k]
            f1s = [x["macro_f1"] for x in r.get("test_runs", [r.get("test")])]
            accs = [x["accuracy"] for x in r.get("test_runs", [r.get("test")])]
            print(f"{k:28s} test macro-F1 {np.mean(f1s):.4f} +/- {np.std(f1s):.4f}  acc {np.mean(accs):.4f}", flush=True)
        RESULTS_PATH.write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
