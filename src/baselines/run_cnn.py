"""1D-CNN multi-task baseline for AMR prediction on the real Afro-TB features.

Each isolate's 157 binary mutation indicators (features.csv column order,
which groups mutations by gene) are treated as a 1-channel sequence. Stacked
Conv1d layers + global pooling feed one linear layer with 9 sigmoid logits
(one per drug). Caveat: the order is gene-grouped, not genomic-position
ordered, so convolutional locality is only a weak inductive bias here.

Same fixed splits, same BCE pos_weight scheme as the GNN, and the same
val-only threshold tuning as src/baselines/run_all.py. Writes
results/cnn_metrics.json.
"""

import json
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import f1_score

from .data_loader import load_aligned_data, split_masks
from .metrics import evaluate_multilabel, tune_per_drug_thresholds

RESULTS_PATH = Path("results/cnn_metrics.json")
SEED = 42
LR = 1e-3
BATCH_SIZE = 256
MAX_EPOCHS = 100
PATIENCE = 15


class MultiTaskCNN1D(nn.Module):
    def __init__(self, n_drugs, channels=(32, 64, 64), dropout=0.3):
        super().__init__()
        c1, c2, c3 = channels
        self.conv = nn.Sequential(
            nn.Conv1d(1, c1, kernel_size=5, padding=2), nn.ReLU(),
            nn.Conv1d(c1, c2, kernel_size=5, padding=2), nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Conv1d(c2, c3, kernel_size=3, padding=1), nn.ReLU(),
        )
        self.head = nn.Sequential(nn.Dropout(dropout), nn.Linear(2 * c3, n_drugs))

    def forward(self, x):
        h = self.conv(x.unsqueeze(1))  # (B, C, L)
        h = torch.cat([h.max(dim=2).values, h.mean(dim=2)], dim=1)
        return self.head(h)  # logits


def predict_proba(model, X):
    model.eval()
    with torch.no_grad():
        return torch.sigmoid(model(X)).numpy()


def main():
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    data = load_aligned_data()
    masks = split_masks(data)
    drug_columns = data["drug_columns"]
    X = torch.tensor(data["X"])
    y = torch.tensor(data["y"])
    X_train, y_train = X[masks["train"]], y[masks["train"]]
    X_val, y_val = X[masks["val"]], y[masks["val"]]
    X_test, y_test = X[masks["test"]], y[masks["test"]]

    n_pos = y_train.sum(dim=0)
    n_neg = y_train.shape[0] - n_pos
    pos_weight = torch.where(n_pos > 0, n_neg / n_pos.clamp(min=1), torch.ones_like(n_pos))
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)

    model = MultiTaskCNN1D(n_drugs=len(drug_columns))
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)

    best_f1, best_state, best_epoch, stale = -1.0, None, -1, 0
    t0 = time.time()
    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        perm = torch.randperm(X_train.shape[0])
        for start in range(0, len(perm), BATCH_SIZE):
            idx = perm[start:start + BATCH_SIZE]
            optimizer.zero_grad()
            criterion(model(X_train[idx]), y_train[idx]).backward()
            optimizer.step()
        val_pred = (predict_proba(model, X_val) >= 0.5).astype(int)
        val_f1 = f1_score(y_val.numpy().astype(int), val_pred, average="macro", zero_division=0)
        if val_f1 > best_f1:
            best_f1, best_epoch, stale = val_f1, epoch, 0
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
        else:
            stale += 1
            if stale >= PATIENCE:
                break
        if epoch % 10 == 0:
            print(f"  epoch {epoch:3d}  val_macro_f1={val_f1:.4f}")
    fit_seconds = time.time() - t0
    model.load_state_dict(best_state)
    print(f"Best epoch {best_epoch} (val macro-F1={best_f1:.4f}), {fit_seconds:.1f}s")

    val_proba, test_proba = predict_proba(model, X_val), predict_proba(model, X_test)
    y_val_np, y_test_np = y_val.numpy(), y_test.numpy()
    tuned = tune_per_drug_thresholds(y_val_np, val_proba, drug_columns)

    results = {
        "seed": SEED,
        "model": "multi_task_cnn1d",
        "n_train": int(masks["train"].sum()),
        "n_val": int(masks["val"].sum()),
        "n_test": int(masks["test"].sum()),
        "n_features": int(X.shape[1]),
        "drug_columns": drug_columns,
        "threshold_for_f1": 0.5,
        "hyperparameters": {"lr": LR, "batch_size": BATCH_SIZE, "max_epochs": MAX_EPOCHS,
                            "patience": PATIENCE, "n_parameters": sum(p.numel() for p in model.parameters())},
        "fit_seconds": round(fit_seconds, 2),
        "best_epoch": best_epoch,
        "default_threshold_0.5": {
            "val": evaluate_multilabel(y_val_np, val_proba, drug_columns, threshold=0.5),
            "test": evaluate_multilabel(y_test_np, test_proba, drug_columns, threshold=0.5),
        },
        "tuned_per_drug_threshold": {
            "thresholds": tuned,
            "tuned_on": "val (F1-maximizing per drug via precision_recall_curve)",
            "val": evaluate_multilabel(y_val_np, val_proba, drug_columns, threshold=tuned),
            "test": evaluate_multilabel(y_test_np, test_proba, drug_columns, threshold=tuned),
        },
    }
    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_PATH, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(f"Wrote {RESULTS_PATH}")


if __name__ == "__main__":
    main()
