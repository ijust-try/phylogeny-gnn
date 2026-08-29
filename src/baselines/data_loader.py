"""Load and align Person 1's processed tabular artifacts for baseline training.

Strictly reuses data/processed/sample_ids.csv, features.csv, y_amr.csv, and
splits.csv as-is -- no relabeling, no re-deriving targets, no changes to the
fixed train/val/test assignment. All four files are validated to share the
exact same row order before anything is returned.
"""

from pathlib import Path

import numpy as np
import pandas as pd

PROCESSED_DIR = Path("data/processed")
SAMPLE_IDS_CSV = PROCESSED_DIR / "sample_ids.csv"
FEATURES_CSV = PROCESSED_DIR / "features.csv"
Y_AMR_CSV = PROCESSED_DIR / "y_amr.csv"
SPLITS_CSV = PROCESSED_DIR / "splits.csv"

DRUG_COLUMNS = ["RIF", "INH", "EMB", "PZA", "STM", "LEV", "CAP", "ETH", "LZD"]
VALID_SPLITS = {"train", "val", "test"}


def load_aligned_data():
    """Load features, y_amr, and split assignment, strictly aligned to canonical sample order.

    Returns a dict with keys: sample_ids (list[str]), feature_names (list[str]),
    X (float32 ndarray, N x 157), y (float32 ndarray, N x 9), split (ndarray[str], N),
    drug_columns (list[str]).

    Raises ValueError if any file's row order or content does not match
    sample_ids.csv, or if values are not binary.
    """
    sample_ids = pd.read_csv(SAMPLE_IDS_CSV)["Name"].tolist()

    features_df = pd.read_csv(FEATURES_CSV)
    y_amr_df = pd.read_csv(Y_AMR_CSV)
    splits_df = pd.read_csv(SPLITS_CSV)

    if features_df["Name"].tolist() != sample_ids:
        raise ValueError("features.csv row order does not match sample_ids.csv")
    if y_amr_df["Name"].tolist() != sample_ids:
        raise ValueError("y_amr.csv row order does not match sample_ids.csv")
    if splits_df["Name"].tolist() != sample_ids:
        raise ValueError("splits.csv row order does not match sample_ids.csv")

    feature_names = [c for c in features_df.columns if c != "Name"]
    X = features_df[feature_names].to_numpy(dtype=np.float32)
    if not np.isin(X, [0.0, 1.0]).all():
        raise ValueError("features.csv contains non-binary values")

    y_amr_drug_columns = list(y_amr_df.columns[1:])
    if y_amr_drug_columns != DRUG_COLUMNS:
        raise ValueError(
            f"y_amr.csv drug columns {y_amr_drug_columns} != expected {DRUG_COLUMNS}"
        )
    y = y_amr_df[DRUG_COLUMNS].to_numpy(dtype=np.float32)
    if not np.isin(y, [0.0, 1.0]).all():
        raise ValueError("y_amr.csv contains non-binary values")

    split = splits_df["split"].to_numpy()
    observed_splits = set(np.unique(split))
    if not observed_splits <= VALID_SPLITS:
        raise ValueError(f"Unexpected split labels: {observed_splits - VALID_SPLITS}")

    return {
        "sample_ids": sample_ids,
        "feature_names": feature_names,
        "X": X,
        "y": y,
        "split": split,
        "drug_columns": DRUG_COLUMNS,
    }


def split_masks(data):
    """Boolean masks for train/val/test over the canonical row order."""
    split = data["split"]
    return {name: (split == name) for name in ("train", "val", "test")}
