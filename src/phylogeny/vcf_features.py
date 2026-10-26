"""Build a VCF-derived feature matrix from the genome-wide SNP matrix.

No curated mutation list is used. SNPs are selected using the TRAIN split only
(so val/test never influence which columns exist):
  1. keep SNPs with train prevalence in [MIN_PREVALENCE, 1 - MIN_PREVALENCE];
  2. per drug (>= MIN_TRAIN_POSITIVES train positives), chi-square association
     with that drug's label on train; keep the TOP_PER_DRUG best SNPs;
  3. final columns = union across drugs, in genomic-position order.
Column identity is (pos, alt) from snp_matrix_sites.csv; no gene annotation is
inferred here. Writes data/processed/vcf_features.npz (+ _sites.csv, _report.json).
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.feature_selection import chi2

PROCESSED = Path("data/processed")
MIN_PREVALENCE = 0.01
TOP_PER_DRUG = 300
MIN_TRAIN_POSITIVES = 5


def main(top_per_drug=TOP_PER_DRUG, tag=""):
    """tag suffixes the output filenames (e.g. '_k1000') so variants don't overwrite the default."""
    sample_ids = pd.read_csv(PROCESSED / "sample_ids.csv")["Name"].tolist()
    y_df = pd.read_csv(PROCESSED / "y_amr.csv")
    splits = pd.read_csv(PROCESSED / "splits.csv")
    if y_df["Name"].tolist() != sample_ids or splits["Name"].tolist() != sample_ids:
        raise ValueError("y_amr.csv / splits.csv order differs from sample_ids.csv")

    snp = sp.load_npz(PROCESSED / "snp_matrix.npz").tocsr()
    sites = pd.read_csv(PROCESSED / "snp_matrix_sites.csv")
    sites = sites[sites["kept_in_core"]].reset_index(drop=True)
    if snp.shape != (len(sample_ids), len(sites)):
        raise ValueError(f"snp_matrix {snp.shape} does not match {len(sample_ids)} samples x {len(sites)} core sites")

    train = (splits["split"] == "train").to_numpy()
    snp_train = snp[train]
    prevalence = np.asarray(snp_train.mean(axis=0)).ravel()
    candidates = np.where((prevalence >= MIN_PREVALENCE) & (prevalence <= 1 - MIN_PREVALENCE))[0]
    snp_cand_train = snp_train[:, candidates]

    selected, per_drug_report = set(), {}
    for drug in y_df.columns[1:]:
        y_train = y_df[drug].to_numpy()[train]
        n_pos = int(y_train.sum())
        if n_pos < MIN_TRAIN_POSITIVES:
            per_drug_report[drug] = {"train_positives": n_pos, "selected": 0, "note": "too few positives"}
            continue
        stat, _ = chi2(snp_cand_train, y_train)
        stat = np.nan_to_num(stat)
        top = candidates[np.argsort(-stat)[:top_per_drug]]
        selected.update(top.tolist())
        per_drug_report[drug] = {"train_positives": n_pos, "selected": len(top)}

    columns = np.array(sorted(selected))  # core-site order == genomic order
    X = snp[:, columns].astype(np.uint8)
    sp.save_npz(PROCESSED / f"vcf_features{tag}.npz", X)
    sites.iloc[columns].to_csv(PROCESSED / f"vcf_features{tag}_sites.csv", index=False)
    report = {
        "n_samples": int(X.shape[0]), "n_features": int(X.shape[1]),
        "n_snps_core": int(len(sites)), "n_candidates_after_prevalence_filter": int(len(candidates)),
        "min_prevalence": MIN_PREVALENCE, "top_per_drug": top_per_drug,
        "min_train_positives": MIN_TRAIN_POSITIVES, "selection_split": "train only",
        "per_drug": per_drug_report,
    }
    (PROCESSED / f"vcf_features{tag}_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps({k: v for k, v in report.items() if k != "per_drug"}, indent=2))
    print(json.dumps(per_drug_report))


if __name__ == "__main__":
    import sys
    main(*(int(sys.argv[1]), sys.argv[2]) if len(sys.argv) == 3 else ())
