"""Terminal demo: raw VCF -> SNP feature vector -> phylogeny-aware GNN -> per-drug resistance + lineage.

Usage (from the project root):
    python -m src.demo.predict_demo SRR1062844 ERR3804168
    python -m src.demo.predict_demo --vcf data/raw/Afro_TB/AFRO_TB_VCF/ERR6362263_MT.vcf.gz

The isolate's 932-SNP feature vector is built directly from the VCF
(src/phylogeny/parse_vcf.py) and swapped into the isolate's node of the
pre-built phylogenetic graph; the saved GraphSAGE weights then run one
forward pass. The graph is transductive, so the isolate must already be a
node (its sample ID must be in sample_ids.csv). Ground truth, split and the
lineage prediction are looked up, not computed live (lineage weights were not saved).
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch

from ..gnn.data import load_gnn_data
from ..gnn.model import MultiTaskGCN
from ..phylogeny.parse_vcf import extract_snp_calls
from ..phylogeny.vcf_paths import vcf_path_for

MODEL_DIR = Path("models")
PROCESSED = Path("data/processed")


def load_model_and_context():
    """Load saved GNN weights, VCF feature columns, graph data and lineage lookup."""
    meta = json.loads((MODEL_DIR / "amr_gnn_meta.json").read_text())
    arch = meta["arch"]
    model = MultiTaskGCN(arch["in_dim"], arch["hidden_dim"], arch["n_drugs"], dropout=arch["dropout"],
                         skip_connection=arch["skip_connection"], conv_type=arch["conv_type"],
                         n_layers=arch["n_layers"], head_hidden=arch["head_hidden"])
    model.load_state_dict(torch.load(MODEL_DIR / "amr_gnn_vcf_k300.pt"))
    model.eval()
    sites = pd.read_csv(PROCESSED / "vcf_features_sites.csv")
    site_index = {(int(p), a): j for j, (p, a) in enumerate(zip(sites["pos"], sites["alt"]))}
    X_stored = sp.load_npz(PROCESSED / "vcf_features.npz").toarray().astype(np.float32)
    data = load_gnn_data()
    splits = pd.read_csv(PROCESSED / "splits.csv").set_index("Name")["split"]
    lineage = pd.read_csv(MODEL_DIR / "lineage_gnn_predictions.csv").set_index("Name")
    return dict(model=model, meta=meta, site_index=site_index, X=torch.tensor(X_stored), data=data,
                splits=splits, lineage=lineage)


def vcf_to_feature_vector(vcf_path, site_index):
    """Parse one VCF and return the 932-dim binary vector plus the number of SNP calls in the VCF."""
    calls, _ = extract_snp_calls(vcf_path)
    vector = np.zeros(len(site_index), dtype=np.float32)
    for call in calls:
        if call in site_index:
            vector[site_index[call]] = 1.0
    return vector, len(calls)


def predict_sample(ctx, sample_id, vcf_path=None):
    """Predict for one isolate from its VCF; print an input/output report."""
    data, meta = ctx["data"], ctx["meta"]
    if sample_id not in data["sample_ids"]:
        raise SystemExit(f"{sample_id} is not a node in the phylogenetic graph (not in sample_ids.csv).")
    node = data["sample_ids"].index(sample_id)
    vcf_path = Path(vcf_path) if vcf_path else vcf_path_for(sample_id)
    vector, n_calls = vcf_to_feature_vector(vcf_path, ctx["site_index"])

    X = ctx["X"].clone()
    matches_stored = bool(np.array_equal(X[node].numpy(), vector))
    X[node] = torch.tensor(vector)
    with torch.no_grad():
        proba = torch.sigmoid(ctx["model"](X, data["edge_index"], data["edge_weight"]))[node].numpy()

    print("=" * 66)
    print(f"INPUT   : {vcf_path.name}")
    print(f"          {n_calls} SNP calls parsed -> {int(vector.sum())} of {len(vector)} model SNPs present")
    print(f"          vector identical to stored feature row: {matches_stored}")
    print(f"          graph node {node}, split = {ctx['splits'][sample_id]}")
    print(f"{'Drug':<5} {'P(resist)':>9} {'thresh':>7}  {'Predicted':<12} {'Truth':<12}")
    for j, drug in enumerate(meta["drugs"]):
        threshold = meta["thresholds"][drug]
        truth = int(data["y"][node, j])
        predicted = "RESISTANT" if proba[j] >= threshold else "susceptible"
        mark = "" if (proba[j] >= threshold) == bool(truth) else "  <-- MISMATCH"
        print(f"{drug:<5} {proba[j]:>9.3f} {threshold:>7.2f}  {predicted:<12} "
              f"{'RESISTANT' if truth else 'susceptible':<12}{mark}")
    lin = ctx["lineage"].loc[sample_id]
    print(f"LINEAGE : predicted {lin['pred']} | true {lin['true']}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sample_ids", nargs="*", help="sample IDs, e.g. SRR1062844")
    parser.add_argument("--vcf", help="explicit VCF path (sample ID is taken from the filename)")
    args = parser.parse_args()
    ctx = load_model_and_context()
    if args.vcf:
        predict_sample(ctx, Path(args.vcf).name.split("_MT")[0], args.vcf)
    for sample_id in args.sample_ids:
        predict_sample(ctx, sample_id)


if __name__ == "__main__":
    main()
