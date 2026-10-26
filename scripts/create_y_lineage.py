"""Generate data/processed/y_lineage.csv and y_lineage_metadata.json from labels.csv.

Lineage strings are used exactly as in labels.csv, with ONE change: the 10
isolates with Lineage == "-" (missing) get lineage_id = -1 and are excluded from
lineage losses/metrics (never guessed). The two spellings BOV_AFRI / BOV-AFRI are
deliberately kept as SEPARATE classes -- merging them is a biological assumption
nobody has confirmed. Class ids are assigned in sorted-name order. Row order =
sample_ids.csv (validated).
"""

import json
from pathlib import Path

import pandas as pd

PROCESSED = Path("data/processed")
MISSING_VALUE = "-"


def main():
    sample_ids = pd.read_csv(PROCESSED / "sample_ids.csv")["Name"].tolist()
    labels = pd.read_csv(PROCESSED / "labels.csv")
    if labels["Name"].tolist() != sample_ids:
        raise ValueError("labels.csv row order differs from sample_ids.csv")

    lineage = labels["Lineage"].astype(str).str.strip()
    classes = sorted(set(lineage) - {MISSING_VALUE})
    class_to_id = {name: i for i, name in enumerate(classes)}
    lineage_id = lineage.map(class_to_id).fillna(-1).astype(int)

    pd.DataFrame({"Name": sample_ids, "Lineage": lineage, "lineage_id": lineage_id}).to_csv(
        PROCESSED / "y_lineage.csv", index=False)
    metadata = {
        "n_samples": len(sample_ids),
        "classes": class_to_id,
        "n_missing_masked": int((lineage_id == -1).sum()),
        "class_counts": {c: int((lineage == c).sum()) for c in classes},
        "notes": "BOV_AFRI and BOV-AFRI kept separate (unconfirmed whether the same class); "
                 "Lineage '-' masked as lineage_id=-1, never imputed.",
    }
    (PROCESSED / "y_lineage_metadata.json").write_text(json.dumps(metadata, indent=2))
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
