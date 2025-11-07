import argparse
import csv
import re
from pathlib import Path


TERM_NAMES = [
    ("gauss1", "gauss(o=0,_w=0.5,_c=8)", -0.035579),
    ("gauss2", "gauss(o=3,_w=2,_c=8)", -0.005156),
    ("repulsion", "repulsion(o=0,_c=8)", 0.840245),
    ("hydrophobic", "hydrophobic(g=0.5,_b=1.5,_c=8)", -0.035069),
    ("hydrogen_bond", "non_dir_h_bond(g=-0.7,_b=0,_c=8)", -0.587439),
    ("torsion", "num_tors_div", 1.923),
]

LOG_RAW_PATTERN = re.compile(r"^##\s+(.+)$")


def parse_log(log_path: Path):
    affinity = None
    intramolecular = None
    raw_values = None

    with log_path.open() as fh:
        for line in fh:
            if line.startswith("Affinity:"):
                affinity = float(line.split()[1])
            elif line.startswith("Intramolecular energy:"):
                intramolecular = float(line.split()[2])
            else:
                match = LOG_RAW_PATTERN.match(line)
                if match:
                    numbers = match.group(1).strip().split()
                    if numbers and all(_is_float(n) for n in numbers):
                        raw_values = [float(n) for n in numbers]
                        break

    if raw_values is None or len(raw_values) < len(TERM_NAMES):
        raise ValueError(
            f"Could not locate all raw term values in log '{log_path}'. "
            "Please ensure the log was generated with '--score_only'."
        )

    return affinity, intramolecular, raw_values[: len(TERM_NAMES)]


def _is_float(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def build_summary(raw_values, intramolecular):
    entries = {}
    total_raw = sum(raw_values)
    total_weighted = 0.0

    for (term_key, _term_label, coeff), raw in zip(TERM_NAMES, raw_values):
        weighted = raw * coeff
        entries[f"{term_key}_raw"] = raw
        entries[f"{term_key}_coeff"] = coeff
        entries[f"{term_key}_weighted"] = weighted
        entries[term_key] = weighted
        total_weighted += weighted

    entries["SMINA_Intramolecular"] = intramolecular
    entries["SMINA_Total_Raw"] = total_raw
    entries["SMINA_Total_Weighted"] = total_weighted
    entries["gauss1_per_atom"] = ""
    entries["gauss2_per_atom"] = ""
    entries["repulsion_per_atom"] = ""
    entries["hydrophobic_per_atom"] = ""
    entries["hydrogen_bond_per_atom"] = ""

    return entries


def write_summary_csv(output_path: Path, summary: dict, affinity: float):
    field_order = [
        "SMINA_Intramolecular",
        "SMINA_Total_Raw",
        "SMINA_Total_Weighted",
        "gauss1_raw",
        "gauss1_coeff",
        "gauss1_weighted",
        "gauss2_raw",
        "gauss2_coeff",
        "gauss2_weighted",
        "repulsion_raw",
        "repulsion_coeff",
        "repulsion_weighted",
        "hydrophobic_raw",
        "hydrophobic_coeff",
        "hydrophobic_weighted",
        "hydrogen_bond_raw",
        "hydrogen_bond_coeff",
        "hydrogen_bond_weighted",
        "torsion_raw",
        "torsion_coeff",
        "torsion_weighted",
        "gauss1_per_atom",
        "gauss2_per_atom",
        "repulsion_per_atom",
        "hydrophobic_per_atom",
        "hydrogen_bond_per_atom",
        "gauss1",
        "gauss2",
        "repulsion",
        "hydrophobic",
        "hydrogen_bond",
        "Affinity",
    ]

    summary_with_affinity = dict(summary)
    summary_with_affinity["Affinity"] = affinity

    # Ensure torsion entries exist even if zero raw value
    for field in ("torsion_raw", "torsion_coeff", "torsion_weighted", "torsion"):
        summary_with_affinity.setdefault(field, 0.0 if field != "torsion_coeff" else TERM_NAMES[-1][2])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=field_order)
        writer.writeheader()
        writer.writerow({key: summary_with_affinity.get(key, "") for key in field_order})


def main():
    parser = argparse.ArgumentParser(
        description="Summarize SMINA score-only output into aggregate component CSV",
    )
    parser.add_argument("--log", required=True, help="Path to the SMINA score log")
    parser.add_argument(
        "--output",
        required=True,
        help="Destination CSV file for aggregate component values",
    )

    args = parser.parse_args()
    log_path = Path(args.log)
    output_path = Path(args.output)

    affinity, intramolecular, raw_values = parse_log(log_path)
    summary = build_summary(raw_values, intramolecular)
    write_summary_csv(output_path, summary, affinity)


if __name__ == "__main__":
    main()

