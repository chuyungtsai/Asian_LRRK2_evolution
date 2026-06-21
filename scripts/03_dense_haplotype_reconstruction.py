#!/usr/bin/env python3
"""Dense SNP haplotype reconstruction for LRRK2 carrier chromosomes.

This script converts phased VCF genotypes into chromosome-level dense SNP
haplotype strings, identifies the most frequent carrier haplotype, tests exact
presence among noncarrier chromosomes, and estimates the shared block over
which major carrier chromosomes are identical.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
from collections import Counter
from pathlib import Path


VARIANTS = {
    "G2385R": {
        "core_id": "rs34778348",
        "core_pos": 40363526,
        "ref": "G",
        "alt": "A",
        "windows": {
            "historical_239kb_G2385R_to_D12S2517_offset": (40363526, 40602829),
            "historical_294kb_D12S2516_to_D12S2517_offset": (40309166, 40602829),
            "239kb_centered": (40363526 - 119500, 40363526 + 119500),
        },
    },
    "R1628P": {
        "core_id": "rs33949390",
        "core_pos": 40320043,
        "ref": "G",
        "alt": "C",
        "windows": {
            "50kb_centered": (40320043 - 50000, 40320043 + 50000),
            "100kb_centered": (40320043 - 100000, 40320043 + 100000),
            "239kb_centered": (40320043 - 119500, 40320043 + 119500),
            "250kb_centered": (40320043 - 250000, 40320043 + 250000),
            "500kb_centered": (40320043 - 500000, 40320043 + 500000),
        },
    },
}


def parse_gt(sample: str) -> tuple[int, int]:
    gt = sample.split(":", 1)[0]
    if "|" not in gt:
        raise ValueError(f"Unphased genotype found: {gt}")
    a, b = gt.split("|", 1)
    if a not in {"0", "1"} or b not in {"0", "1"}:
        raise ValueError(f"Unexpected genotype found: {gt}")
    return int(a), int(b)


def load_vcf(path: Path):
    samples = []
    variants = []
    alleles_by_variant = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = fields[9:]
                continue
            chrom, pos, vid, ref, alt = fields[:5]
            hap_alleles = bytearray()
            for sample_field in fields[9:]:
                a, b = parse_gt(sample_field)
                hap_alleles.append(ord("0") + a)
                hap_alleles.append(ord("0") + b)
            variants.append({"chrom": chrom, "pos": int(pos), "id": vid, "ref": ref, "alt": alt})
            alleles_by_variant.append(bytes(hap_alleles))
    return samples, variants, alleles_by_variant


def variant_indices(variants, start: int, end: int) -> list[int]:
    return [i for i, v in enumerate(variants) if start <= int(v["pos"]) <= end]


def hap_key(alleles_by_variant, idxs: list[int], hap_index: int) -> bytes:
    return bytes(alleles_by_variant[i][hap_index] for i in idxs)


def summarize_window(label, start, end, variants, alleles_by_variant, core_idx, carrier_haps, noncarrier_haps):
    idxs = variant_indices(variants, start, end)
    if not idxs:
        raise ValueError(f"No variants in window {label}: {start}-{end}")

    carrier_keys = [hap_key(alleles_by_variant, idxs, h) for h in carrier_haps]
    counts = Counter(carrier_keys)
    major_key, major_count = counts.most_common(1)[0]
    exact_noncarrier_count = sum(1 for h in noncarrier_haps if hap_key(alleles_by_variant, idxs, h) == major_key)
    major_haps = [h for h, key in zip(carrier_haps, carrier_keys) if key == major_key]

    left_shared = core_idx
    for i in range(core_idx, -1, -1):
        if len({alleles_by_variant[i][h] for h in major_haps}) == 1:
            left_shared = i
        else:
            break

    right_shared = core_idx
    for i in range(core_idx, len(variants)):
        if len({alleles_by_variant[i][h] for h in major_haps}) == 1:
            right_shared = i
        else:
            break

    return {
        "window": label,
        "start": start,
        "end": end,
        "variant_count": len(idxs),
        "carrier_chromosomes": len(carrier_haps),
        "noncarrier_chromosomes": len(noncarrier_haps),
        "exact_major_carrier_hap_count": major_count,
        "exact_major_carrier_hap_freq": major_count / len(carrier_haps),
        "exact_same_hap_noncarrier_count": exact_noncarrier_count,
        "exact_same_hap_noncarrier_freq": exact_noncarrier_count / len(noncarrier_haps),
        "major_hap_members": ",".join(str(h) for h in major_haps),
        "shared_left_pos": variants[left_shared]["pos"],
        "shared_right_pos": variants[right_shared]["pos"],
        "shared_size_bp": variants[right_shared]["pos"] - variants[left_shared]["pos"] + 1,
    }


def write_tsv(path: Path, rows: list[dict]):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant", choices=sorted(VARIANTS), required=True)
    parser.add_argument("--vcf", default=os.environ.get("TWB_LRRK2_PHASED_VCF"))
    parser.add_argument("--out-dir", default=os.environ.get("RESULT_DIR", "results/publication_lrrk2"))
    args = parser.parse_args()

    if not args.vcf:
        raise SystemExit("Provide --vcf or set TWB_LRRK2_PHASED_VCF.")

    spec = VARIANTS[args.variant]
    out_dir = Path(args.out_dir) / "dense_haplotype"
    out_dir.mkdir(parents=True, exist_ok=True)

    samples, variants, alleles_by_variant = load_vcf(Path(args.vcf))
    core_idx = next(i for i, v in enumerate(variants) if v["id"] == spec["core_id"] or v["pos"] == spec["core_pos"])
    core = variants[core_idx]
    if core["pos"] != spec["core_pos"] or core["ref"] != spec["ref"] or core["alt"] != spec["alt"]:
        raise RuntimeError(f"Unexpected core variant metadata: {core}")

    core_alleles = alleles_by_variant[core_idx]
    carrier_haps = [i for i, allele in enumerate(core_alleles) if chr(allele) == "1"]
    noncarrier_haps = [i for i, allele in enumerate(core_alleles) if chr(allele) == "0"]

    rows = [
        {"variant": args.variant, "core_id": spec["core_id"], "core_pos": spec["core_pos"], **summarize_window(
            label, start, end, variants, alleles_by_variant, core_idx, carrier_haps, noncarrier_haps
        )}
        for label, (start, end) in spec["windows"].items()
    ]
    write_tsv(out_dir / f"{args.variant}_dense_haplotype_window_summary.tsv", rows)
    print(out_dir / f"{args.variant}_dense_haplotype_window_summary.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

