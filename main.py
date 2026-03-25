import os
import gzip
import shutil
import urllib.request
import re
import GEOparse
import pandas as pd

# ── Configuration ─────────────────────────────────────────────────────────────
MANIFEST_PATH = "manifests/expression_studies.tsv"
GEO_CACHE_DIR = "resources/raw/geo_soft"
RESULTS_BASE = "results/transcriptome"
os.makedirs(GEO_CACHE_DIR, exist_ok=True)

# ── Load manifest ─────────────────────────────────────────────────────────────
manifest = pd.read_csv(MANIFEST_PATH, sep="\t")
# Expected columns: study_id, source, accession
STUDIES = {
    row["study_id"]: row["accession"]
    for _, row in manifest.iterrows()
    if row["source"].strip().upper() == "GEO"
}
print(f"Loaded {len(STUDIES)} GEO studies from manifest: {STUDIES}")

GEO_CACHE_DIR = "resources/raw/geo_soft"
RESULTS_BASE = "results/transcriptome"
os.makedirs(GEO_CACHE_DIR, exist_ok=True)

# ── Helpers ───────────────────────────────────────────────────────────────────


def geo_ftp_base(gse_id: str) -> str:
    """GEO FTP base URL for a given GSE accession."""
    prefix = gse_id[: len(gse_id) - 3] + "nnn"  # e.g. GSE104nnn
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse_id}"


def download_soft(gse_id: str, cache_dir: str) -> str:
    """
    Manually download the SOFT .gz file via urllib and decompress it locally.
    Returns the path to the decompressed SOFT file.
    Bypasses GEOparse's broken Windows temp-file downloader entirely.
    """
    gz_filename = f"{gse_id}_family.soft.gz"
    soft_filename = f"{gse_id}_family.soft"
    gz_path = os.path.join(cache_dir, gz_filename)
    soft_path = os.path.join(cache_dir, soft_filename)

    if os.path.exists(soft_path):
        print(f"  [CACHE] {soft_filename}")
        return soft_path

    url = f"{geo_ftp_base(gse_id)}/soft/{gz_filename}"
    print(f"  [DL] SOFT → {url}")
    urllib.request.urlretrieve(url, gz_path)

    print(f"  [DECOMPRESS] {gz_filename}")
    with gzip.open(gz_path, "rb") as f_in, open(soft_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(gz_path)  # free space; soft file is all we need

    return soft_path


def list_suppl_files(gse_id: str) -> list:
    """List supplementary filenames from GEO FTP HTML index."""
    url = f"{geo_ftp_base(gse_id)}/suppl/"
    try:
        with urllib.request.urlopen(url) as r:
            html = r.read().decode()
        return re.findall(r'href="([^"]+\.(?:txt|tsv|csv|gz)[^"]*)"', html)
    except Exception as e:
        print(f"  [WARN] Could not list supplementary files: {e}")
        return []


def download_file(url: str, dest: str):
    if os.path.exists(dest):
        print(f"  [CACHE] {os.path.basename(dest)}")
        return
    print(f"  [DL] {url}")
    urllib.request.urlretrieve(url, dest)


def decompress_gz(gz_path: str) -> str:
    """Decompress .gz in-place and return path to decompressed file."""
    out_path = gz_path[:-3]
    if not os.path.exists(out_path):
        with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    return out_path


def load_table(filepath: str) -> pd.DataFrame:
    sep = "," if filepath.endswith(".csv") else "\t"
    return pd.read_csv(filepath, sep=sep, index_col=0)


def extract_metadata(gse) -> pd.DataFrame:
    records = []
    for gsm_name, gsm in gse.gsms.items():
        row = {"sample_id": gsm_name}
        for key, val in gsm.metadata.items():
            row[key] = "; ".join(val) if isinstance(val, list) else val
        records.append(row)
    return pd.DataFrame(records).set_index("sample_id")


# ── Main pipeline ─────────────────────────────────────────────────────────────

for study_label, gse_id in STUDIES.items():
    out_dir = os.path.join(RESULTS_BASE, study_label)
    os.makedirs(out_dir, exist_ok=True)
    expr_out = os.path.join(out_dir, "expression_raw.tsv")
    meta_out = os.path.join(out_dir, "metadata.tsv")

    print(f"\n{'=' * 60}")
    print(f"Processing {study_label}: {gse_id}")
    print(f"{'=' * 60}")

    # ── Step 1: Download SOFT manually → parse with GEOparse via filepath ─────
    # Using filepath= skips GEOparse's internal downloader (broken on Windows)
    soft_path = download_soft(gse_id, GEO_CACHE_DIR)
    print("  Parsing SOFT file...")
    gse = GEOparse.get_GEO(filepath=soft_path, silent=True)

    # ── Step 2: Save metadata ─────────────────────────────────────────────────
    meta_df = extract_metadata(gse)
    meta_df.to_csv(meta_out, sep="\t")
    print(f"  Metadata saved → {meta_out}  ({meta_df.shape[0]} samples)")

    # # ── Step 3: Try series matrix expression pivot ────────────────────────────
    # expr_df = gse.pivot_samples("VALUE")
    #
    # if expr_df is not None and expr_df.shape[0] > 10:
    #     print(f"  Expression from series matrix: {expr_df.shape}")
    #     expr_df.to_csv(expr_out, sep="\t")
    #
    # else:
    #     # ── Step 4: Fall back to supplementary files (standard for RNA-seq) ──
    #     print("  Series matrix empty → scanning supplementary files on FTP...")
    #     suppl_files = list_suppl_files(gse_id)
    #
    #     priority = ["count", "expression", "expr", "tpm", "fpkm", "rpkm",
    #                 "raw", "gene", "matrix"]
    # ── Step 3: Try series matrix expression pivot ────────────────────────────
    expr_df = None

    # Check if the first sample even has a table to pivot
    first_gsm = list(gse.gsms.values())[0]

    if not first_gsm.table.empty:
        try:
            # Try the standard pivot
            expr_df = gse.pivot_samples("VALUE")

            # Validation: If it's just a tiny table, it might be junk
            if expr_df is not None and expr_df.shape[0] > 10:
                print(f"  Expression from series matrix: {expr_df.shape}")
                expr_df.to_csv(expr_out, sep="\t")
            else:
                expr_df = None  # Force fallback if table is too small

        except KeyError as e:
            print(f"  [INFO] Pivot failed (missing {e}). Likely RNA-Seq/HTS data.")
            expr_df = None
    else:
        print("  [INFO] SOFT file contains no data tables (Virtual Platform).")

    # ── Step 4: Fall back to supplementary files ──────────────────────────────
    if expr_df is None:
        print("  Scanning supplementary files on FTP as fallback...")
        suppl_files = list_suppl_files(gse_id)

        # Your existing ranking and download logic starts here...
        priority = [
            "count",
            "expression",
            "expr",
            "tpm",
            "fpkm",
            "rpkm",
            "raw",
            "gene",
            "matrix",
        ]
        ranked = sorted(
            suppl_files,
            key=lambda f: any(k in f.lower() for k in priority),
            reverse=True,
        )

        expr_df = None
        suppl_base = f"{geo_ftp_base(gse_id)}/suppl/"

        for fname in ranked:
            fname_clean = fname.split("/")[-1]
            dl_path = os.path.join(GEO_CACHE_DIR, fname_clean)
            download_file(suppl_base + fname_clean, dl_path)

            if dl_path.endswith(".gz"):
                dl_path = decompress_gz(dl_path)

            try:
                candidate = load_table(dl_path)
                if candidate.shape[1] >= 2:
                    expr_df = candidate
                    print(
                        f"  Expression from supplementary '{fname_clean}': "
                        f"{expr_df.shape}"
                    )
                    break
            except Exception as e:
                print(f"  [SKIP] {fname_clean}: {e}")
                continue

        if expr_df is not None:
            expr_df.to_csv(expr_out, sep="\t")
        else:
            print(f"  [ERROR] No parseable expression file found for {gse_id}.")

    print(f"  Expression saved → {expr_out}")

print("\nDone.")
