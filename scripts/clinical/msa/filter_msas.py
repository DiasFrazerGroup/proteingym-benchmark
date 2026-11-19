"""Final Genebass archives: proteins <= 2,700 aa, with the EVE-protocol filtered alignment stored directly.

User decisions (2026-11-19): the benchmark is capped at 2,700 aa (one full-length AlphaFold DB structure per protein;
APOB, LRP2 and PKD1 are dropped), and the archives carry only the alignment rows the models use, instead of the full
ColabFold search (up to 1.1M rows) that every alignment model then filtered itself.

The filter is the model side's own code, vendored verbatim: `msa_filter.py` here is a byte copy of
`models/<model>_filtered/src/proteingym/models/<model>_filtered/msa_filter.py` on branch clinical/model-api (git blob
MSA_FILTER_BLOB below; the script refuses to run if the file differs). Parameters = the `*_filtered` card defaults
(FILTER_DEFAULTS in msa_filter.py). The filter needs evedesign's Sequence/Sequences containers only as holders; this
script passes minimal stand-ins (`_install_sequence_standins`), so the vendored code runs unchanged without evedesign.
The original row index of every kept row is recovered through the `weights` field (the filter carries weights of kept
rows through), so the output is the exact original a3m record (header with bit score, insertions, line layout) of
every kept row, in the original order, query first.

For every archive in the index with seq_length <= --max-length: the MSA entry is replaced by the filtered a3m;
every other entry of the inner dataset (manifest, assays, sequences, structures) and of the outer archive (slices)
is copied byte for byte. Input archives are not modified. Writes the archives, a per-protein filter table and the
capped index (archive paths relative to the dataset folder).

  filter_msas.py --index datasets/clinical/genebass/genebass_fig6b.csv --archives $WORK/genebass/splits_layers \
      --msa-dir $WORK/genebass_msas --out $WORK/genebass/final_le2700 \
      --out-index datasets/clinical/genebass/genebass_le2700.csv \
      --out-filter datasets/clinical/genebass/genebass_le2700_msa_filter.tsv [--max-length 2700] [--workers N]
"""
import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import tomllib
import types
import zipfile
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
MSA_FILTER_BLOB = "d9bd2d4bf23134488010600f6ff3e96035df0b87"   # clinical/model-api f7b70ee, all *_filtered copies


def _check_vendored():
    blob = subprocess.run(["git", "hash-object", str(HERE / "msa_filter.py")], capture_output=True, text=True,
                          check=True).stdout.strip()
    if blob != MSA_FILTER_BLOB:
        sys.exit(f"msa_filter.py differs from the model side's copy (blob {blob} != {MSA_FILTER_BLOB})")


def _install_sequence_standins():
    """Minimal evedesign.sequence.{Sequence,Sequences} (the filter only stores and reads these fields)."""
    if "evedesign.sequence" in sys.modules:
        return

    class Sequence:
        def __init__(self, seq, id=None, type=None):
            self.seq, self.id, self.type_ = seq, id, type

    class Sequences:
        def __init__(self, seqs, aligned=True, type=None, weights=None, format=None):
            self.seqs, self.aligned, self.type_, self.weights, self.format_ = seqs, aligned, type, weights, format

    pkg = types.ModuleType("evedesign")
    mod = types.ModuleType("evedesign.sequence")
    mod.Sequence, mod.Sequences = Sequence, Sequences
    pkg.sequence = mod
    sys.modules["evedesign"], sys.modules["evedesign.sequence"] = pkg, mod


def records(raw: bytes):
    """[(record bytes as in the file, header line without '>', sequence with line breaks removed)] in file order."""
    out = []
    pos, n = 0, len(raw)
    if not raw.startswith(b">"):
        raise ValueError("a3m does not start with '>'")
    while pos < n:
        nxt = raw.find(b"\n>", pos)
        end = n if nxt < 0 else nxt + 1
        rec = raw[pos:end]
        head, _, body = rec.partition(b"\n")
        out.append((rec, head[1:], body.replace(b"\n", b"").replace(b"\r", b"")))
        pos = end
    return out


def filter_a3m(raw: bytes, L: int):
    """Filtered a3m bytes (exact original records of the kept rows, original order) and the filter stats."""
    _install_sequence_standins()
    sys.path.insert(0, str(HERE))
    import msa_filter as mf
    from evedesign.sequence import Sequence, Sequences

    recs = records(raw)
    seqs = Sequences(seqs=[Sequence(seq=r[2].decode("ascii")) for r in recs], aligned=True, type="protein",
                     weights=list(range(len(recs))), format="a3m")
    bits, query = mf.bitscores_from_a3m(raw)
    if len(bits) != len(recs):
        raise ValueError(f"{len(bits)} headers for {len(recs)} records")
    seqs.pg_bitscores, seqs.pg_query_name, seqs.pg_bitscore_source = bits, query, "a3m headers"
    params = {k: v for k, v in mf.FILTER_DEFAULTS.items() if k not in ("msa_filter", "a3m_headers_dir", "stats_dir")}
    new, stats = mf.filter_sequences(seqs, L, **params)
    idx = [int(w) for w in new.weights]
    if idx[0] != 0 or any(b <= a for a, b in zip(idx, idx[1:])):
        raise AssertionError("kept rows are not the query followed by increasing original indices")
    stats["params"] = params
    return b"".join(recs[i][0] for i in idx), stats, idx


def build_protein(args):
    """Filter one protein's a3m once; rewrite all its archives."""
    uniprot, L, archives, msa_dir, out_dir = args
    raw = (msa_dir / f"{uniprot}.a3m").read_bytes()
    filtered, stats, idx = filter_a3m(raw, L)
    stats.pop("sweep", None)
    res = []
    for arch in archives:
        out = out_dir / arch.name
        with zipfile.ZipFile(arch) as z, tempfile.TemporaryDirectory(dir=out_dir / ".tmp") as td:
            inner_in = Path(td) / "in.pgdata"
            with z.open("dataset.pgdata") as src, open(inner_in, "wb") as dst:
                shutil.copyfileobj(src, dst, 1 << 24)
            with zipfile.ZipFile(inner_in) as zi:
                man = tomllib.loads(zi.read("manifest.lock").decode())
                msas = man.get("msas", [])
                assert len(msas) == 1 and msas[0]["reference_sequence_name"] == uniprot, arch.name
                assert len(man.get("structures", [])) == 1, f"{arch.name}: {len(man.get('structures', []))} structures"
                entry = msas[0]["path"]
                # the archive's MSA must be the original a3m (restore_msa_headers.py): same CRC as the file
                if zi.getinfo(entry).CRC != (zipfile.crc32(raw) & 0xFFFFFFFF):
                    raise AssertionError(f"{arch.name}: archive MSA differs from {uniprot}.a3m")
                inner_out = Path(td) / "out.pgdata"
                with zipfile.ZipFile(inner_out, "w") as zo:
                    for info in zi.infolist():
                        new = zipfile.ZipInfo(info.filename, date_time=info.date_time)
                        new.compress_type = info.compress_type
                        if info.filename == entry:
                            with zo.open(new, "w", force_zip64=True) as dst:
                                dst.write(filtered)
                        else:
                            with zi.open(info) as src, zo.open(new, "w", force_zip64=True) as dst:
                                shutil.copyfileobj(src, dst, 1 << 24)
            tmp_out = out.with_suffix(".tmp")
            with zipfile.ZipFile(tmp_out, "w") as zo:
                for info in z.infolist():
                    new = zipfile.ZipInfo(info.filename, date_time=info.date_time)
                    new.compress_type = info.compress_type
                    src_f = open(inner_out, "rb") if info.filename == "dataset.pgdata" else z.open(info)
                    with src_f as src, zo.open(new, "w", force_zip64=True) as dst:
                        shutil.copyfileobj(src, dst, 1 << 24)
            os.replace(tmp_out, out)
        res.append(arch.name)
    stats.update(uniprot=uniprot, a3m_rows=len(idx), a3m_bytes_in=len(raw), a3m_bytes_out=len(filtered),
                 a3m_sha256=hashlib.sha256(filtered).hexdigest(), archives=len(res))
    print(f"{uniprot} L={L} rows {stats['rows_in']} -> {stats['N']} (t={stats['threshold']}, "
          f"{stats['relax_rule']}, {stats['quality']}) {len(res)} archives", flush=True)
    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", type=Path, required=True)
    ap.add_argument("--archives", type=Path, required=True)
    ap.add_argument("--msa-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--out-index", type=Path, required=True)
    ap.add_argument("--out-filter", type=Path, required=True)
    ap.add_argument("--max-length", type=int, default=2700)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    a = ap.parse_args()
    _check_vendored()
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / ".tmp").mkdir(exist_ok=True)
    idx = pd.read_csv(a.index)
    keep = idx[idx.seq_length <= a.max_length].copy()
    print(f"index {len(idx)} archives / {idx.uniprot_id.nunique()} proteins; <= {a.max_length} aa: {len(keep)} / "
          f"{keep.uniprot_id.nunique()}; dropped: {sorted(set(idx.entry_name) - set(keep.entry_name))}", flush=True)
    jobs = []
    for (u, L), g in keep.groupby(["uniprot_id", "seq_length"]):
        jobs.append((u, int(L), [a.archives / Path(p).name for p in g.archive], a.msa_dir, a.out))
    jobs.sort(key=lambda j: -(a.msa_dir / f"{j[0]}.a3m").stat().st_size)
    with Pool(a.workers) as pool:
        stats = list(pool.imap_unordered(build_protein, jobs))
    (a.out / ".tmp").rmdir()
    st = pd.DataFrame(stats)
    cols = ["uniprot", "L", "rows_in", "rows_fragments", "rows_indeterminate", "rows_bad_length",
            "rows_after_row_filters", "threshold", "N", "N_over_L", "Lcov", "Lcov_frac", "relax_step", "relax_rule",
            "subsampled", "quality", "a3m_bytes_in", "a3m_bytes_out", "a3m_sha256"]
    st = st[cols].sort_values("uniprot")
    a.out_filter.parent.mkdir(parents=True, exist_ok=True)
    st.to_csv(a.out_filter, sep="\t", index=False)
    m = st.set_index("uniprot")
    keep["msa_n_sequences_colabfold"] = keep.msa_n_sequences
    keep["msa_n_sequences"] = keep.uniprot_id.map(m.N).astype(int)
    keep["msa_bits_threshold"] = keep.uniprot_id.map(m.threshold)
    keep["msa_N_over_L"] = keep.uniprot_id.map(m.N_over_L)
    keep["msa_Lcov_frac"] = keep.uniprot_id.map(m.Lcov_frac)
    keep["msa_filter_rule"] = keep.uniprot_id.map(m.relax_rule)
    keep["msa_quality"] = keep.uniprot_id.map(m.quality)
    keep["archive"] = [f"{a.out.name}/{Path(p).name}" for p in keep.archive]
    keep.to_csv(a.out_index, index=False)
    print(f"written {len(keep)} archives to {a.out}, index {a.out_index}, filter table {a.out_filter}")
    print(f"a3m bytes: {st.a3m_bytes_in.sum() / 1e9:.2f} GB -> {st.a3m_bytes_out.sum() / 1e9:.2f} GB")


if __name__ == "__main__":
    main()
