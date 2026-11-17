"""Put the original ColabFold a3m (with its headers) back into Genebass archives.

proteingym-base `MSA.dump` rewrites every header as `seq_<i>` (and wraps lines at 80 characters), so the archives
lose the MMseqs2 hit annotations of the ColabFold a3m: name, bit score, identity, e-value, query/target ranges. The
alignment models need the bit scores to follow the EVE protocol (bit-score threshold selection; see
docs/clinical/genebass.md "MSA headers" and the model branch's docs/clinical/models_long_proteins.md §3.5-3.6).

For every archive: the archive's MSA and the original `<msa-dir>/<UNIPROT>.a3m` are streamed side by side and must
contain the same sequences in the same order (residues, insertions and gaps, exact). Headers must be unique (the
loader keys sequences by header). Then the archive is rewritten with the original file in place of the dumped one;
every other entry (manifest, assays, sequences, structures, slices) is copied byte for byte. The input archives are
not modified.

  restore_msa_headers.py <archives...> --msa-dir $WORK/genebass_msas --out $WORK/genebass/splits_layers_hdr
"""
import argparse
import os
import shutil
import tempfile
import tomllib
import zipfile
from multiprocessing import Pool
from pathlib import Path


def records(lines):
    """(header, sequence) from a FASTA/a3m line stream (bytes), sequences joined across wrapped lines."""
    head, seq = None, []
    for line in lines:
        line = line.rstrip(b"\r\n")
        if line.startswith(b">"):
            if head is not None:
                yield head, b"".join(seq)
            head, seq = line[1:], []
        elif line:
            seq.append(line)
    if head is not None:
        yield head, b"".join(seq)


def restore(args):
    path, msa_dir, out_dir, tmp_dir = args
    out = out_dir / path.name
    with zipfile.ZipFile(path) as z:
        outer_infos = z.infolist()
        with tempfile.TemporaryDirectory(dir=tmp_dir) as td:
            inner_in = Path(td) / "in.pgdata"
            with z.open("dataset.pgdata") as src, open(inner_in, "wb") as dst:
                shutil.copyfileobj(src, dst, 1 << 24)
            with zipfile.ZipFile(inner_in) as zi:
                man = tomllib.loads(zi.read("manifest.lock").decode())
                msas = man.get("msas", [])
                if len(msas) != 1:
                    return dict(archive=path.name, status=f"skipped: {len(msas)} msas")
                entry, uniprot = msas[0]["path"], msas[0]["reference_sequence_name"]
                orig = msa_dir / f"{uniprot}.a3m"
                n, heads = 0, set()
                with zi.open(entry) as fa, open(orig, "rb") as fo:
                    ra, ro = records(fa), records(fo)
                    for (ha, sa), (ho, so) in zip(ra, ro):
                        if sa != so:
                            raise AssertionError(f"{path.name}: sequence {n} differs from {orig}")
                        hh = hash(ho)
                        if hh in heads:
                            raise AssertionError(f"{path.name}: duplicate header {ho[:80]!r} at row {n}")
                        heads.add(hh)
                        n += 1
                    if next(ra, None) is not None or next(ro, None) is not None:
                        raise AssertionError(f"{path.name}: different number of sequences (stopped at {n})")
                inner_out = Path(td) / "out.pgdata"
                with zipfile.ZipFile(inner_out, "w") as zo:
                    for info in zi.infolist():
                        if info.filename == entry:
                            zi_info = zipfile.ZipInfo(info.filename, date_time=info.date_time)
                            zi_info.compress_type = info.compress_type
                            with open(orig, "rb") as src, zo.open(zi_info, "w", force_zip64=True) as dst:
                                shutil.copyfileobj(src, dst, 1 << 24)
                        else:
                            with zi.open(info) as src, zo.open(info, "w", force_zip64=True) as dst:
                                shutil.copyfileobj(src, dst, 1 << 24)
            tmp_out = out.with_suffix(".tmp")
            with zipfile.ZipFile(tmp_out, "w") as zo:
                for info in outer_infos:
                    src_f = open(inner_out, "rb") if info.filename == "dataset.pgdata" else z.open(info)
                    new = zipfile.ZipInfo(info.filename, date_time=info.date_time)
                    new.compress_type = info.compress_type
                    with src_f as src, zo.open(new, "w", force_zip64=True) as dst:
                        shutil.copyfileobj(src, dst, 1 << 24)
            os.replace(tmp_out, out)
    return dict(archive=path.name, status="ok", uniprot=uniprot, n_sequences=n)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("archives", nargs="+", type=Path)
    ap.add_argument("--msa-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    tmp = a.out / ".tmp"
    tmp.mkdir(exist_ok=True)
    jobs = sorted(((p, a.msa_dir, a.out, tmp) for p in a.archives), key=lambda j: -j[0].stat().st_size)
    with Pool(a.workers) as pool:
        for r in pool.imap_unordered(restore, jobs):
            print(r, flush=True)
    tmp.rmdir()


if __name__ == "__main__":
    main()
