"""ColabFold MSA search on UniRef100 with all post-alignment filtering disabled (PG2 MSA protocol).

PG2 MSAs: "ColabFold MSA built on UniRef100 (January 2026 release) with no post-alignment filtering". This runs the
monomer search of ColabFold's `colabfold_search` (colabfold 1.5.5, colabfold/mmseqs/search.py,
mmseqs_search_monomer; MIT licence) with `--filter 0`, against a plain MMseqs2 database of UniRef100 2026_01:

  1. search       qdb vs UniRef100: 3 profile iterations, -e 0.1, -s 8, -a, --max-seqs MAX_SEQS
  2. align        realign the hits with the final query profile: -e 10, --max-accept 1e6, --alt-ali 10, -a
  3. filterresult filter=0 values: --qid 0 --qsc -20 --diff 0 --max-seq-id 1.0 --filter-min-enable 100
  4. result2msa   a3m (--msa-format-mode 6), --filter-msa 0 (the ColabFold MSA filter switched off)
  5. unpackdb     one a3m per query

Differences from colabfold_search, forced by the database:
  - no `expandaln`: ColabFold expands UniRef30 cluster hits to their members (`_aln` cluster DB); UniRef100 is not
    clustered, so every member is searched directly;
  - no environmental database (colabfold_envdb): the PG2 protocol names UniRef100 only;
  - --max-seqs: ColabFold keeps 10,000 prefilter hits and gets depth from the cluster expansion. With no expansion
    that would cap the MSA at ~10k, while PG2 ColabFold MSAs reach 1,000,000 sequences (msa_survey: median 168k),
    so the default here is 1,000,000. Set with --max-seqs.

Input: FASTA of query sequences (header = protein ID). Output: <out>/<ID>.a3m, query first, header renamed to the ID.
Run: python colabfold_search_uniref100.py queries.fasta UNIREF100_DB OUT_DIR --mmseqs <bin> --threads 64
"""
import argparse
import shutil
import subprocess
import time
from pathlib import Path


def run(mmseqs, *args):
    cmd = [str(mmseqs)] + [str(a) for a in args]
    print(f"[{time.strftime('%H:%M:%S')}] $ {' '.join(cmd)}", flush=True)
    subprocess.check_call(cmd)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("queries", type=Path)
    ap.add_argument("db", type=Path, help="MMseqs2 database prefix of UniRef100")
    ap.add_argument("out", type=Path)
    ap.add_argument("--mmseqs", default="mmseqs")
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--max-seqs", type=int, default=1_000_000)
    ap.add_argument("--s", type=float, default=8.0, help="sensitivity (ColabFold default 8)")
    ap.add_argument("--split-memory-limit", default="0", help="mmseqs --split-memory-limit, e.g. 500G (0 = auto)")
    ap.add_argument("--work", type=Path, help="work dir (default <out>/work, removed at the end)")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    w = a.work or a.out / "work"
    w.mkdir(parents=True, exist_ok=True)
    m, t = a.mmseqs, ["--threads", a.threads]

    ids = [line[1:].split()[0] for line in a.queries.read_text().splitlines() if line.startswith(">")]
    run(m, "createdb", a.queries, w / "qdb", "--shuffle", 0)
    # filter=0 settings of colabfold mmseqs_search_monomer
    run(m, "search", w / "qdb", a.db, w / "res", w / "tmp", *t, "--num-iterations", 3, "--db-load-mode", 0, "-a",
        "-e", 0.1, "--max-seqs", a.max_seqs, "-s", f"{a.s:.1f}", "--split-memory-limit", a.split_memory_limit)
    run(m, "mvdb", w / "tmp/latest/profile_1", w / "prof_res")
    run(m, "lndb", w / "qdb_h", w / "prof_res_h")
    run(m, "align", w / "prof_res", a.db, w / "res", w / "res_realign", *t, "--db-load-mode", 0, "-e", 10,
        "--max-accept", 1_000_000, "--alt-ali", 10, "-a")
    run(m, "filterresult", w / "qdb", a.db, w / "res_realign", w / "res_filt", *t, "--db-load-mode", 0, "--qid", 0,
        "--qsc", -20.0, "--diff", 0, "--max-seq-id", 1.0, "--filter-min-enable", 100)
    run(m, "result2msa", w / "qdb", a.db, w / "res_filt", w / "msa", *t, "--msa-format-mode", 6, "--db-load-mode", 0,
        "--filter-msa", 0, "--filter-min-enable", 1000, "--diff", 3000, "--qid", "0.0,0.2,0.4,0.6,0.8,1.0",
        "--qsc", 0, "--max-seq-id", 0.95)
    run(m, "unpackdb", w / "msa", w / "unpacked", "--unpack-name-mode", 0, "--unpack-suffix", ".a3m")

    for i, pid in enumerate(ids):                       # unpackdb names files by query index (createdb order)
        src = w / "unpacked" / f"{i}.a3m"
        text = src.read_text().lstrip("\x00")
        lines = text.split("\n")
        assert lines[0].startswith(">"), (pid, lines[0][:50])
        lines[0] = f">{pid}"                            # query first; header = protein ID
        (a.out / f"{pid}.a3m").write_text("\n".join(lines).rstrip("\x00\n") + "\n")
        n = sum(1 for x in lines if x.startswith(">"))
        print(f"{pid}\t{n} sequences", flush=True)
    shutil.rmtree(w)


if __name__ == "__main__":
    main()
