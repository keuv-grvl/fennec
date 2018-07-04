#!/usr/bin/env python3

# source activate fennec2-dev

import gc
import os
import subprocess
import sys

import h5py

import fennec

if os.environ["CONDA_DEFAULT_ENV"] != "fennec2-dev":
    raise Warning("Not in correct conda environment")

if not fennec._utils.isinteractive():
    try:
        _, DATASET, min_length, chunk_size, overlap, n_jobs = sys.argv
        min_length, chunk_size, n_jobs = (int(min_length), int(chunk_size), int(n_jobs))
        try:
            overlap = int(overlap)
        except:  # `overlap` may be "auto"
            assert overlap == 'auto', "`overlap` must be 0+ or 'auto'"
    except:
        print(
            f"usage: {sys.argv[0]} <S,M,L> <min_length> <chunk_size> <overlap> <n_jobs>"
        )
        sys.exit(1)
else:
    DATASET = "S"
    min_length = 1000
    chunk_size = 10000
    overlap = "auto"
    n_jobs = os.cpu_count()

print(f"== Processing {DATASET} ==")

# -- variable defintions
fastafile = f"DATA/{DATASET}.Scaffolds.fasta"
h5file = f"DATA/{DATASET}.completedata.l{min_length}c{chunk_size}o{overlap}.h5"
cov_file = f"DATA/{DATASET}.Scaffolds.l{min_length}c{chunk_size}o{overlap}.csv"
force_gc = True

# -- load sequences
if os.path.exists(h5file):
    seqdb = fennec.DNASequenceBank.read_hdf(h5file)
else:
    seqdb = fennec.DNASequenceBank(
        min_length=min_length, chunk_size=chunk_size, overlap=overlap, verbose=2
    )
    seqdb.read_fasta(fastafile)
    print(seqdb)
    print(f"Saving to {h5file}")
    seqdb.to_hdf(h5file)

print(f"DNASequenceBank has {len(seqdb)} sequences.")

# -- feature model definitions
models_definitions = {
    # "label": fennec.xxxModel(param)
    "kmers110010011": fennec.MaskedKmerModel(
        mask="110010011", n_jobs=n_jobs, verbose=3
    ),
    "kmers110011": fennec.MaskedKmerModel(mask="110011", n_jobs=n_jobs, verbose=3),
    "kmers4": fennec.MaskedKmerModel(mask="1111", n_jobs=n_jobs, verbose=3),
    "kmers5": fennec.MaskedKmerModel(mask="11111", n_jobs=n_jobs, verbose=3),
    "ind15": fennec.InterNucleotideDistanceModel(K=15, n_jobs=n_jobs, verbose=3),
    "contig2vec4": fennec.Contig2VecModel(k=4, verbose=3),
    "contig2vec6": fennec.Contig2VecModel(k=6, verbose=3),
}

# -- Get coverage using GATTACA (Popic et al, 2017, doi: 10.1101/130997)
if not os.path.exists(cov_file) and os.access("./bin/gattaca", os.X_OK):
    fastafile = seqdb.to_fasta()
    cov_file = fastafile.replace(".fasta", ".csv")

    print(f"[INFO] Indexing {fastafile}")
    _ = subprocess.check_output(["./bin/gattaca", "index", "-k", "31", "-i", fastafile])

    print(f"[INFO] Running GATTACA on {fastafile}")
    _ = subprocess.check_output(
        [
            "./bin/gattaca",
            "lookup",
            "-c",
            fastafile,
            "-i",
            f"{fastafile}.k31.gatc",
            "-o",
            cov_file,
            "-m",  # median coverage of each contig in each sample (recommended)
        ]
    )

models_definitions["cov_gattaca31"] = fennec.SequenceCoverageModel(
    (cov_file,), verbose=3
)

# -- list models to apply
with h5py.File(h5file, "r") as hf:
    try:
        available_features = set(hf["rawmodel"].keys())
    except:
        available_features = set()

models_to_apply = set(models_definitions.keys()).difference(available_features)
print(
    f"Models to apply: {set(models_definitions.keys()).difference(available_features)}"
)

# -- extract features
for model in models_to_apply:
    print(f" ø {model}")
    X = models_definitions[model].fit_transform(seqdb)
    if model.startswith("kmers"):  # if kmer count
        X = X.astype(int)
    print(f"   shape: {X.shape}")
    print(f"   saving to: {h5file}")
    X.to_hdf(h5file, f"rawmodel/{model}")
    if force_gc:
        del X
        gc.collect()

print("Bye.")
