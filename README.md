# fennec


## Installation

```bash
git clone --recurse-submodules https://github.com/keuv-grvl/fennec.git
cd fennec/
pip install numpy  # scikit-bio requires numpy to be already installed
pip install .
```

## Usage

```python
import os, fennec

fastafile = '/path/to/file.fasta'
h5file = '/path/to/file.h5'  # every data will be stored here

# load sequences
if os.path.exists(h5file):
    seqdb = fennec.DNASequenceBank.read_hdf(h5file)
else:
    seqdb = fennec.DNASequenceBank(min_lenght=1000, chunk_size=10000, verbose=2)
    seqdb.read_fasta(fastafile)
    seqdb.to_hdf(h5file)

# define models
models_to_apply = {
    'raw_kmers4':       fennec.MaskedKmerModel(mask="1111", n_jobs=160, verbose=3),
    'raw_kmers110011':  fennec.MaskedKmerModel(mask="110011", n_jobs=160, verbose=3),
    'raw_ind15':        fennec.InterNucleotideDistanceModel(K=15, n_jobs=160, verbose=2),
    'raw_contig2vec4':  fennec.Contig2VecModel(n_jobs=160, verbose=2)
}

# apply models
for model in models_to_apply.keys():
    print(f" ø {model}")
    X = models_to_apply[model].fit_transform(seqdb)
    if model.startswith("raw_kmers"):  # if kmer count
        X = X.astype(int)
    print(f"{model} loaded (shape={X.shape})")
    X.to_hdf(h5file, model)
```

## Dependencies

The Conda environment is provided.

### System dependencies

- gcc
- libgsl-dev

### Python libraries

- scikit-bio
- bcbio-gff
- numpy
- scipy
- pandas
- gensim

### External software

- Prodigal
- FragGeneScan
- [MetaGeneAnnotator](http://metagene.cb.k.u-tokyo.ac.jp/metagene/mga_x86_64.tar.gz)
