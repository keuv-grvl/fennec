# fennec

## Installation

```bash
git clone --recurse-submodules https://github.com/keuv-grvl/fennec.git
cd fennec/
pip install numpy  # scikit-bio requires numpy to be already installed
pip install .
```

By default, VBGMM will use 64 CPU or your number of CPU. You may override this with:

```bash
N_RTHREADS=24 pip install .
```

## Usage

```python
import os, fennec

fastafile = '/path/to/file.fasta'  # your contigs
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
    'raw_kmers4':       fennec.MaskedKmerModel(mask="1111"),
    'raw_kmers110011':  fennec.MaskedKmerModel(mask="110011"),
    'raw_ind15':        fennec.InterNucleotideDistanceModel(K=15),
    'raw_contig2vec4':  fennec.Contig2VecModel(k=4, modelfile='urq')
}

# apply models
for model in models_to_apply.keys():
    print(f" - applying {model}")
    X = models_to_apply[model].fit_transform(seqdb)
    print(f"{model} loaded (shape={X.shape})")
    X.to_hdf(h5file, model)
```

## Dependencies

The Conda environment is provided.

```bash
conda env create --name fennec-env --file conda-env.yml
source activate fennec-env
```

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
