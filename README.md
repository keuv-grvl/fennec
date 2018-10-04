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
    # the HDF5 file already exists, load data from it
    seqdb = fennec.DNASequenceBank.read_hdf(h5file)
else:
    # otherwise, parse the FASTA file
    seqdb = fennec.DNASequenceBank(min_length=1000, chunk_size=10000, overlap=0, verbose=2)
    seqdb.read_fasta(fastafile)
    seqdb.to_hdf(h5file)

# define models
models_to_apply = {
    'kmers4':       fennec.MaskedKmerModel(mask="1111"),
    'kmers110011':  fennec.MaskedKmerModel(mask="110011"),
    'ind15':        fennec.InterNucleotideDistanceModel(K=15),
    'contig2vec4':  fennec.Contig2VecModel(k=4, modelfile='urq')
}

# apply models
for model in models_to_apply.keys():
    print(f" - applying {model}")
    X = models_to_apply[model].fit_transform(seqdb)
    print(f"{model} loaded (shape={X.shape})")
    X.to_hdf(h5file, model)

# load model from HDF5 file
raw_models, id_list, mustlink_mat = load_models(h5file, wanted_models)

# print number of dimension per model
print([(i, d.shape[1]) for i, d in raw_models.items()])

# merge models
kpca_params = {
    "inertia": 0.85,  # kernel PCA inertia to keep
    "n_jobs": 8,      # number of jobs
    "verbose": 3,     # verbosity level
    "t": 0.33         # proportion of data to be sampled for training
}

D, pca_components, pca_explained_variance_ratio, n_comp = merge_models(
    raw_models, id_list, kpca_params
)

# first component origins
pcacomp_to_model(D[0], raw_models, 0, outfile="pca_components_origin_comp0.csv")
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

- [Prodigal](https://github.com/hyattpd/Prodigal)
- [FragGeneScan](https://sourceforge.net/projects/fraggenescan/)
- [MetaGeneAnnotator](http://metagene.cb.k.u-tokyo.ac.jp/metagene/mga_x86_64.tar.gz)
