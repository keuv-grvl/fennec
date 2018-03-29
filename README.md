# fennec


## Installation

```bash
git clone --recurse-submodules https://github.com/keuv-grvl/fennec.git
cd fennec/
for F in $(cat requirements.txt); do pip install $F; done  # install pip dependencies one by one
pip install .
```

## Usage

```python
import fennec
fastafile = '/home/prof/gravouil/fennec/sequences/S.Scaffolds.1kb+.split10kb.fasta'
seqdb = fennec.DNASequenceBank(fastafile, 1000, verbose=2, nb_seq=2500)

kc = fennec.MaskedKmerModel(mask="1111", n_jobs=160, verbose=2)
X = kc.fit(seqdb).transform(seqdb)
X.shape

ind = fennec.InterNucleotideDistanceModel(K=15, n_jobs=160, verbose=2)
Y = ind.fit(seqdb).transform(seqdb)
Y.shape

cd = fennec.CodingDensityModel(force=True, n_jobs=160, verbose=2)
cd.tools  # list of tools that will be used
Z = cd.fit(seqdb).transform(seqdb)
Z.shape

dv = fennec.Contig2VecModel(verbose=2)
dv.available_models  # list of available models
W = dv.fit(seqdb).transform(seqdb)
W.shape
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
