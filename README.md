# fennec


## Installation

```bash
git clone xxx
cd fennec
pip install .
```

## Usage

```python
import fennec
fastafile = '/home/prof/gravouil/fennec/sequences/S.Scaffolds.1kb+.split10kb.fasta'
seqdb = fennec.DNASequenceBank(fastafile, 1000, verbose=2, nb_seq=25)

kc = fennec.MaskedKmerModel(masks=["1111", "111", "11", "110011"], n_jobs=160, verbose=2)
X = kc.fit(seqdb).transform(seqdb)
X = kc.fit_transform(seqdb)
X.shape

ind = fennec.InterNucleotideDistanceModel(K=15, n_jobs=160, verbose=2)
Y = ind.fit(seqdb).transform(seqdb)
Y.shape

cd = fennec.CodingDensityModel(force=True, verbose=3)
cd.tools  # list of tools
Z = cd.fit(seqdb).transform(seqdb)
Z.shape

```


## Dependencies

The Conda environment is provided.

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
