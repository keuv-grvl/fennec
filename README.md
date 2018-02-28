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
fastafile = '/home/prof/gravouil/fennec/sequences/XS.c10k.all.fna'
seqdb = fennec.DNASequenceBank(fastafile, 1000, verbose=2, nb_seq=0)
kc = fennec.MaskedKmerModel(masks=["1111", "111", "11", "110011"], n_jobs=160, verbose=2)
X = kc.fit(seqdb).transform(seqdb)
X = kc.fit_transform(seqdb)
X.shape
```
