
import logging

from sklearn.base import BaseEstimator, TransformerMixin

from ._utils import _print_progressbar

#-------------------------------------------------------------------------------

# Functions called by multiprocessing.Pool().map must be declared out of the class
# See https://stackoverflow.com/questions/32321324/pool-within-a-class-in-python

# src: http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
_BASE_COMPLEMENT = {
    "A": "T", "T": "A", "G": "C", "C": "G",
    "Y": "R", "R": "Y", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "D": "H", "H": "D", "V": "B",
    "N": "N",
    "x": "x" }

def _revcomp(dna):
    '''
    Get the reference sequence for a dna sequence.
    Reference sequence is the lowest (alphabetically) between the sequence
    and its reversed-complemented counterpart.
    Example:
      "CTAA" returns ('C, 'T', 'A', 'A')
      "TTAG" returns ('C, 'T', 'A', 'A')
    '''
    dna = tuple(dna)
    rcdna = tuple(map(lambda x: x.replace(x, _BASE_COMPLEMENT[x]), dna[::-1]))
    if dna < rcdna:
        return(dna)
    else:
        return(rcdna)

def _maskify(seq, mask):
    '''
    Apply a mask to a DNA sequences.

    From "ATGCGT" & "110011" to "ATXXGT"
    From "GTG" & "101" to "CXC"
    '''
    if len(seq) != len(mask):
        raise Error(
            "sequence ('%s') and mask ('%s') have not the same length"
            % (seq, mask) )
    return("".join(
        _revcomp(
            [ x if y=='1' else 'x' for x,y in zip(seq, mask) ]
            )
        )
    )

def _hurricane_death_megatron_300(kargs):
    '''
    Count features (defined by `mask`) from a sequence.

    Input args:
        ((sequenceid, sequence), mask) = kargs

    Output:
        (sequence_id, mask, norm_feature_count)

        - sequence_id : the input sequence id
        - mask : the apply maks
        - norm_feature_count : Pandas DataFrame of count

    Example:
        (sid,msk,nfc) = self._hurricane_death_megatron_300( [ ('SEQ1', 'ATGCGTA'), '101' ] )

    '''
    import pandas as pd
    from collections import Counter
    ((sid, seq), mask) = kargs
    ftlist = list()
    for i in range(0, len(seq) - len(mask) + 1):
        masked = _maskify(str(seq[i:i+len(mask)]), mask)
        if not 'N' in masked:
            ftlist.append( masked )
    ftcnt = Counter(ftlist)
    del ftlist
    ftdf = pd.DataFrame(data=ftcnt, index=[sid]).fillna(0).astype(int).T
    del ftcnt
    ftnorm = ftdf
    # ftnorm = ftdf / ftdf.sum() # vraiment pas sûr de mon coup
    # ftnorm = np.log(ftdf.divide(ftdf.sum(axis=0),axis=1)) # CONCOCT-like normalization
    return(sid, mask, ftnorm)

def _contiguous_kmer_pool(kargs):
    '''
    Count features of length `k` from a sequence.

    Input args:
        ((sequenceid, sequence), k) = kargs

    Output:
        (sequence_id, mask, feature_count)

        - sequence_id : the input sequence id
        - mask : the apply maks
        - feature_count : Pandas DataFrame of count

    Example:
        (sid,msk,nfc) = _contiguous_kmer_pool( [ ('SEQ1', 'ATGCGTATG'), 3 ] )

    '''
    import pandas as pd
    from collections import Counter
    ((sid, seq), k) = kargs
    ftlist = list()
    for i in range(0, len(seq) - k + 1):
        masked = "".join(_revcomp(seq[i:i+k]))
        if not 'N' in masked:
            ftlist.append(masked)
    ftcnt = Counter(ftlist)
    del ftlist
    ftdf = pd.DataFrame(data=ftcnt, index=[sid]).fillna(0).astype(int).T
    del ftcnt
    return(sid, k, ftdf)


class MaskedKmerModel(BaseEstimator, TransformerMixin):
    def __init__(self, mask, mode='strict', verbose=0, n_jobs=1):
        '''Extract (masked) k-mer composition from a DNASequenceBank.

        Usual 4-mers are represented by the mask "1111" but you can set
        spaced-seed (eg: "110011"). Mask must validated the regex "^1[01]*?1$".

        Parameters
        ----------
        
        mask: str
            List of mask to apply (eg: "111", "1111", "110011").

        mode: str (default: 'strict')
            Either 'strict' or 'iupac'. 'strict' will skip non-ATGC kmers while
            'iupac' will consider all kmers, even if they contain 'N'.
            Currently, all non ATGC nucleotide will be considered as 'N'

        verbose:  int (default: 0)
            Verbosity level.

        n_jobs: int (default: 1)
            Number of parallel jobs to run for Kmer extraction.

        '''
        import logging
        # parameters
        if not mode in ['iupac', 'strict']:
            raise ValueError("mode must be 'strict' or 'iupac'. Given: '%s'" % mode)

        self.mask = mask
        self.mode = mode
        self.verbose = verbose
        self.n_jobs = n_jobs
        # get functions to use with Pool
        self._hurricane_death_megatron_300 = _hurricane_death_megatron_300
        self._contiguous_kmer_pool = _contiguous_kmer_pool
        self.k = None

    def fit(self, X):
        import re
        if not re.search('^1[01]*?1$', self.mask):
            raise ValueError("Mask '%s' is not valid" % (self.mask))

        if not re.search('0', self.mask):
            self.k = len(self.mask)

        return self

    def transform(self, X):
        import itertools
        from time import sleep
        import pandas as pd
        from multiprocessing import Pool

        if self.verbose >= 1:
            print("[INFO] Extracting k-mers using mask '%s' (1/3)" % self.mask)

        tmpNORM = {}
        p = Pool(self.n_jobs)

        if self.k:  # we have a contiguous kmer, do not use masking (faster)
            if self.verbose >= 3:
                print("[INFO] Using contiguous kmer method")
            z = itertools.product(list(X.items()), (self.k,))
            result = p.map_async(self._contiguous_kmer_pool, z)
        else:  # we have spaced seed, we have to apply maksing
            if self.verbose >= 3:
                print("[INFO] Using spaced kmer method")
            z = itertools.product(list(X.items()), (self.mask,))
            result = p.map_async(self._hurricane_death_megatron_300, z)

        maxjob = result._number_left

        if self.verbose >= 2:
            while not result.ready():
                _print_progressbar(maxjob - result._number_left, maxjob,
                    msg="Characterizing sequences")
                sleep(1)
            _print_progressbar(maxjob - result._number_left, maxjob,
                msg="Characterizing sequences")
            print()
        else:
            result.wait()
        p.terminate()

        #-- aggregate the results
        for (s, m, df) in result.get():
            # print("ft['%s']['%s']=0" % (s, m))
            try:
                tmpNORM[s][m] = df
            except:
                tmpNORM[s]={}
                tmpNORM[s][m] = df
        del result

        # build index first (AA, AT, AG, AC, ...) to speed up joining
        if self.verbose >= 1:
            print("[INFO] Building feature index (2/3)")

        idx = set()
        i = 0
        for s in tmpNORM.keys():
            i += 1
            if self.verbose >= 2:
                _print_progressbar(i, len(tmpNORM), msg=s)
            for m in tmpNORM[s].keys():
                for k in tmpNORM[s][m].index:
                    idx.add(k)

        if self.verbose >= 2:
            print()

        #-- merge all compositional feature count into one dataframe
        if self.verbose >= 1:
            print("[INFO] Joining features (3/3)")

        i=0
        lstdt=[]
        for s in tmpNORM.keys():
            i+=1
            if self.verbose >= 2:
                _print_progressbar(i, len(tmpNORM.keys()), msg=s)
            tmpX = pd.concat( tmpNORM[s].values() )
            lstdt.append(tmpX)

        if self.verbose >= 2:
            print()

        X = pd.DataFrame(index=idx).join(lstdt)
        del lstdt
        del tmpNORM # clear old data
        X = X.T.fillna(0)
        return(X)


#-------------------------------------------------------------------------------


def get_ind_profiles(seq, K=15):
    '''
    doi: 10.1016/j.physa.2017.04.064

    2.1. Inter-nucleotide distance in a DNA sequence
    '''
    import itertools
    import numpy as np

    # convert 'ATGC' sequences to '0123' sequences
    n2i = { "A":0, "T":1, "G":2, "C":3 }
    seqi = []
    for s in seq:
        try:
            seqi.append(n2i[s])
        except:
            seqi.append(None)

    lastPos = np.zeros(len(n2i), dtype=int)  # last position of nucleotide
    save = np.empty((len(n2i), len(n2i)), object)  # should I calculate the distance?
    count = np.zeros((len(n2i), len(n2i)), dtype=bool)  # distance between 2 nucl
    cpt = 0

    for c in seqi:
        if c is None:  # if not ATGC (mostly N)
            cpt += 1
            continue
        for e in n2i.values():
            if c is None or e is None:
                continue
            if count[e][c]:
                if save[e][c] is None:
                    save[e][c] = []
                save[e][c].append( cpt - int(lastPos[e]) )
                count[e][c] = False
            count[c][e] = True
        lastPos[c] = cpt
        cpt += 1

    f0 = np.zeros((len(n2i), len(n2i), K), float)

    for a,b,k in itertools.product(n2i, n2i, range(K)):
        if save[n2i[a]][n2i[b]] is None:
            f0[ n2i[a] ][ n2i[b] ][k] = 0
        else:
            f0[ n2i[a] ][ n2i[b] ][k] = \
                save[n2i[a]][n2i[b]].count(k+1) / len(save[n2i[a]][n2i[b]]) 

    return f0.reshape(f0.size)


def get_nearest_dissimilar_distance(seq, K=15):
    '''
    doi: 10.1016/j.physa.2017.04.064

    2.1. Inter-nucleotide distance in a DNA sequence
    '''
    import itertools
    import numpy as np

    nucl = ("A", "T", "G", "C")
    W = dict.fromkeys(nucl)
    i = 0
    l = 1

    while i in range(len(seq) - 1):
        if not seq[i] in nucl:
            i+=1
            continue
        if seq[i] != seq[i+1]:
            try:
                W[ seq[i] ].append(l)
            except:
                W[ seq[i] ] = []
                W[ seq[i] ].append(l)
            l = 1
        else:
            l += 1
        i += 1

    if seq[i] in nucl:
        W[seq[i]].append(1)  # count last nucleotide as 1
    
    f1 = np.empty((len(nucl), K), float)

    for a, k in itertools.product(nucl, range(K)):
        f1[ nucl.index(a) ][k] = W[a].count(k+1) / len(W[a])

    return f1.reshape(f1.size)


def inter_nucleotide_distance_profile(kargs):
    (sid, seq), K = kargs
    import numpy as np
    f0 = get_ind_profiles(seq, K=K)
    f1 = get_nearest_dissimilar_distance(seq, K=K)
    X = np.concatenate( [f1, f0] )
    return (sid, X)



class InterNucleotideDistanceModel(BaseEstimator, TransformerMixin):
    '''
    Compute the inter-nucleotide distance profiles as described by Xie et al. (2017)
    '''
    def __init__(self, K=15, verbose=0, n_jobs=1):
        '''
        Compute the inter-nucleotide distance profiles as described by
        Xie et al. (2017) DOI: 10.1016/j.physa.2017.04.064

        Parameters
        ----------
            
        K: int (default: 15)
            Number of distance to consider.

        verbose:  int
            Verbosity level.

        n_jobs: int
            Number of parallel jobs to run for IND profile extraction.

        '''
        self.K = K
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.inter_nucleotide_distance_profile = inter_nucleotide_distance_profile

    def fit(self, X):
        return self

    def transform(self, X):
        import itertools
        from time import sleep
        import pandas as pd
        from multiprocessing import Pool
        if self.verbose >= 1:
            print("[INFO] Extracting IND profiles")
        p = Pool(self.n_jobs)
        result = p.map_async(self.inter_nucleotide_distance_profile,
            zip (list(X.items()), itertools.cycle( (self.K,) )))
        maxjob = result._number_left

        if self.verbose >= 2:
            while not result.ready():
                _print_progressbar(maxjob - result._number_left, maxjob,
                    msg="Characterizing sequences")
                sleep(1)
            _print_progressbar(maxjob - result._number_left, maxjob,
                msg="Characterizing sequences")
            print()
        else:
            result.wait()

        p.terminate()
        # tmpDict = {}
        # for sid, vec in result.get():
        #     tmpDict[sid] = vec
        tmpDict = { k:v for k,v in result.get() }
        return pd.DataFrame(data=tmpDict).T


#-------------------------------------------------------------------------------

class CodingDensityModel(BaseEstimator, TransformerMixin):
    '''
    Compute coding density.
    '''
    def __init__(self, tools=['metageneannotator', 'prodigal', 'fraggenescan'],
                 tmp_fasta="tmp.fennec.fna", force=False, verbose=0, n_jobs=1):
        '''
        Compute coding density.

        Coding density of a sequence is defined as:
            "number of nucleotide included in CDS / total number of nucleotide"

        Software must be located it the PATH. Their respective executables are
        'prodigal', 'run_FragGeneScan.pl' and 'mga'. Both Prodigal and
        FragGeneScane can be installed using Bioconda.


        Parameters
        ----------
        
        tools: list[str] (default: ['metageneannotator', 'prodigal', 'fraggenescan'])
            List of gene prediction tools to use.

        tmp_fasta: str (default: "tmp.fennec.fna")
            File where the DNASequenceBank will be stored (as Fasta)

        force: bool (default: False)
            Force the gene prediciton, even if the resultat are already available

        verbose:  int (default: 0)
            Verbosity level.

        n_jobs: int (default: 1)
            Number of parallel jobs to run for some 
        '''
        import os
        from shutil import which
        from ._utils import run_prodigal, run_metageneannotator, run_fraggenescan

        os.environ["PATH"] = os.environ["PATH"] + ":" + os.path.realpath("bin")

        self.force = force
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.tmp_fasta = tmp_fasta
        self._execpaths = {
            # 'tool identifier', : 'executable path'
            'metageneannotator': which('mga'),
            'prodigal':          which('prodigal'),
            'fraggenescan':      which('run_FragGeneScan.pl')
        }
        self._execfunctions = {
            'metageneannotator': run_metageneannotator,
            'prodigal':          run_prodigal,
            'fraggenescan':      run_fraggenescan
        }
        self._outputs = {
            'metageneannotator': None,
            'prodigal': None,
            'fraggenescan':  None
        }
        self.tools = self._check_tools(tools)

    def _check_tools(self, tools):
        import warnings
        validated = set()
        for t in tools:
            if self._execpaths[t]:
                validated.add(t)
            else:
                warnings.warn("Cannot find '%s'. Will be ignored." % t)
        return validated
        
    def fit(self, X): 
        '''
        Run gene prediction software

        X: DNASequenceBank
        '''
        X._save_sequences(self.tmp_fasta)

        for t in self.tools:
            (retcode, outputfile) = self._execfunctions[t](
                    self._execpaths[t], self.tmp_fasta, force=self.force,
                    verbose=self.verbose, n_jobs=self.n_jobs
                )
            self._outputs[t] = outputfile

        # todo: check if output files exist
        return self

    def transform(self, X=None):
        '''
        Compute coding density

        X: DNASequenceBank, unused
        '''
        import subprocess
        import pandas as pd
        from BCBio import GFF

        tmpX = []
        for t in self.tools:
            if self.verbose >= 1:
                print("[INFO] Computing coding density from %s" % t)

            if self.verbose >= 2:
                # get number of sequences in GFF file
                cmd = 'cat ' + self._outputs[t] + ' | grep -v "^#" | cut -f1 | sort -u | wc -l'
                nbentry = int(subprocess.Popen(cmd, shell=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0])
                i = 0

            ftdensity = {}
            # compute coding density
            for seq in GFF.parse(self._outputs[t]):  # use GFF3? (https://pypi.python.org/pypi/gff3)
                if self.verbose >= 2:
                    i += 1
                    _print_progressbar(i, nbentry, msg=seq.id)
                l = len(seq)  # sequence length
                cl = 0  # number of coding nucleotides
                for f in seq.features:
                    cl += 1 + abs(f.location.start - f.location.end)
                ftdensity[seq.id] = cl/l

            if self.verbose >= 2:
                print()

            tmpX.append(pd.DataFrame(ftdensity, index=["CD_" + t]).T)
            del ftdensity
        ret = pd.DataFrame().join(tmpX, how='outer').fillna(0)
        return ret


#-------------------------------------------------------------------------------


def _get_kmers(seq, k):
    # kmers_list = []
    for i in range(0, len(seq) - k + 1):
        km = seq[i:i + k]
        if not 'N' in km:
            # kmers_list.append(km)
            yield km
    # return kmers_list

def _get_sentence(kargs):
    from ._sentence2vec.sentence2vec import Word, Sentence
    ((sid, seq), k, model) = kargs
    kmers = _get_kmers(seq, k)
    vectors = [model(x) for x in kmers]
    words = [Word(None, v) for v in vectors]
    s = Sentence(words)
    return (sid, s)


class Contig2VecModel(BaseEstimator, TransformerMixin):
    def __init__(self,
                 k=4,
                 modelfile='urq',
                 verbose=0,
                 n_jobs=1
                ):
        '''
        Extract k-mers from sequences, apply a pretrained Dna2Vec model then model
        sequences using sentence2vec.

        NOTE: Dna2vec produces "legacy" .w2v file. You may want to convert it
        using:
            from gensim.models.word2vec import Word2Vec
            model = Word2Vec.load_word2vec_format("dna2vec.output.legacy-w2v")
            model.save("dna2vec.output.w2v")


        References
        ----------
         - https://github.com/pnpnpn/dna2vec
         - https://github.com/peter3125/sentence2vec

        Parameter
        ---------

        k: int (defaut: 4)
            k-mer size

        modelfile: str (default: "urq")
            Either a pretrained Dna2Vec model label (see self.available_models)
            or an model file path.

        verbose:  int (default: 0)
            Verbosity level.

        n_jobs: int (default: 1)
            Number of jobs.

        Example
        -------

        W = Contig2VecModel().fit_transform(...)
        kpca = KernelPCA(kernel='cosine')
        W_kpca = kpca.fit_transform(W)
        tsne = TSNE(n_components=2, verbose=2, n_iter=1500, n_jobs=160)
        W_tsne = tsne.fit_transform(W_kpca)
        plt.scatter(W_tsne[:,0], W_tsne[:,1], s=5)
        plt.show()

        '''
        from glob import glob
        from os.path import dirname
        self.k = k
        self.available_models = {
                f.split("-")[-1].split('.')[0]: f
                    for f in glob(dirname(__file__)+"/pretrained/*w2v")
            }
        self.modelfile = modelfile
        if self.modelfile in self.available_models.keys():
            self.modelfile = self.available_models[self.modelfile]
        self.model = None  # word2vec model will be loaded by .fit()
        self.verbose = verbose
        self.n_jobs = n_jobs
        self._get_sentence = _get_sentence
        # self._word_vec = self.model.wv.word_vec  # shorcut to get vector, 7x faster than `self.model[word]`

    def fit(self, X):
        '''
        If dna2vec model is not available, train it. Otherwise, just load it.

        Parameter
        ---------

        X: DNASequenceBank
        '''
        from gensim.models import word2vec
        if self.modelfile:
            # load the pretrained model
            if self.verbose >= 1:
                print("[INFO] loading model from '%s'" %(self.modelfile))
            self.model = word2vec.Word2Vec.load(self.modelfile)
        else:
            # not supported yet + it's super long to run!
            # but we could train on `X` only
            pass
        return self


    def transform(self, X):
        '''
        Parameter
        ---------
        X: DNASequenceBank
        '''
        import pandas as pd
        from ._sentence2vec.sentence2vec import Word, Sentence, sentence_to_vec

        if not self.model:
            raise Error("[ERROR] Model has not been loaded")

        sentences = {}
        step = 0

        if self.verbose >= 1:
            print("[INFO] Extracting sentences (1/2)")

        for sid, seq in X.items():
            if self.verbose >= 2:
                step += 1
                _print_progressbar(step, len(X), msg=sid)
            sentences[sid] = Sentence([Word(x, self.model[x])
                for x in _get_kmers(seq, self.k)])

        if self.verbose >= 2:
            print()

        # convert Sentences to vectors
        if self.verbose >= 1:
            print("[INFO] Converting sentences to vectors (2/2)")

        seqids = list(sentences.keys())
        vectors = sentence_to_vec(list(sentences.values()), self.model.vector_size)
        d = dict(zip(seqids, vectors))
        del sentences, seqids, vectors

        return pd.DataFrame(data=d).T


    def transform_par(self, X):
        '''

        Parameter
        ---------

        X: DNASequenceBank
        '''
        import itertools
        import pandas as pd
        from multiprocessing import Pool
        from time import sleep
        from ._sentence2vec.sentence2vec import sentence_to_vec

        if not self.model:
            raise Error("Model has not been loaded")

        if self.verbose >= 1:
            print("[INFO] Extracting sentences (1/2)")

        p = Pool(self.n_jobs)
        z = zip(X.items(),
                itertools.repeat(self.k),
                itertools.repeat(self.model.wv.word_vec)
               )
        result = p.map_async(_get_sentence, z)

        maxjob = result._number_left

        if self.verbose >= 2:
            while not result.ready():
                _print_progressbar(maxjob - result._number_left, maxjob,
                    msg="Characterizing sequences")
                sleep(1)
            _print_progressbar(maxjob - result._number_left, maxjob,
                msg="Characterizing sequences")
            print()
        else:
            result.wait()
        p.terminate()

        sentences = {k:v for k,v in result.get()}
        del result

        if self.verbose >= 2:
            print()

        # convert Sentence to vector
        if self.verbose >= 1:
            print("[INFO] Converting sentences to vectors (2/2)")

        seqids = list(sentences.keys())
        vectors = sentence_to_vec(list(sentences.values()), self.model.vector_size)
        d = dict(zip(seqids, vectors))
        del sentences, seqids, vectors

        return pd.DataFrame(data=d).T


#-------------------------------------------------------------------------------

class SequenceAbundanceModel(BaseEstimator, TransformerMixin):
    def __init__(self,
            abdfiles,
            verbose=0,
            n_jobs=1,
            prefix="cov_"
        ):
        self.abdfiles = { x: self._get_file_format(x) for x in abdfiles }
        self.verbose = verbose
        self.n_jobs = n_jobs
        self._available_format = ['CSV', 'TSV']  #, 'BAM', 'SAM', 'FASTA']
        self.prefix = prefix
        '''
        Return per-sequence abundance from `abdfiles`.

        Parameter
        ---------

        abdfiles: list[string]
            Files contains average coverage

        verbose:  int (default: 0)
            Verbosity level.

        n_jobs: int (default: 1)
            Ingored.
        '''

    def _get_file_format(self, file):
        return file.split(".")[-1].upper()
        # if ext not in self._available_format:
        #     raise Error("File format '%s' is not supported")
        # return ext

    def fit(self, X):
        '''

        Parameter
        ---------

        X: DNASequenceBank
        '''
        import pandas as pd

        self.data = pd.DataFrame(index=list(X.keys()))
        i = 0
        for file, ext in self.abdfiles.items():
            if ext not in self._available_format:
                print("Skipping '%s': format '%s' is not supported" % (file, ext))
                continue
            if ext == "CSV":
                tmp = pd.read_csv(file, index_col=0)
                tmp.columns = (self.prefix + str(i))
                self.data.join(tmp, how='outer')
            elif ext == "TSV":
                tmp = pd.read_csv(file, index_col=0, sep="\t", columns=None)
                tmp.columns = (self.prefix+str(i))
                self.data = self.data.join(tmp, how='outer')
            else:
                print("Do not know what to do...")  # this should not happen
            # elif ext == "SAM":
            #     self.data = self._get_coverage_from_sam(file)
            # elif ext == "BAM":
            #     self.data = self._get_coverage_from_bam(file)
            # elif ext == "FASTA":
            #     self.data = self._map_fasta(file, self.ref_contigs)
            i += 1

        return self

    def transform(self, X):
        '''

        Parameter
        ---------

        X: DNASequenceBank
        '''
        return self.data
