
import logging

from sklearn.base import BaseEstimator, TransformerMixin

from ._utils import _print_progressbar


#-------------------------------------------------------------------------------

# Functions called by multiprocessing.Pool().map must be declared out of the class
# See https://stackoverflow.com/questions/32321324/pool-within-a-class-in-python
_BASE_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", "x": "x"}

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
    Count features (defined by `mask`) from a sequence and normalize feature count.

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


class MaskedKmerModel(BaseEstimator, TransformerMixin):
    """
    kc = MaskedKmerModel(["1111"", "110011"], verbose=2, n_jobs=32)
    X = fit_transform( DNASequencebank(...))
    """
    def __init__(self, masks, verbose=0, n_jobs=1):
        '''Extract k-mer composition from a Fasta file to model a DNASequenceBank.

        Usual 4-mers are represented by the mask "1111" but you can set
        spaced-seed (eg: "110011").

        Parameters
        ----------
        
        masks: list[str]
            List of mask to apply (eg: ["111", "1111"]).

        verbose:  int
            Verbosity level.

        n_jobs: int
            Number of parallel jobs to run for Kmer extraction.

        '''
        import logging
        self._logger = logging.getLogger(__name__)
        # paramters
        self.masks = masks
        self.verbose = verbose
        self.n_jobs = n_jobs
        # get functions to use with Pool
        self._hurricane_death_megatron_300 = _hurricane_death_megatron_300

        self._check_masks()

    def _check_masks(self):
        '''
        Check if masks are valid
        '''
        import re
        for m in self.masks:
            if re.search('[^01]', m):
                raise ValueError(
                    "Mask '%s' in [ %s ] is not valid" 
                    % (m, ', '.join(self.masks)))

    def fit(self, X):
        return self

    def transform(self, X):
        import itertools
        import time
        import pandas as pd
        from multiprocessing import Pool

        if self.verbose >= 1:
            print("[INFO] Extracting masked k-mers")

        tmpNORM = {}
        p = Pool(self.n_jobs)
        z = itertools.product( list(X.items()), self.masks)
        result = p.map_async(self._hurricane_death_megatron_300, z)
        maxjob = result._number_left

        if self.verbose >= 2:
            while not result.ready():
                _print_progressbar(maxjob - result._number_left, maxjob,
                    msg="Characterizing sequences")
                time.sleep(0.5)
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
            print("[INFO] Building feature index")

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
            print("[INFO] Joining features")

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
        X = X.T
        X = X.fillna(0)
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
    DOI: 10.1016/j.physa.2017.04.064
    '''
    def __init__(self, K=15, verbose=0, n_jobs=1):
        self.K = K
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.inter_nucleotide_distance_profile = inter_nucleotide_distance_profile

    def fit(self, X):
        return self

    def transform(self, X):
        import itertools
        import time
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
                time.sleep(0.5)
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

