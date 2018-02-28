
import logging

from sklearn.base import BaseEstimator, TransformerMixin

from ._utils import _print_progressbar

#-------------------------------------------------------------------------------

class DNASequenceBank(dict):
    '''Store Fasta file as dict {sequence_ID: sequence }'''
    def __init__(self, fastafile, min_length, nb_seq=0, verbose=0):
        '''
        Parameters
        ----------

        fastafile: str
            Path to a Fasta file

        min_length: int
            Minimum sequence length, shorter sequences will be skipped

        nb_seq: int (default: 0)
            Load the `nb_seq` first sequences, all if 0

        verbose: int (default: 0)
            Verbosity level
        '''
        from skbio.io import read as FastaReader
        self.fastafile = fastafile
        self.min_length = min_length
        self.verbose = verbose
        if nb_seq > 0:
            self.nb_seq = nb_seq
        else:
            import skbio
            self.nb_seq = sum(
                    1 for x in FastaReader(self.fastafile, format='fasta')
                )
        self._load_sequences()

    def _load_sequences(self):
        '''Load sequences longer than `min_length` from a FASTA file to a dict.
        '''
        from skbio.io import read as FastaReader
        if self.verbose >= 1:
            print("[INFO] Reading input")
        i = 0
        for s in FastaReader(self.fastafile, format='fasta', verify=True):
            if (len(s) < self.min_length):
                next
            i += 1
            sid = s.metadata['id']
            if self.verbose >= 2:
                _print_progressbar(i, self.nb_seq, msg=sid)
            # print("%d/2305" % i)
            self[sid] = str(s)
            if i >= self.nb_seq:
                break
        if self.verbose >=2:  print()
        if self.verbose >= 1:
            print("[INFO] %d sequences loaded" % i)
        # s  = pandas.Series(my_dict,index=my_dict.keys())



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
            print("[INFO] Characterizing sequences")

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
        else:
            result.wait()

        if self.verbose >= 2:
            _print_progressbar(maxjob - result._number_left, maxjob,
                msg="Characterizing sequences")
            print()

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

