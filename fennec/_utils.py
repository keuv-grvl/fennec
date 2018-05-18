
def _max_proc(min_cpu=1):
    '''
    Return the maximum number of available CPU defined as:
        number of cpu - roundup(average load over the last minute)

    Will return at least `min_cpu` (1 by default)
    '''
    import os
    return max(os.cpu_count() - int(os.getloadavg()[0] + 1), min_cpu)


def _print_progressbar(step, maxi, msg="", char="=", width=50):
    '''
    Print a progress bar then place the cursor at the begging of the line.
    Display can be really messy if `maxi` is set incorrectly.

    import time
    n=32
    for i in range(n):
        time.sleep(0.1)
        _print_progressbar(i+1, n, msg="Test", char='=', width=50)
    print()

    '''
    # rows, columns = os.popen('stty size', 'r').read().split()
    p = int(100 * step / maxi)
    print("[%s>%s] %d%% (%d/%d) %-20s" %
        (char * int(p * width/100),
        (" " * (width-int(p*width/100))),
        p, step, maxi, msg),
        end="\r", flush=True)


def run_prodigal(
        execpath,
        inputfile,
        outputfile="tmp.prodigal.gff",
        force=False,
        verbose=0,
        n_jobs=1,
    ):
    '''
    Predict genes from `inputfile` (FASTA format) using Prodigal.

    Returns the Prodigal return code, or -1 if `outputfile` already exists.


    Parameters
    ----------

    execpath: str
        Path of the Prodigal executable.

    inputfile: str
        Fasta file to predict genes from.

    outputfile: str (default: "tmp.prodigal.gff")
        Name of the GFF output file.

    force: bool (default: False)
        Choose to run Prodigal or not if `outputfile` already exists.

    verbose: int (default: 0)
        Verbosity level.

    n_jobs: int (default: 1)
        Ignored.

    '''
    import os
    import subprocess
    import sys
    from skbio.io import read as FastaReader

    retcode = -1
    if force or not os.path.isfile(outputfile):
        nb_seq = sum(1 for x in FastaReader(inputfile, format='fasta', verify=True))
        i = 0
        cmd = [execpath, "-q", "-f gff", "-i", inputfile]
        if verbose >= 1:
            print("[INFO] Predicting genes with Prodigal")
            print("[INFO] Running '%s'" % (" ".join(cmd)))
        with open(outputfile, "w") as outfile:
            p = subprocess.Popen(" ".join(cmd), shell=True,
                stdout=subprocess.PIPE) #, stderr= subprocess.DEVNULL)
            seqid="NULL"
            for x in p.stdout:
                xx = x.decode(sys.getdefaultencoding()).rstrip()
                print(xx, file=outfile)
                if xx.startswith("# Sequence Data:"):
                    seqid=xx.split("\"")[1]
                    i += 1
                    if verbose >= 2:
                        _print_progressbar(i, nb_seq, msg=seqid)
            p.wait()
            p.terminate()
            if verbose >= 2:
                print()
            retcode = p.returncode

    return (retcode, outputfile)



def run_fraggenescan(
        execpath,
        inputfile,
        outputfile="tmp.fraggenescan.gff",
        force=False,
        verbose=0,
        n_jobs=1,
    ):
    '''
    Predict genes from `inputfile` (FASTA format) using FragGeneScan.

    Returns the FragGeneScan return code, or -1 if `outputfile` already exists.


    Parameters
    ----------

    execpath: str
        Path of the 'run_FragGeneScan.pl' script.

    inputfile: str
        Fasta file to predict genes from.

    outputfile: str (default: "tmp.fraggenescan.gff")
        Name of the GFF output file.

    force: bool (default: False)
        Choose to run FragGeneScan or not if `outputfile` already exists.

    verbose: int (default: 0)
        Verbosity level.

    n_jobs: int (default: 1)
        Number of CPU FragGeneScan will use.

    '''
    import os
    import subprocess

    retcode = -1
    outputlabel = os.path.splitext(outputfile)[0]
    if force or not os.path.isfile(outputfile):
        cmd = [execpath,
            '-genome='+str(inputfile),
            '-out='+str(outputlabel),
            '-thread='+str(n_jobs),
            '-complete=1', '-train=complete' ]
        if verbose >= 1:
            print("[INFO] Predicting genes with FragGeneScan")
            print("[INFO] Running '%s'" % (" ".join(cmd)))
        res = subprocess.run(cmd)
        retcode = res.returncode

    return (retcode, outputfile)



def run_metageneannotator(
        execpath,
        inputfile,
        outputfile="tmp.metageneannotator.gff",
        force=False,
        verbose=0,
        n_jobs=1,
    ):
    '''
    Predict genes from `inputfile` (FASTA format) using MetaGeneAnnotator.

    MetaGeneAnnotator does not output a GFF3 file. Thus the output is parsed
    and formatted to comply with the GFF3 standard.

    Returns the MetaGeneAnnotator return code, or -1 if `outputfile` already exists.


    Parameters
    ----------

    execpath: str
        Path of the MetaGeneAnnotator executable.

    inputfile: str
        Fasta file to predict genes from.

    outputfile: str (default: "tmp.metageneannotator.gff")
        Name of the GFF output file.

    force: bool (default: False)
        Choose to run MetaGeneAnnotator or not if `outputfile` already exists.

    verbose: int (default: 0)
        Verbosity level.

    n_jobs: int (default: 1)
        Ignored.
    '''
    import os
    import subprocess
    import sys
    from skbio.io import read as FastaReader

    retcode = -1
    if force or not os.path.isfile(outputfile):
        if verbose >= 1:
            print("[INFO] Predicting genes with MetaGeneAnnotator")
        cmd = [execpath, '-m', inputfile]
        if verbose >= 1:
            print("[INFO] Running '%s'" % (" ".join(cmd)))
        nb_seq = sum(1 for x in FastaReader(inputfile, format='fasta', verify=True))
        i = 0
        with open(outputfile, "w") as outfile:
            p = subprocess.Popen(" ".join(cmd), shell=True,
                stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            print("##gff-version 3", file=outfile)  # GFF3 header
            seqid = "null"
            for x in p.stdout:
                xx = x.decode(sys.getdefaultencoding()).rstrip()
                if xx.startswith("#"):
                    if not xx.startswith("# gc") and not xx.startswith("# self"):
                        seqid = xx[2:]
                        i += 1
                        if verbose >= 2:
                            _print_progressbar(i, nb_seq, msg=seqid)
                else:
                    (geneid, start, end, strand, frame, _, score,
                        _, _, _, _) =  xx.split("\t")
                    print(seqid, "MGA", "gene", start, end, score, strand, frame, "-", sep="\t", file=outfile)
            p.wait()
            p.terminate()
            if verbose >= 2:
                print()
            retcode = p.returncode

    return (retcode, outputfile)



#-------------------------------------------------------------------------------

class DNASequenceBank(dict):
    '''
    Load a FASTA file, chunks long sequences and store chunks as dict
    {chunk_ID: chunk} while tracking must-link relationships between chunks.
    '''
    def __init__(self, min_length=1000, chunk_size=10000, overlap="auto",
                    max_nb_seq=0, mode='strict', verbose=0):
        '''
        Parameters
        ----------

        min_length: int (default: 1000)
            Minimum sequence length, shorter sequences will be skipped

        chunk_size: int (default: 10000)
            DNA chunk size. If 0, sequences are not splitted.

        overlap: int (default: 'auto')
            Size of overlap between two chunks.
            If "auto", overlap is set to `chunk_size // 7`

        max_nb_seq: int (default: 0)
            Load the `max_nb_seq` first sequences, all if 0. This is NOT the
            final number of chunks.

        mode: str (default: 'strict')
            Either 'strict' or 'iupac'. 'strict' will replace non-ATGC nucleotides
            with Ns while 'iupac' will replace non-ATGCYRSWKMBDHVN by Ns

        verbose: int (default: 0)
            Verbosity level
        '''
        import re
        self._patterns = {
            "strict": re.compile(r"[^ATGC]", flags=re.IGNORECASE),
            "iupac": re.compile(r"[^ATGCYRSWKMBDHVN]", flags=re.IGNORECASE)
        }

        if not mode in self._patterns.keys():
            raise ValueError(f"mode must be 'strict' or 'iupac'. Given: {mode}")
        
        from skbio.io import read as FastaReader
        self.min_length = min_length
        self.max_nb_seq = max_nb_seq
        self.mode = mode
        self.verbose = verbose
        self.chunk_size = chunk_size
        self.overlap = overlap
        if self.overlap == "auto":
            self.overlap_ = self.chunk_size // 7  # this is completely arbitrary
        else:
            self.overlap_ = self.overlap

    def __repr__(self):
        '''Display object in a sklearn-style'''
        fields = ['fastafile', 'min_length', 'max_nb_seq', 'mode', 'verbose', 'chunk_size', 'overlap', 'overlap_', 'nb_chunks']
        return "%s(%s)" % (self.__class__.__name__, ', '.join([f+"="+self.__dict__[f].__repr__() for f in fields if f in self.__dict__]))

    def _chunks(self, s, n, o, merge_last):
        '''
        Yield successive `n`-sized chunks from DNA sequence `s` with given
        overlap `o` between the chunks. If `n` is 0, do not split `s`.

        Highly inspired from: https://github.com/BinPro/CONCOCT/blob/master/scripts/cut_up_fasta.py
        '''
        assert n >= 0, f"Chunk size must be 0 or more, got {n}."
        if n >= len(s) or n == 0:
            yield s
        else:
            if not n > o:
                raise ValueError(f"Chunk size ({n}) must be higher than the overlap ({o}).")
            if not merge_last:
                for i in range(0, len(s), n - o):
                    yield s[i:i + n]
            else:
                for i in range(0, len(s) - n + 1, n - o):
                    yield s[i:i + n] if i + n + n - o <= len(s) else s[i:]

    def read_fasta(self, path, format='fasta'):
        '''
        Load sequences from a FASTA file to a dict.

        Parameters
        ----------

        path: str
            Path to a sequence file

        format: str (default: 'fasta')
            The format of the file if known. If `None`, the format will be
            inferred from the file.

        '''
        import os
        import numpy as np
        from skbio.io import read as FastaReader
        from scipy.sparse import dok_matrix
        from itertools import product
        from functools import reduce

        assert os.path.isfile(path), f"'{path}' is not a file"
        
        self.fastafile = path

        if self.max_nb_seq <= 0:
            self.max_nb_seq = sum(
                [1 for x in FastaReader(self.fastafile, format=format, verify=False) if len(x) >= self.min_length]
            )

        if self.verbose >= 1:
            print("[INFO] Reading input")

        i = 0
        mustlink = []
        for s in FastaReader(self.fastafile, format=format, verify=True):
            if (len(s) < self.min_length):
                continue

            i += 1
            # merge last chunck if shorter than `min_length`
            merge_last = len(s) % self.chunk_size < self.min_length
            tmpML = set()

            for j, split_seq in enumerate(self._chunks(str(s), self.chunk_size, self.overlap_, merge_last)):
                fid = f"{s.metadata['id']}__{j}"
                tmpML.add(fid)
                self[fid] = self._patterns[self.mode].sub("N", str(split_seq))

            mustlink.append(tmpML)

            if self.verbose >= 2:
                _print_progressbar(i, self.max_nb_seq, msg=s.metadata['id'])
            
            if i >= self.max_nb_seq:
                break

        if self.verbose >= 2:
            print()

        if self.verbose >= 1:
            print("[INFO] Formatting must-link matrix")

        self.mustlink_matrix = dok_matrix((len(self), len(self)), dtype=np.bool)
        self._seqid2numid = {k:i for i,k in enumerate(self.keys())}

        for i, mlset in enumerate(mustlink):
            if self.verbose >= 2:
                _print_progressbar(i+1, len(mustlink))
            for id1, id2 in product(mlset, mlset):
                self.mustlink_matrix[self._seqid2numid[id1], self._seqid2numid[id2]] = True

        if self.verbose >= 2:
            print()

        del mustlink ; del tmpML
        self.nb_chunks = len(self)
        if self.verbose >= 1:
            print("[INFO] %d sequence chuncks loaded" % self.nb_chunks)
            print("[INFO] %d must-link relationships found (%.4f%%)" % (
                    self.mustlink_matrix.nnz,
                    100 * self.mustlink_matrix.nnz / reduce(lambda x, y: x*y, self.mustlink_matrix.shape)
               ))

    def to_fasta(self, ids=None, compress=False, append=False):
        '''
        Export splitted sequences to a FASTA file.
        One of the valid extensions will be added to `path` if not found.
        If `path` already exists, it will be overrided without warning.

        Note: the output file name is based on the input file name. For example,
        sequences from `/path/to/file.fasta` loaded using DNASequenceBank using
        default parameters will be saved in `/path/to/file.l1000c10000o1428.fasta`
        denoting `min_length`, `chunk_size` and `overlap`.

        Parameters
        ----------

        ids: list of str (default: None)
            List of IDs to export.

        compress: bool (default: False)
            Compress the output file using gzip.

        append: bool (default: False)
            Append sequences if the output file already exists.
        '''
        import re
        import gzip

        if not self.fastafile:
            raise Exception("No file loaded")

        tmpnameparts = self.fastafile.split('.')
        tmpnameparts.append(tmpnameparts[-1])
        tmpnameparts[-2] = f"l{self.min_length}c{self.chunk_size}o{self.overlap_}"
        path = ".".join(tmpnameparts)
        if compress:
            path += ".gz"
        
        if self.verbose >= 1:
            print("[INFO] Writing sequences to '%s'" % path)

        if append:
            opt_open = 'a'
        else:
            opt_open = 'w'

        if compress:
            f = gzip.open(path, opt_open + 't')  # open in text mode
        else:
            f = open(path, opt_open)

        i = 0
        if ids is None:
            ids = self.keys()

        for sid, seq in self.items():
            if sid not in ids:
                next
            i += 1
            if self.verbose >= 2:
                _print_progressbar(i, len(ids), msg=sid)

            _ = f.write(">" + sid + "\n")
            _ = f.write(re.sub(r'(.{,80})',r"\1\n", seq))  # 80 nt per line

        f.close()

        if self.verbose >= 2:
            print()
        if self.verbose >= 1:
            print(f"[INFO] {i} sequences written to '{path}'")

    @staticmethod
    def _store_sparse_mat(M, filename, name):
        '''
        Store a csr matrix in HDF5 (https://stackoverflow.com/a/44282655)

        Parameters
        ----------
        M : scipy.sparse.csr.csr_matrix
            sparse matrix to be stored

        name: str
            node prefix in HDF5 hierarchy

        filename: str
            HDF5 filename
        '''
        import numpy as np
        import tables
        from scipy.sparse.csr import csr_matrix

        M = csr_matrix(M)
        
        assert(M.__class__ == csr_matrix), "M must be a csr matrix"
        
        with tables.open_file(filename, 'a') as f:
            for attribute in ('data', 'indices', 'indptr', 'shape'):
                full_name = f'{name}_{attribute}'
                # remove existing nodes
                try:
                    n = getattr(f.root, full_name)
                    n._f_remove()
                except AttributeError:
                    pass
                # add nodes
                arr = np.array(getattr(M, attribute))
                atom = tables.Atom.from_dtype(arr.dtype)
                ds = f.create_carray(f.root, full_name, atom, arr.shape)
                ds[:] = arr

    @staticmethod
    def _load_sparse_mat(filename, name):
        '''
        Load a csr matrix from HDF5 (https://stackoverflow.com/a/44282655)

        Parameters
        ----------
        name: str
            node prefix in HDF5 hierarchy

        filename: str
            HDF5 filename

        Returns
        ----------
        M : scipy.sparse.csr.csr_matrix
            loaded sparse matrix
        '''
        import numpy as np
        import tables
        from scipy import sparse

        with tables.open_file(filename) as f:
            # get nodes
            attributes = []
            for attribute in ('data', 'indices', 'indptr', 'shape'):
                attributes.append(getattr(f.root, f'{name}_{attribute}').read())

        # construct sparse matrix
        M = sparse.csr_matrix(tuple(attributes[:3]), shape=attributes[3])
        return M

    def to_hdf(self, outfile):
        import h5py, sys
        try:
            hf = h5py.File(outfile)
            fields = ['fastafile', 'min_length', 'max_nb_seq', 'mode', 'verbose', 'chunk_size', 'overlap', 'overlap_', 'nb_chunks']

            for f in fields:
                hf.attrs[f] = self.__dict__[f]
            #- encode ids and sequences to UTF-8 (https://github.com/h5py/h5py/issues/289) and store them separately
            hf.create_dataset('_data_dnabank/_useqids', data=[s.encode('utf-8') for s in self.keys()])
            hf.create_dataset('_data_dnabank/_useqs', data=[s.encode('utf-8') for s in self.values()])
            hf.close()
            #- store must-link data as a dense matrix
            DNASequenceBank._store_sparse_mat(self.mustlink_matrix, outfile, '_data_mustlink_matrix')
        except:
            print(f"[ERROR] while writing to {outfile}")
            sys.exit(1)

    @staticmethod
    def read_hdf(path):
        '''
        Read HDF5 file to load DNASequenceBank data (squences and must-link matrix)

        Parameter:
        ----------

        path: str
            Path to the HDF5 file
        '''
        import h5py, sys
        try:
            hf = h5py.File(path, 'r')
            #- build full DNASequenceBank object to return
            ret = DNASequenceBank()
            for k, v in hf.attrs.items():
                setattr(ret, k, v)
            if not '_data_dnabank' in hf.keys():
                raise Exception("Cannot find appropriate data group")
            #- rebuild dict from keys and values
            for k, v in zip([x.decode('utf-8') for x in hf.get('_data_dnabank/_useqids')],
                           [x.decode('utf-8') for x in hf.get('_data_dnabank/_useqs')]):
                ret[k] = v
            #- load mustlink_matrix
            setattr(ret, 'mustlink_matrix', DNASequenceBank._load_sparse_mat(path, '_data_mustlink_matrix').todok())
            hf.close()
        except:
            print(f"[ERROR] while writing to {outfile}")
            sys.exit(1)

        return ret
