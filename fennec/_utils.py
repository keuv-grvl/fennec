
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
    '''Store Fasta file as dict {sequence_ID: sequence}'''
    def __init__(self, fastafile, min_length=1000, mode='strict', nb_seq=0, verbose=0):
        '''
        Parameters
        ----------

        fastafile: str
            Path to a Fasta file

        min_length: int (default: 1000)
            Minimum sequence length, shorter sequences will be skipped

        mode: str (default: 'strict')
            Either 'strict' or 'iupac'. 'strict' will replace non-ATGC nucleotides
            with Ns while 'iupac' will replace non-ATGCYRSWKMBDHVN by Ns

        nb_seq: int (default: 0)
            Load the `nb_seq` first sequences, all if 0

        verbose: int (default: 0)
            Verbosity level
        '''
        import re
        self._patterns = {
            "strict": re.compile(r"[^ATGC]", flags=re.IGNORECASE),
            "iupac": re.compile(r"[^ATGCYRSWKMBDHVN]", flags=re.IGNORECASE)
        }

        if not mode in self._patterns.keys():
            raise ValueError("mode must be 'strict' or 'iupac'. Given: '%s'" % mode)
        
        from skbio.io import read as FastaReader
        self.fastafile = fastafile
        self.min_length = min_length
        self.mode = mode
        self.verbose = verbose

        if nb_seq > 0:
            self.nb_seq = nb_seq
        else:
            self.nb_seq = sum(
                    1 for x in FastaReader(self.fastafile, format='fasta', verify=False)
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
            if (len(s) <= self.min_length):
                continue

            i += 1
            sid = s.metadata['id']
            
            if self.verbose >= 2:
                _print_progressbar(i, self.nb_seq, msg=sid)
    
            self[sid] = self._patterns[self.mode].sub("N", str(s))
            
            if i >= self.nb_seq:
                break

        if self.verbose >= 2:  print()
        if self.verbose >= 1:
            print("[INFO] %d sequences loaded" % i)

    def _save_sequences(self, outfile, ids=None, append=False):
        '''Save sequences to a FASTA file.
        One of the valid extensions will be added to `outfile` if not found.
        If `outfile` already exists, it will be overrided without warning.
        '''
        import re
        
        opt_open = "w"
        if append:
            opt_open = "a"

        extensions = [ 'fasta', 'fa', 'fna' ]
        if not outfile.split('.')[-1] in extensions:
            outfile += '.fna'

        if self.verbose >= 1:
            print("[INFO] Writing sequences to '%s'" % outfile)

        with open(outfile, opt_open) as f:
            i = 0
            if ids is None:
                for sid, seq in self.items():
                    i += 1
                    if self.verbose >= 2:
                        _print_progressbar(i, len(self), msg=sid)
                    _ = f.write(">" + sid + "\n")
                    _ = f.write(re.sub(r'(.{,80})',r"\1\n", seq))
            else:
                for sid, seq in self.items():
                    if sid not in ids:
                        next
                    i += 1
                    if self.verbose >= 2:
                        _print_progressbar(i, len(ids), msg=sid)
                    _ = f.write(">" + sid + "\n")
                    _ = f.write(re.sub(r'(.{,80})',r"\1\n", seq))


        if self.verbose >= 2:  print()
        if self.verbose >= 1:
            print("[INFO] %d sequences written to '%s'" % (i, outfile))
