
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

    import timeg    n=32
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

    def _save_sequences(self, outfile):
        '''Save sequences to a FASTA file.
        One of the valid extensions will be added to `outfile` i not found.
        If `outfile` already exists, it will be overrided.
        '''
        import re

        extensions = [ 'fasta', 'fa', 'fna' ]
        if not outfile.split('.')[-1] in extensions:
            outfile += '.fna'

        with open(outfile, "w") as f:
            for sid, seq in self.items():
                _ = f.write(">" + sid + "\n")
                _ = f.write(re.sub(r'(.{,80})',r"\1\n", seq))

