def isinteractive():
    import sys

    return bool(getattr(sys, "ps1", sys.flags.interactive))


def _max_proc(min_cpu=1):
    """
    Return the maximum number of available CPU defined as:
        number of cpu - roundup(average load over the last minute)

    Will return at least `min_cpu` (1 by default)
    """
    import os

    return max(os.cpu_count() - int(os.getloadavg()[0] + 1), min_cpu)


def _print_progressbar(step, maxi, msg="", char="=", width=50):
    """
    Print a progress bar then place the cursor at the begging of the line.
    Display can be really messy if `maxi` is set incorrectly.

    import time
    n=32
    for i in range(n):
        time.sleep(0.1)
        _print_progressbar(i+1, n, msg="Test", char='=', width=50)
    print()

    """
    # rows, columns = os.popen('stty size', 'r').read().split()
    p = int(100 * step / maxi)
    print(
        "[%s>%s] %d%% (%d/%d) %-20s"
        % (
            char * int(p * width / 100),
            (" " * (width - int(p * width / 100))),
            p,
            step,
            maxi,
            msg,
        ),
        end="\r",
        flush=True,
    )


def run_prodigal(
    execpath, inputfile, outputfile="tmp.prodigal.gff", force=False, verbose=0, n_jobs=1
):
    """
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

    """
    import os
    import subprocess
    import sys
    from skbio.io import read as FastaReader

    retcode = -1
    if force or not os.path.isfile(outputfile):
        nb_seq = sum(1 for x in FastaReader(inputfile, format="fasta", verify=True))
        i = 0
        cmd = [execpath, "-q", "-f gff", "-i", inputfile]
        if verbose >= 1:
            print("[INFO] Predicting genes with Prodigal")
            print("[INFO] Running '%s'" % (" ".join(cmd)))
        with open(outputfile, "w") as outfile:
            p = subprocess.Popen(
                " ".join(cmd), shell=True, stdout=subprocess.PIPE
            )  # , stderr= subprocess.DEVNULL)
            seqid = "NULL"
            for x in p.stdout:
                xx = x.decode(sys.getdefaultencoding()).rstrip()
                print(xx, file=outfile)
                if xx.startswith("# Sequence Data:"):
                    seqid = xx.split('"')[1]
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
    """
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

    """
    import os
    import subprocess

    retcode = -1
    outputlabel = os.path.splitext(outputfile)[0]
    if force or not os.path.isfile(outputfile):
        cmd = [
            execpath,
            "-genome=" + str(inputfile),
            "-out=" + str(outputlabel),
            "-thread=" + str(n_jobs),
            "-complete=1",
            "-train=complete",
        ]
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
    """
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
    """
    import os
    import subprocess
    import sys
    from skbio.io import read as FastaReader

    retcode = -1
    if force or not os.path.isfile(outputfile):
        if verbose >= 1:
            print("[INFO] Predicting genes with MetaGeneAnnotator")
        cmd = [execpath, "-m", inputfile]
        if verbose >= 1:
            print("[INFO] Running '%s'" % (" ".join(cmd)))
        nb_seq = sum(1 for x in FastaReader(inputfile, format="fasta", verify=True))
        i = 0
        with open(outputfile, "w") as outfile:
            p = subprocess.Popen(
                " ".join(cmd),
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
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
                    (
                        geneid,
                        start,
                        end,
                        strand,
                        frame,
                        _,
                        score,
                        _,
                        _,
                        _,
                        _,
                    ) = xx.split("\t")
                    print(
                        seqid,
                        "MGA",
                        "gene",
                        start,
                        end,
                        score,
                        strand,
                        frame,
                        "-",
                        sep="\t",
                        file=outfile,
                    )
            p.wait()
            p.terminate()
            if verbose >= 2:
                print()
            retcode = p.returncode

    return (retcode, outputfile)


def boxplot(x, g, n, min_nb_seq=25):
    """
    Plot a violon plot using `x` values into classes defined by `g`.
    """
    import numpy as np
    from matplotlib import pyplot as plt
    from collections import Counter

    fig, axes = plt.subplots(figsize=(11, 8))
    axes.axhline(y=0, color="red", linewidth=1)
    categ = np.unique(g)
    nbseq = Counter(g)
    d = []
    for c in categ:
        if nbseq[c] >= min_nb_seq:
            d.append(x[g == c])
    legends = [
        str(k) + " (" + str(v) + ")" for k, v in {k: nbseq[k] for k in categ}.items()
    ]
    legends = [
        str(k) + " (" + str(v) + ")"
        for k, v in {i: k.size for i, k in enumerate(d)}.items()
    ]
    # axes.axhline(y=x.mean(), color='green', linewidth=1, alpha=.5)
    axes.axhline(y=x.median(), color="blue", linewidth=1, alpha=.5)
    axes.boxplot(d, sym=".", whis=[10, 90])
    axes.set_title("Violin plot of silhouette scores per sample - iteration %d" % n)
    axes.set_xlabel("Cluster ID (nb sequences)")
    # plt.xticks(np.arange(len(categ)) + 1, legends, rotation='vertical')
    plt.xticks(np.arange(len(d)) + 1, legends, rotation="vertical")
    axes.set_ylabel("Silhouette scores")
    plt.tight_layout()
    # plt.show()


def pcacomp_to_model(Pcomp, models, n, min_coeff=0.3, max_pval=0.05, outfile=None):
    """
    Search correlation between a component (usually the first principal component)
    and attributes from a model dict.

    Pearson correlation shows monotonic AND linear relationship between 2
    vectors while Spearman correlation shows monotonic relationship only (see:
    http://www.statisticshowto.com/monotonic-relationship/).

    Example:
        pcacomp_to_model(D_pca[0], {'4mer':X_cpca, 'IND': Y_cpca, 'C2V': Z_cpca},
            min_coeff=0.3, max_pval=0.05, outfile="correl.csv")
    """
    import sys
    import pandas as pd
    from scipy.stats import zscore, pearsonr, spearmanr

    if outfile is not None:
        f = open(outfile, "a")
    else:
        f = sys.stdout
    print(
        ",".join(
            ["#iter,model,modelIdx,SpearmanCoeff,SpearmanPval,PearsonCoeff,PearsonPval"]
        ),
        file=f,
    )
    for m, d in models.items():
        d = pd.DataFrame(zscore(d), index=Pcomp.index).fillna(
            0
        )  # zscore might produce NaN
        nrow, ncol = d.shape
        for i in range(ncol):
            coeff1, pval1 = spearmanr(Pcomp, d[i])
            coeff2, pval2 = pearsonr(Pcomp, d[i])
            print(
                ",".join([str(x) for x in [n, m, i, coeff1, pval1, coeff2, pval2]]),
                file=f,
            )
            # if (abs(coeff) >= min_coeff) and (pval <= max_pval) :
            #     print(','.join([str(x) for x in [m,i,coeff1,pval1,coeff2,pval2]]), file=f)
    if outfile is not None:
        f.close()


def load_models(h5file, models):
    """
    Load sequence models listed in `models` then returns a dict formatted as
    `{"model_name": model}` and the associated must-link matrix.

    If one of the provided `models` is not found, dies.

    Example:
        D_raw, D_ml = load_models(h5file, ['contig2vec4', 'kmers4'])
    """
    import sys
    import pandas as pd
    from fennec import DNASequenceBank

    M = {}
    with pd.HDFStore(h5file) as hdf:
        available_models = hdf.keys()
    for m in models:
        key = "/rawmodel/{}".format(m)
        if key in available_models:
            M[m] = pd.read_hdf(h5file, key)
            idx = M[m].index
        else:
            raise Exception("[ERROR] CANNOT LOAD '{}'".format(key))
            sys.exit(1)
    mlmat = DNASequenceBank._load_sparse_mat(h5file, "_data_mustlink_matrix").todok()
    return M, idx, mlmat
