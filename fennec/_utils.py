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


def boxplot(x, g, n, min_nb_seq=25):
    """
    Plot a violon plot using `x` values into classes defined by `g`.
    """
    import numpy as np
    from matplotlib import pyplot as plt
    from collections import Counter

    _, axes = plt.subplots(figsize=(11, 8))
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


def draw_2D_plot(
    D,
    n,
    title="",
    labels=None,
    n_iter=1500,
    metric="cosine",
    cmap="tab20",
    figsize=(9, 8),
    outfile="tmp.tsne.data.csv",
    figfmt="png",
    verbose=0,
    force=False,
    sklearn_tsne=False,
):
    """
    Visualization
    """
    assert hasattr(D, "index"), "D is not a pandas.DataFrame?"

    import pandas as pd
    from matplotlib import pyplot as plt
    from os import cpu_count

    tsne_args = {
        "n_components": 2,
        "init": "random",
        "n_iter": n_iter,
        "metric": metric,
        "verbose": verbose,
    }

    try:
        if force:
            raise Exception("Forcing!")
        D_tsne = pd.read_csv(outfile, index_col=0)
    except:
        if sklearn_tsne:
            from sklearn.manifold import TSNE
        else:
            from MulticoreTSNE import MulticoreTSNE as TSNE

            tsne_args["n_jobs"] = cpu_count()

        from scipy.stats import zscore
        from sklearn.decomposition import PCA

        if D.shape[1] > 30:  # reduce dimensionality to 30 if needed
            if verbose >= 2:
                print(f"[PCA] reducing from {D.shape[1]} to 30 attributes")
            D_pca = PCA(30).fit_transform(zscore(D))
        else:
            D_pca = D.as_matrix()

        tsne = TSNE(**tsne_args)
        D_tsne = pd.DataFrame(
            tsne.fit_transform(D_pca), index=D.index, columns=("x", "y")
        )
        D_tsne.to_csv(outfile)  # write data for next call

    D_tsne.shape
    plt.figure(figsize=figsize, dpi=300)
    plt.title(f"{title} - step {n}")
    plt.scatter(D_tsne.x, D_tsne.y, s=3, cmap=cmap, c=labels)


def pcacomp_to_model(Pcomp, models, n, outfile=None):
    """
    Search correlation between a component (usually the first principal component)
    and attributes from a model dict.

    Pearson correlation shows monotonic AND linear relationship between 2
    vectors while Spearman correlation shows monotonic relationship only (see:
    http://www.statisticshowto.com/monotonic-relationship/).

    Example:
        pcacomp_to_model(D_pca[0], {'4mer': X_cpca, 'IND': Y_cpca, 'C2V': Z_cpca},
            outfile="correl.csv")
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
    for m, dd in models.items():
        d = dd.copy().reindex(Pcomp.index)
        d = pd.DataFrame(zscore(d), index=Pcomp.index).fillna(
            0
        )  # zscore might produce NaN
        for i in range(d.shape[1]):
            coeff1, pval1 = spearmanr(Pcomp, d[i])
            coeff2, pval2 = pearsonr(Pcomp, d[i])
            print(
                ",".join([str(x) for x in [n, m, i, coeff1, pval1, coeff2, pval2]]),
                file=f,
            )
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
            # sys.exit(1)
    mlmat = DNASequenceBank._load_sparse_mat(h5file, "_data_mustlink_matrix").todok()
    return M, idx, mlmat


def myKernelPCA(
    X,
    inertia,
    kernel="cosine",
    index=None,
    t=5,
    min_comp=5,
    n_jobs=1,
    verbose=False,
):
    """
    Perform KernelPCA on `X` then keep only `inertia`.

    Parameters:
    -----------
    X: numpy.ndarray or pandas.DataFrame
        Values to transform

    inertia: float
        Percentage of inertia to keep
        0.0 < inertia <= 1.0

    kernel: string (default: "cosine")
        Kernel to use

    index: pandas.core.indexes.base.Index (default: None)
        If `None` all individuals from X are considered. Otherwise, only the
        subset from `index` is considered

    min_comp: int (default: 5)
        Minimum number of components to keep

    t: int (default: 4)
        KernelPCA will be fitted on 1/`t` of `X`

    n_jobs: int (default: 1)
        Number of jobs

    verbose: bool (default: False)
        Toggle verbosity

    Output:
    -------
    X_kpca: numpy.ndarray

    """
    import numpy as np
    import pandas as pd

    assert 0.0 < inertia <= 1.0, (
        "Must be: 0.0 < inertia <= 1.0 (got inertia=" + str(inertia) + ")"
    )
    from scipy.stats import zscore
    from sklearn.decomposition import KernelPCA

    if verbose:
        print(f"[INFO] fitting KernelPCA using 1/{t}th of the data")
    kpca = KernelPCA(kernel=kernel, n_jobs=n_jobs, copy_X=False)
    kpca.fit(X.sample(len(X) // t))
    if verbose:
        print(f"[INFO] filtering to keep {100 * inertia}% of inertia")
    explained_variance_ratio_ = (kpca.lambdas_ / kpca.lambdas_.sum()).cumsum()
    nb_tokeep = max(np.count_nonzero(explained_variance_ratio_ < inertia), min_comp)
    if verbose:
        print(f"[INFO] will keep {nb_tokeep} components")
    X_kpca = kpca.transform(X)[:, 0:nb_tokeep]
    if verbose:
        print(f"[INFO] scaling with zscore")
    X_kpca = zscore(X_kpca)
    return pd.DataFrame(X_kpca, index=X.index)


def merge_models(models, index, kpca_params={"n_jobs": 1}):
    """
    Merge mutliple DNA sequence models. First, each model is processed using
    myKernelPCA with `kpca_params`, then they are concatenated.
    Principal components representing 99.99% of inertia are then extracted using
    PCA to discard "duplicate" attributes between models.

    If the first component represents more 99.99% of inertia, the PCA will not
    return any component. Thus the PCA will be forced to return 5 PC.
    """
    import pandas as pd
    from sklearn.decomposition import PCA

    kmodels = {}
    # jobs = {}
    # #- apply KernelPCA(0.9) to each raw sequence model (asynchronously)
    # with ProcessPoolExecutor(max_workers=kpca_params['n_jobs']) as pe:
    #     for i, d in models.items():
    #         d = d.reindex(index)
    #         jobs[i] = pe.submit(myKernelPCA, d, **kpca_params)
    # #- gather results
    # for i, j in jobs.items():
    #     kmodels[i] = j.result()
    # del jobs
    for i, d in models.items():
        d = d.reindex(index)
        kmodels[i] = myKernelPCA(d, **kpca_params)
    # - concatenate all kernelmodels
    D = pd.concat(list(kmodels.values()), ignore_index=True, axis=1)
    D.as_matrix().mean(), D.as_matrix().std()  # should be almost 0 and 1 resp.
    samplesize = min(max(len(D) // 4, 3333), len(D))
    print(
        "[INFO] Selecting components using PCA using %d individuals (%.2f%%)."
        % (samplesize, (100 * samplesize / len(D)))
    )
    try:
        # - Pick principal components (99.99% of inertia) to discard "duplicate" attributes between models
        pca = PCA(0.9999)
        D_pca = pca.fit(D.sample(samplesize)).transform(D)
        _, n_comp = D_pca.shape  # may raise an Exception
    except Exception as ee:
        print(f"Got exception: '{ee}'. Will force PCA to produce 5 components.")
        pca = PCA(5)
        D_pca = pca.fit(D.sample(samplesize)).transform(D)
        _, n_comp = D_pca.shape

    D_pca = pd.DataFrame(D_pca, index=index)
    return D_pca, pca.components_, pca.explained_variance_ratio_, n_comp


def run_vbgmm(
    D_pca,
    pca_components,
    vbgmm_input_dir,
    max_cluster,
    min_length,
    n,
    seed=None,
    maxiter=500,
    epsilon=1e-6,
    write_converg=False,
    verbose=0,
):
    """
    Run CONCONCT's VBGMM clustering algorithm.

    Parameters:
    -----------

    D_pca: pandas.DataFrame
        Results of PCA(0.75).fit_transform(D)

    pca_components: numpy.ndarray
        Components of PCA. Can be obtained with `PCA(0.75).fit(D).components_`

    vbgmm_input_dir: str
        Directory where data will be written

    max_cluster: int
        Max number of cluster

    min_length: int
        Minimum sequence length

    n: int
        Current iteration number

    seed: int (default: None)
        VBGMM seed (random if None)

    maxiter: int (default: 500)
        Maximum number of iteration for VBGMM

    epsilon: float (default: 1e-6)
        Minimum epsilon between 2 VBGMM iterations

    write_converg: bool (default: False)
        Write convergence info to files

    verbose: int (default: 0)
        Verbose level
    """
    import datetime, os, time, vbgmm, numpy as np, pandas as pd
    from random import randint

    PCA_FILE_BASE = f"{vbgmm_input_dir}PCA_transformed_data_gt{min_length}.csv"
    PCA_COMPONENTS_FILE_BASE = (
        f"{vbgmm_input_dir}PCA_components_data_gt{min_length}.csv"
    )
    FLOAT_FORMAT = "%1.8e"
    seed = seed or randint(2, 22222)
    if verbose >= 1:
        print("[INFO] Writing data before clustering")
    # - save PCA data
    D_pca.to_csv(PCA_FILE_BASE, float_format=FLOAT_FORMAT, index_label="contig_id")
    np.savetxt(
        PCA_COMPONENTS_FILE_BASE, pca_components, fmt=FLOAT_FORMAT, delimiter=","
    )
    # - clustering on the selected features
    if verbose >= 1:
        print(
            "[INFO] VBGMM clustering will be run %d times in parallel"
            % vbgmm.get_n_jobs()
        )
    start = time.time()
    vbgmm.fit(
        vbgmm_input_dir,  # folder contained PCA data
        max_cluster,  # max cluster
        min_length,  # minimum length (for filenames only)
        seed,  # seed
        maxiter,  # max iter
        epsilon,  # epsilon (=tolerance)
        write_converg,  # write convergeance data?
    )
    end = time.time()
    if verbose >= 1:
        print("[INFO] Clustering done in: %s" % datetime.timedelta(seconds=end - start))
    # - load clustering results
    vbgmm_clus = (
        pd.read_csv(
            f"{vbgmm_input_dir}clustering_gt{min_length}.csv", header=None, index_col=0
        )
        .reindex(index=D_pca.index)[1]
        .as_matrix()
    )
    # - rename VBGMM output files for traceability
    os.rename(PCA_FILE_BASE, f"{PCA_FILE_BASE}_step{n}")
    os.rename(PCA_COMPONENTS_FILE_BASE, f"{PCA_COMPONENTS_FILE_BASE}_step{n}")
    return vbgmm_clus


def reassign_tiny_cluster_mustlink(vbgmm_clus, D_ml, verbose=False):
    """
    Reassign tiny cluster individuals to "bigger" clusters according to must-link
    relationship from `D_ml`.
    """
    import numpy as np
    from collections import Counter

    res = vbgmm_clus.copy()
    # np.unique(vbgmm_clus); np.unique(vbgmm_clus).size
    # - step 1: identify small clusters
    nbseq = Counter(res)
    min_nb_seq = max(
        np.percentile(list(nbseq.values()), 15), 50
    )  # 15th percentile of number of seq, or 25
    # cluster_to_keep = np.unique([c for i, c in enumerate(res) if (nbseq[c] > min_nb_seq)])
    # dist_densities = {c: plot_density(pdist(D_pca[res==c], metric="cosine")) for c in cluster_to_keep}
    # - step 2: identify sequences to reassign
    seq_to_reassign = [i for i, c in enumerate(res) if (nbseq[c] <= min_nb_seq)]
    print(len(seq_to_reassign), min_nb_seq)
    for i, seqid in enumerate(seq_to_reassign):
        # - step 3: reassign according to must-links
        if verbose:
            print([i, seqid, res[seqid], nbseq[res[seqid]]], end=" ")
        if D_ml[:, seqid].count_nonzero() > 1:  # if `seqid` links to other than itself
            if verbose:
                print("\tø", end=" ")
            grp, _ = np.where(D_ml[res != -1, seqid].toarray())
            cnt = Counter(res[grp])
            res[seqid] = max(cnt, key=cnt.get)
        elif False:  # - step 4: reassign according to proximity with cluster centers
            pass
            ### REASSIGN A SEQUENCE TO THE CLOSEST CLUSTER
            # Does word because:
            # - sequence are already far from other clusters
            # - the KDTree does not accept cosine distance!!
            # - very few sequences were reassigned using this method during tests on S (3/15)
            # -
            # from sklearn.neighbors import KDTree
            # kdt = KDTree(D_pca[res != -1])
            # seq_to_reassign = np.where(res == -1)[0]
            # len(seq_to_reassign)
            # nb_reassigned = 0
            # for seqid in seq_to_reassign:
            #     print(f" -- Sequence {seqid} {'-'*40}")
            #     dist, ind = kdt.query(D_pca.iloc[seqid].values.reshape(-1,1).T, k=15+1)
            #     dist = dist[:,1:] ; ind = ind[:,1:]  # drop itself from neighborhood
            #     dist = dist.reshape(15) ; ind = ind.reshape(15)
            #     closeclust = res[ind]
            #     # for elts in zip(ind, dist, closeclust):
            #     #     print("%10s  %-20s  %s" % elts)
            #     cnt = Counter(closeclust)
            #     cnt
            #     if cnt.most_common()[0][1] >= 6: # 3/5 même cluster
            #         new_clus = cnt.most_common()[0][0]
            #         nb_reassigned += 1
            #     else:
            #         new_clus = -1
            #     print(f" {' '*35} ~> reassigned to {new_clus}")
            # nb_reassigned / len(seq_to_reassign)
            # - step 5: just discard
        else:
            if verbose:
                print("\t.", end=" ")
            res[seqid] = -1
        if verbose:
            print(f"\t{res[seqid]}")
    return res


def extract_cluster_silhouette(
    curated_vbgmm_clus, D, max_cluster, n, min_nb_seq, pdf, plot=True, verbose=False
):
    """
    Keep a cluster if:
     - silhouette Q1 is higher than global sithoulette median
     - silhouette 10th percentile is higher than 0
     - silhouette median is high than global silhouette median
    """
    assert hasattr(D, "index"), "D is not a pandas.DataFrame?"
    import numpy as np
    import pandas as pd
    from collections import Counter
    from matplotlib import pyplot as plt
    from sklearn.metrics import silhouette_samples

    validatedclus = {}
    # -- silhouette analysis
    silsamp = silhouette_samples(D, curated_vbgmm_clus, metric="cosine")

    if verbose:
        print("[INFO] mean silhouette = %.4f" % silsamp.mean())
        print("[INFO] median silhouette = %.4f" % np.median(silsamp))

    CLUS = pd.DataFrame(data=list(zip(curated_vbgmm_clus, silsamp)), index=D.index)
    CLUS.columns = ("cluster", "silhouette")

    # x = curated_vbgmm_clus
    # g = silsamp
    # - group sequence by cluster ID
    nbseq = Counter(curated_vbgmm_clus)
    d = []

    for c in nbseq.keys():
        d.append(CLUS.silhouette[CLUS.cluster == c])

    global_median = CLUS.silhouette.median()
    treated_ids = []
    nb_good = 0
    for i in range(len(d)):
        print("[INFO] Cluster %i;" % i, end="")
        q1, p10 = np.percentile(d[i], (25, 10))
        if q1 < global_median:
            print("Q1 (%.2f) is lower than global median (%.2f)" % (q1, global_median))
        elif p10 < 0:
            print("P10 (%.2f) is lower than 0" % (p10))
        else:
            print(
                "GOOD (nb_seq=%d (%d); q1=%.2f (%.2f); p10=%.2f (%.2f)"
                % (d[i].size, min_nb_seq, q1, global_median, p10, 0)
            )
            validatedclus[len(validatedclus)] = d[i].index
            treated_ids.extend(d[i].index)
            max_cluster -= 1  # 1 less cluster to search for!
            nb_good += 1

    if plot:
        boxplot(CLUS.silhouette, CLUS.cluster, n)
        # plt.show()
        pdf.savefig()
        plt.close()
    remaining_ids = [x for x in D.index if not x in treated_ids]
    color = []
    for i in D.index:
        if i in remaining_ids:
            color.append(False)
        else:
            color.append(True)

    return validatedclus, remaining_ids, max_cluster
