#!/usr/bin/env python3

import gc
import os
import sys

# import matplotlib
# matplotlib.use('Agg')  # https://stackoverflow.com/questions/37604289

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.spatial.distance import pdist
from sklearn import preprocessing

import fennec
from fennec._utils import (
    boxplot,
    draw_2D_plot,
    extract_cluster_silhouette,
    isinteractive,
    list_models,
    load_models,
    merge_models,
    myKernelPCA,
    pcacomp_to_model,
    reassign_tiny_cluster_mustlink,
    run_vbgmm,
)

# - check is we are in the correct conda env
if os.environ["CONDA_DEFAULT_ENV"] != "fennec2-dev":
    raise Exception("Not in correct conda environment")


# -- functions in development ---------------------------------------------------


def _nonredondant_pairwise_index(x):
    for i in range(len(x)):
        for j in range(i, len(x)):
            yield i, j


def _count_must_link_np(ci, cj):
    """Count number of mustlink for cluster `ci` and `cj`"""
    import numpy as np

    global D_ml, curated_vbgmm_clus
    return D_ml[
        np.ix_(
            np.where(curated_vbgmm_clus == ci)[0], np.where(curated_vbgmm_clus == cj)[0]
        )
    ].count_nonzero()


def _mysquareform(x, l):
    """Thanks to @pidupuis"""
    import numpy as np

    m = np.zeros([l, l]).astype(int)
    xs, ys = np.triu_indices(l)
    m[xs, ys] = m[ys, xs] = x
    return m


def get_nb_mustlink_per_cluster(curated_vbgmm_clus, D_ml, n_jobs=1, verbose=False):
    """
    Return number of mustlink for each pair of cluster
    """
    ## (incredibly ugly) list comprehension
    # return pd.DataFrame(data=_mysquareform([
    #         D_ml[
    #             np.ix_(
    #                 np.where(curated_vbgmm_clus == ci)[0],
    #                 np.where(curated_vbgmm_clus == cj)[0]
    #             )
    #         ].count_nonzero()
    #         for ci, cj in _nonredondant_pairwise_index(categ)],
    #         l=len(categ)).reshape(len(categ), len(categ)),
    #     index=categ, columns=categ)
    import numpy as np
    import pandas as pd
    from time import sleep
    from concurrent.futures import ProcessPoolExecutor

    categ = np.unique(curated_vbgmm_clus)
    jobs = []
    with ProcessPoolExecutor(max_workers=n_jobs) as pe:
        for ci, cj in _nonredondant_pairwise_index(categ):
            jobs.append(pe.submit(_count_must_link_np, categ[ci], categ[cj]))
        if verbose:
            while True:
                nb_f = 0
                for j in jobs:
                    if j.done():
                        nb_f += 1
                if nb_f >= len(jobs):
                    print("100% finised")
                    break
                else:
                    perc_f = nb_f / len(jobs) * 100
                    print("%05.2f%%" % perc_f, end="\r")
                sleep(1)
    cnt_ml = pd.DataFrame(
        data=_mysquareform([j.result() for j in jobs], l=len(categ)).reshape(
            len(categ), len(categ)
        ),
        index=categ,
        columns=categ,
    )
    return cnt_ml


def extract_unlink_clusters(curated_vbgmm_clus, D_ml, tol=0.9, verbose=True):
    """
    If 2 clusters from `curated_vbgmm_clus` have at least 90% of must-link link from
    `D_ml`, they are dropped from the cluster list and set as `remaining_ids`.
    """
    import numpy as np

    cnt_ml = get_nb_mustlink_per_cluster(
        curated_vbgmm_clus, D_ml, verbose=verbose, n_jobs=8
    )
    perc_ml = np.diag(cnt_ml) / cnt_ml.sum()
    pd.concat([cnt_ml, perc_ml])
    cnt_ml.to_csv(f"{vbgmm_input_dir}/iter{n}_mustlink.csv")
    clustokeep = cnt_ml.index[perc_ml >= tol]
    remaining_ids = cnt_ml.index[perc_ml < tol]

    # return values
    # validated clusters
    retA = curated_vbgmm_clus[np.isin(curated_vbgmm_clus, clustokeep)]
    # remainings ids
    retB = curated_vbgmm_clus[np.isin(curated_vbgmm_clus, remaining_ids)]
    # max nb of cluster
    retC = len(retA)
    return retA, retB, retC


# -------------------------------------------------------------------------------#
# -- SCRIPT STARTS HERE ---------------------------------------------------------#
# -------------------------------------------------------------------------------#

if isinteractive():  # debug args if script is run in python shell
    DATASET = "XS"
    models_str = "kmers4,contig2vec4,contig2vec6,cov_gattaca31,kmers110011"  # kmers5,
else:
    if len(sys.argv) != 3:
        raise Exception(
            'usage: python3 fennec_VBGMM_cluster_extraction.py <DATASET> "contig2vec4,ind15,kmers4,kmers110011,kmers1100110011,cov_gattaca31"'
        )
    else:
        _, DATASET, models_str = sys.argv

# -- user input
vbgmm_input_dir = f"run.{DATASET}.output/"
min_length = 1000  # minimum sequence length
max_cluster = 300  # maximum number of cluster to extract
max_iter = 15  # maximum number of iteration

# -- variables
wanted_models = models_str.split(",")
print(f"Models: {wanted_models}")
h5file = f"DATA/{DATASET}.completedata.h5"

if not os.path.isfile(h5file):
    # TODO: extract feature from fasta file if h5file does not exist
    print("[ERROR] can not find file '%s'. Exiting." % h5file)
    sys.exit(1)
else:
    # Is the input a fasta file? In this case, extract features and set 'h5file'
    # See: `fennec_sequence_characterization.py`
    pass

# -- load data
raw_models, remaining_ids, D_ml = load_models(h5file, wanted_models)
print([(i, d.shape[1]) for i, d in raw_models.items()])

# -- set some parameters
kpca_params = {
    "inertia": 0.85,
    "n_jobs": os.cpu_count(),
    "verbose": 3,
}  # see: help(myKernelPCA)
min_nb_seq = 50  # default value, will be updated later
max_pca_dim = min(250, sum([d.shape[1] for d in raw_models.values()]) // 3)

final_cluster = {}  # final cluster (id: [seqid])
HISTORY = {}  # store data, pca, clusterings, filtered clustered, etc.
n = 0  # current iteration

# -- open report files
os.makedirs(vbgmm_input_dir, exist_ok=True)  # create output dir
pdf = PdfPages(
    vbgmm_input_dir + "/vbgmm_iterative_extraction_" + DATASET + ".pdf",
    keep_empty=False,
)

# -- dev parameters
devmode = False  # if True, stop at each iteration
force_gc = True  # force garbage collection at the end of each iteration
draw_plot = True  # draw 2D plot of each clustering (take some times)

# --------------------------- #
# -- main loop starts here -- #
# --------------------------- #

while True:
    print(f"[INFO] --- ITERATION {n} {'-'*60}")

    # -- check if we have to continue
    if n >= max_iter:
        print("[END] It's already been %d iterations! Exiting." % n)
        final_cluster["maxiter_" + str(n)] = remaining_ids
        break

    if max_cluster <= 1:
        print("[END] I will search for only %d clusters Exiting." % max_cluster)
        final_cluster["lastcluster_" + str(n)] = remaining_ids
        break

    # -- select features
    print("[INFO] Merging models")
    D, pca_components, pca_explained_variance_ratio, n_comp = merge_models(
        raw_models, remaining_ids, kpca_params
    )
    print("[INFO] Merging models produced %d components" % n_comp)

    with open(f"{vbgmm_input_dir}/pca_explvarratio.txt", "a") as outfile:
        print(n, pca_explained_variance_ratio, file=outfile)

    # -- check if we have to continue
    if len(remaining_ids) < min_nb_seq:
        print(
            "[END] There is only %d sequences remaining (minimum %d per cluster). Exiting."
            % (len(D), min_nb_seq)
        )
        final_cluster["notenough_" + str(n)] = remaining_ids
        break

    if n_comp > max_pca_dim:
        print(
            "[END] PCA produced %d components for %.2f %% of inertia (maximum is %d). Exiting."
            % (n_comp, 100 * 0.9999, max_pca_dim)
        )
        final_cluster["unbinned_" + str(n)] = remaining_ids
        break

    # TODO: Check for clusterability before clustering (see: https://github.com/alimuldal/diptest)

    # -- which models are used?
    os.makedirs(f"{vbgmm_input_dir}/comporigins", exist_ok=True)
    for c in range(min(20, n_comp)):  # 20 first components
        pcacomp_to_model(
            D[c],
            raw_models,
            n,
            outfile=f"{vbgmm_input_dir}/comporigins/pca_components_origin_comp{c}.csv",
        )

    # -- clustering
    vbgmm_clus = run_vbgmm(
        D,
        pca_components,
        D_ml,
        vbgmm_input_dir,
        max_cluster,
        min_length,
        n,
        # seed=666,
        epsilon=1e-4,
        verbose=2,
    )
    assert vbgmm_clus.shape[0] == len(
        remaining_ids
    ), "[ERROR] Not all sequences were clustered!"

    # -- clustering post processing
    # - drop tiny cluster (reassign sequences eventually)
    curated_vbgmm_clus = reassign_tiny_cluster_mustlink(vbgmm_clus, D_ml, verbose=True)
    np.unique(curated_vbgmm_clus).size

    # - TODO: merge clusters with enough must-link relationships
    # linked_vbgmm_clus = extract_unlink_clusters(
    #     curated_vbgmm_clus, D_ml, verbose=True
    # )  ## NOTE: CURRENTLY DOES NOTHING

    # - extract cluster with high enough silhouette scores
    validatedclus, remaining_ids, max_cluster = extract_cluster_silhouette(
        curated_vbgmm_clus, D, max_cluster, n, min_nb_seq, pdf, verbose=True
    )

    # -- draw plot using t-SNE
    if draw_plot:
        draw_2D_plot(
            D,
            n,
            labels=vbgmm_clus,
            title="raw_vbgmm_clus",
            n_iter=1200,
            force=True,
            verbose=3,
        )
        plt.savefig(f"{vbgmm_input_dir}/iter{n}.tsne.rawclus.png")
        draw_2D_plot(
            D,
            n,
            labels=curated_vbgmm_clus,
            title="curated_vbgmm_clus",
            n_iter=1200,
            force=False,
            verbose=3,
        )
        plt.savefig(f"{vbgmm_input_dir}/iter{n}.tsne.curatedclus.png")

    CLUS = pd.DataFrame(
        data=list(curated_vbgmm_clus), index=D.index, columns=("cluster",)
    )

    HISTORY[n] = (
        # DATASET
        D.copy(),
        # pca ratio
        pca_explained_variance_ratio.copy(),
        # raw vbgmm clustering
        vbgmm_clus.copy(),
        # curated clustering
        curated_vbgmm_clus.copy(),
        # extracted clusters
        validatedclus,
        remaining_ids,
        max_cluster,
        #
        CLUS.copy(),
    )

    if len(validatedclus) == 0:
        print("[END] No good cluster have been extracted Exiting.")
        # final_cluster['unbinned_' + str(n)] = remaining_ids
        # - accept clustering as is
        for seqid in remaining_ids:
            try:
                final_cluster["badbinned_" + str(CLUS.cluster[seqid])].append(seqid)
            except:
                final_cluster["badbinned_" + str(CLUS.cluster[seqid])] = list()
                final_cluster["badbinned_" + str(CLUS.cluster[seqid])].append(seqid)
        break

    # - TODO: stop if more than 90~95% sequence untreated?

    # - prepare for next iteration
    n += 1
    kpca_params["index"] = remaining_ids
    if force_gc:
        gc.collect()
    if devmode:
        break


# ---------------------- #
# -- end of main loop -- #
# ---------------------- #

print("[INFO] done in %d iterations" % (n + 1))
pdf.close()

# -- write clusters to CSV file
output_cluter_file = f"{vbgmm_input_dir}/final_clustering.csv"
cnt = 0

with open(output_cluter_file, "w") as f:
    print("%s,%s" % ("#SEQID", "CLUSTER"), file=f)
    # - print extracted clusters from HISTORY
    for i, (_, _, _, _, vc, _, _, _) in HISTORY.items():
        for k, v in vc.items():
            for seqid in v:
                print(f"{seqid},{cnt}", file=f)
            cnt += 1
    # - print unbinned or badbinned seqids
    for i, seqids in final_cluster.items():
        for seqid in seqids:
            print(f"{seqid},{i}", file=f)
        cnt += 1


print("[INFO] Bye.")
sys.exit(0)

# ----------------------------- #
# -- Fin ---------------------- #
# ----------------------------- #

# -- crawl history
(
    D,
    pca_explained_variance_ratio,
    vbgmm_clus,
    curated_vbgmm_clus,
    validatedclus,
    remaining_ids,
    max_cluster,
    CLUS,
) = HISTORY[0]
