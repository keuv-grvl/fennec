#!/usr/bin/env bash

# ensure parallel is installed
which parallel || exit 1

# enter a screen session
screen -dmS fen.bench
screen -r fen.bench

# force environment reactivation
source activate fennec2-dev

# reinstall fennec
rm -r build/ fennec.egg-info/ vbgmm.cpython-36m-x86_64-linux-gnu.so
N_RTHREADS=8 pip install -e .  # can run 20 instances simultaneously on totoro

# get path to GNU time
TIME=$(for A in $(whereis time); do test -x $A && echo $A; done|tail -n1)
test -z $TIME && exit 1 || echo "GNU time path: $TIME" ;  # exit if $TIME is not defined

# parameters to test
ROOT="FENNEC_RESULTS"  # main output folder
L0="S M L"  # datasets
L1="auto 0"  # overlap
L2="kmeans mustlink"  # clustering initialization
L3="nopostprocessing reassigntiny fullpipeline"  # pipeline postprocessing
L4="kmers4 contig2vec4 kmers110011 ind15 kmers4,contig2vec4,contig2vec6,cov_gattaca31,kmers110011,ind15"  # models1

# generate list of commands to run and folders
for A in $(echo $L0) ; do
    for B in $(echo $L1) ; do
        F="/dev/shm/gravouil/${A}.completedata.l1000c10000o${B}.h5"
        test -e "$F" || continue  # skip if input file does not exist
        for C in $(echo $L2) ; do
            for D in $(echo $L3) ; do
                for E in $(echo $L4) ; do
                    test -e "$ROOT/$A/$B/$C/$D/$E/DONE" && continue  # skip if 'DONE' file already exists
                    mkdir -p "$ROOT/$A/$B/$C/$D/$E/"
                    echo $TIME -v -o \"$ROOT/$A/$B/$C/$D/$E/time.log\" python3 fennec_VBGMM_cluster_extraction.py \"$F\" \"$A\" \"$B\" \"$C\" \"$D\" \"$E\" \> \"$ROOT/$A/$B/$C/$D/$E/exec.log\"
                done
            done
        done
    done
done > benchmark.cmds

# run command in parallel
parallel -j16 < benchmark.cmds
echo "pfiou... done"
