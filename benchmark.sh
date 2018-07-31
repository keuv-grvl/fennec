#!/usr/bin/env bash

# enter a screen session
screen -ls "fen.bench" || screen -dmS "fen.bench"
screen -r "fen.bench"

# increase max proc limit to the maximum
ulimit -Su $(ulimit -Hu)
prlimit

# force environment reactivation
source deactivate "fennec2-dev"
source activate "fennec2-dev"

# ensure parallel is installed
which parallel || exit 1 ;

# reinstall fennec
rm -r build/ fennec.egg-info/ vbgmm.cpython-36m-x86_64-linux-gnu.so
N_RTHREADS=8 pip install -e .  # can run 20 instances simultaneously on totoro

# check fennec version
python --version
python -c "import fennec; print(fennec.__version__)"
python -c "import vbgmm; print(vbgmm.get_n_jobs())"

# get path to GNU time
TIME=$(for A in $(whereis time); do test -x $A && echo $A; done|tail -n1)
test -z $TIME && exit 1 || echo "GNU time path: $TIME" ;  # exit if $TIME is not defined

# parameters to test
ROOT="FENNEC_RESULTS"  # main output folder
ERROR_DIR="EXIT1"
L0="S M L"  # datasets
L1="auto 0"  # overlap
L2="kmeans mustlink"  # clustering initialization
L3="nopostprocessing reassigntiny"  # pipeline postprocessing
L4="kmers4 contig2vec4 kmers110011 ind15 kmers4,contig2vec4,contig2vec6,cov_gattaca31,kmers110011,ind15"  # models

mkdir -p "$ERROR_DIR/"

# generate list of commands to run and folders
for A in $(echo $L0) ; do
    for B in $(echo $L1) ; do
        F="/dev/shm/${USER}/${A}.completedata.l1000c10000o${B}.h5"
        test -e "$F" || continue  # skip if input file does not exist
        for C in $(echo $L2) ; do
            for D in $(echo $L3) ; do
                for E in $(echo $L4) ; do
                    DIR=$ROOT/$A/$B/$C/$D/$E/
                    mkdir -p "$DIR/"
                    # skip if 'DONE' file already exists
                    test -e "$DIR/DONE" && continue
                    # skip if exit status is 0, otherwise archive existing results
                    grep -sq "Exit status: 0" "${DIR}/time.log"  \
                        && continue  \
                        || (DDIR=$(sed 's/\//_/g' <<< $DIR) && tar zcf "${ERROR_DIR}/${DDIR}.tgz" $DIR && rm -r $DIR/*)
                    # build command to run
                    echo $TIME -v -o \"$DIR/time.log\" python3 fennec_cluster_extraction_pipeline.py \"$F\" \"$A\" \"$B\" \"$C\" \"$D\" \"$E\" 2\> \"$DIR/exec.err\" \> \"$DIR/exec.log\"
                done
            done
        done
    done
done > benchmark.cmds

# run commands in parallel
parallel -j16 < benchmark.cmds
echo "pfiou... done"
