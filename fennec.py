#!/usr/bin/env python3

import argparse
import sys

from fennec import __version__ as VERSION, __citation__ as CITATION

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog="fennec",
        description="Fennec",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    parser.add_argument("--citation", action="version", version=CITATION)

    subparsers = parser.add_subparsers()

    # - "model" subparser
    m_parser = subparsers.add_parser(
        "model",
        help="Extract features from sequences",
        add_help=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    m_req = m_parser.add_argument_group("mandatory arguments")
    m_req.add_argument(
        "--input",
        default=argparse.SUPPRESS,
        help="Input file",
        required=True,
        metavar="FASTA",
    )
    m_req.add_argument("--_PROG", default="model", help=argparse.SUPPRESS)

    m_opt = m_parser.add_argument_group("optionnal arguments")
    m_opt.add_argument(
        "--min_length", type=int, default=1000, help="Minimum sequence length to consider"
    )
    m_opt.add_argument("--chunk_size", type=int, default=10000, help="Chunk size")
    m_opt.add_argument(
        "--overlap",
        type=str,
        default="auto",
        help="Overlap between chunks. Must be 'auto' or 0+",
    )
    m_opt.add_argument("--outfile", default="<input.h5>", help="Output file")
    m_opt.add_argument(
        "--verbosity", type=int, default=3, choices=[0, 1, 2, 3, 4], help="Verbosity level"
    )
    m_opt.add_argument("--n_jobs", type=int, default=1, help="Number of CPU to use")
    m_opt.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )

    # - "describe" subparser
    d_parser = subparsers.add_parser(
        "describe",
        help="Describe modelled sequences",
        add_help=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    d_req = d_parser.add_argument_group("mandatory arguments")
    d_req.add_argument(
        "--input", required=True, default=argparse.SUPPRESS, help="Input HDF5 file"
    )
    d_req.add_argument("--_PROG", default="describe", help=argparse.SUPPRESS)

    d_opt = d_parser.add_argument_group("optionnal arguments")
    d_opt.add_argument(
        "--db_size",
        action="store_true",
        help="Print number of sequence fragements in the database",
        default=argparse.SUPPRESS,
    )
    d_opt.add_argument(
        "--list_models",
        action="store_true",
        help="List available models in the database",
        default=argparse.SUPPRESS,
    )
    d_opt.add_argument(
        "--repack",
        action="store_true",
        help=argparse.SUPPRESS,
        # help="Repack the HDF5 file",
        default=argparse.SUPPRESS,
    )
    d_opt.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )

    # - "extract" subparser
    e_parser = subparsers.add_parser(
        "extract",
        help="Extract bins from modelled sequences",
        add_help=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    e_req = e_parser.add_argument_group("mandatory arguments")
    e_req.add_argument(
        "--input", required=True, default=argparse.SUPPRESS, help="Input HDF5 file"
    )
    e_req.add_argument(
        "--models",
        required=True,
        default=["kmers4", "kmers110011", "contig2vec4", "cov_gattaca31"],
        nargs="+",
        help="List of models to use",
        metavar="MODEL",
    )
    e_req.add_argument("--_PROG", default="extract", help=argparse.SUPPRESS)

    e_opt = e_parser.add_argument_group("optionnal arguments")
    e_opt.add_argument("--label", type=str, default="fennec", help="Label")
    e_opt.add_argument(
        "--max_iter", type=int, default=25, help="Maximum number of iteration"
    )
    e_opt.add_argument(
        "--max_cluster", type=int, default=600, help="Maximum number of cluster"
    )
    e_opt.add_argument(
        "--kpca_inertia",
        type=float,
        default=0.85,
        help="Inertia to keep after kernel PCA",
        metavar="[0.0-1.0]",
    )
    e_opt.add_argument(
        "--kpca_t",
        type=float,
        default=0.33,
        help="Proportion of data to use to fit kernel PCA",
        metavar="[0.0-1.0]",
    )
    e_opt.add_argument(
        "--ppmode",
        type=str,
        default="reassigntiny",
        choices=["nopostprocessing", "reassigntiny"],
        help="Postprocessing mode",
    )
    e_opt.add_argument(
        "--verbosity", type=int, default=3, choices=[0, 1, 2, 3, 4], help="Verbosity level"
    )
    e_opt.add_argument(
        "--min_cluster_size",
        type=int,
        default=50,
        help=argparse.SUPPRESS,
        # help="Minimum number of sequence per cluster",
    )
    e_opt.add_argument("--n_jobs", type=int, default=1, help="Number of CPU to use")
    e_opt.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )

    args = parser.parse_args()
    print(args)

    if not args.__dict__:  # print usage if there is no args
        parser.error("No argument given")
    elif args._PROG == "model":
        print("== model")
        sys.exit(0)
    elif args._PROG == "describe":
        print("== describe")
        sys.exit(0)
    elif args._PROG == "extract":
        print("== extract")
        sys.exit(0)
    else:
        print("== ERROR ==")
        sys.exit(1)
