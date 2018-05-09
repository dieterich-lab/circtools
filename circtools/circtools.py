#!/usr/bin/env python3

# Copyright (C) 2017 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import sys
import os.path

# global settings
version = "1.2.0-beta"
program_name = "circtest"


# samtools/git like parsing from http://chase-seibert.github.io/blog/2014/03/21/python-multilevel-argparse.html


def main():
    CircTools()


def product_range(p):
    try:
        a, b = map(int, p.split(','))
        return a, b
    except:
        raise argparse.ArgumentTypeError("Product must be a range x->y as in x,y")


class CircTools(object):

    def __init__(self):
        parser = argparse.ArgumentParser(
            description="circtools: a modular, python-based framework for circRNA-related tools that unifies "
                        "several functions in single command line driven software.",
            usage="""circtools [-V] <command> [<args>]
            
            Available commands:

               enrich:       circular RNA RBP enrichment scan
               primex:       circular RNA primer design tool
               detect:       circular RNA detection with DCC
               reconstruct:  circular RNA reconstruction with FUCHS
               circtest:     circular RNA statistical testing
               exon:         circular RNA alternative exon analysis
               quickcheck:   circular RNA sequencing library quick checks               
            """)
        parser.add_argument("command", help="Command to run")

        parser.add_argument("-V",
                            "--version",
                            action="version",
                            version=version
                            )

        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print("The supplied command is unknown")
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    @staticmethod
    def enrich():
        # build the argument list
        parser = argparse.ArgumentParser(
            description="circular RNA RBP enrichment tools")

        # REQUIRED ARGUMENTS
        group = parser.add_argument_group("Required options")

        group.add_argument("-c",
                           "--circ-file",
                           dest="circ_rna_input",
                           help="Path to the CircRNACount file generated by DCC",
                           required=True
                           )

        group.add_argument("-b",
                           "--bed-input",
                           dest="bed_input",
                           help="One or more BED files containing features to overlap",
                           required=True
                           )

        group.add_argument("-a",
                           "--annotation",
                           dest="annotation",
                           help="Genome reference annotation file used to not shuffle into intragenic regions",
                           required=True
                           )

        group.add_argument("-g",
                           "--genome",
                           dest="genome_file",
                           help="Genome file for use with bedtools shuffle. See bedtools man page for details.",
                           required=True
                           )

        # OPTIONAL ARGUMENTS
        group = parser.add_argument_group("Additional options")

        group.add_argument("-o",
                           "--output",
                           dest="output_directory",
                           default="./",
                           help="The output folder for files created by " + program_name + " [default: .]",
                           )

        group.add_argument("-i",
                           "--iterations",
                           dest="num_iterations",
                           help="Number of iterations for CLIP shuffling [default: 1000]",
                           type=int,
                           default=1000
                           )

        group.add_argument("-p",
                           "--processes",
                           dest="num_processes",
                           help="Number of threads to distribute the work to",
                           type=int,
                           default=1
                           )

        group.add_argument("-t",
                           "--temp",
                           dest="tmp_directory",
                           help="Temporary directory used by pybedtools",
                           default="/tmp/"
                           )

        group.add_argument("-T",
                           "--threshold",
                           dest="threshold",
                           help="p-value cutoff",
                           type=int,
                           default=2
                           )

        group.add_argument("-P",
                           "--pval",
                           dest="pval",
                           help="p-value cutoff",
                           type=float,
                           default=0.05
                           )

        # group.add_argument("-H",
        #                    "--header",
        #                    dest="has_header",
        #                    help="Defines if the circRNA input file has a header line [default: no]",
        #                    type=bool,
        #                    default=False
        #                    )

        group.add_argument("-F",
                           "--output-filename",
                           dest="output_filename",
                           help="Defines the output file prefix [default: output]",
                           default="output"
                           )

        group.add_argument("-I",
                           "--include-features",
                           dest="include_features",
                           help="Defines the the features which should be used for shuffling. "
                                "May be specified multiple times. [default: all - shuffle over the whole genome]",
                           # May be used as list: e.g. -I exon -I UTR
                           action='append',
                           )

        group.add_argument("-k",
                           "--keep-temp",
                           dest="keep_temp",
                           help="Keep temporary files created by circtools/bedtools [default: no]",
                           type=bool,
                           default=False
                           )

        args = parser.parse_args(sys.argv[2:])

        # make sure we can load the sub module
        sys.path.append(os.path.join(os.path.dirname(__file__)))

        # start the enrichment module
        import enrichment.enrichment_check
        enrich = enrichment.enrichment_check.EnrichmentModule(args, program_name, version)
        enrich.run_module()

    @staticmethod
    def primex():
        parser = argparse.ArgumentParser(
            description="circular RNA primer design")
        # NOT prefixing the argument with -- means it"s not optional

        group = parser.add_argument_group("Input")

        group.add_argument("-d",
                           "--dcc-file",
                           dest="dcc_file",
                           help="CircCoordinates file from DCC / detect module",
                           required=True
                           )

        group.add_argument("-g",
                           "--gtf-file",
                           dest="gtf_file",
                           help="GTF file of genome annotation e.g. ENSEMBL",
                           required=True
                           )

        group.add_argument("-f",
                           "--fasta",
                           dest="fasta_file",
                           help="FASTA file with genome sequence (must match annotation)",
                           required=True
                           )

        group.add_argument("-O",
                           "--organism",
                           dest="organism",
                           help="Organism of the study (used for primer BLASTing), mm = Mus musculus, hs = Homo sapiens",
                           choices=("mm", "hs"),
                           default="hs"
                           )

        group = parser.add_argument_group("Output options")

        group.add_argument("-o",
                           "--output",
                           dest="output_dir",
                           help="Output directory (must exist)",
                           default="./"
                           )

        group.add_argument("-T",
                           "--title",
                           dest="experiment_title",
                           help="Title of the experiment for HTML output and file name",
                           default="circtools_primer_design"
                           )

        group = parser.add_argument_group("Additional options")

        group.add_argument("-t",
                           "--temp",
                           dest="global_temp_dir",
                           help="Temporary directory (must exist)",
                           default="/tmp/"
                           )

        group.add_argument("-G",
                           "--genes",
                           dest="gene_list",
                           help="Space-separated list of host gene names. Primers for CircRNAs of those genes will be "
                                "designed."
                                "E.g. -G \"CAMSAP1\" \"RYR2\"",
                           required=False,
                           nargs='+'
                           )

        group.add_argument("-p",
                           "--product-size",
                           dest="product_size",
                           help="Space-separated range for the desired PCR product. E.g. -p 80 160 [default]",
                           required=False,
                           default=[80, 160],
                           type=int,
                           nargs='+'
                           )

        group.add_argument("-i",
                           "--id-list",
                           dest="id_list",
                           help="Space-separated list of circRNA IDs."
                                " E.g. -i \"CAMSAP1_9_135850137_135850461_-\" \"CAMSAP1_9_135881633_135883078_-\"",
                           required=False,
                           nargs='+'
                           )

        group.add_argument("-j",
                           "--junction",
                           dest="junction",
                           help="Should the forward [f] or reverse [r] primer be located on the BSJ? [Default: n]",
                           choices=("r", "n", "f"),
                           default="n"
                           )

        group.add_argument("-b",
                           "--no-blast",
                           dest="blast",
                           help="Should primers be BLASTED? Even if selected yes here, not more than 50 primers will"
                                "be sent to BLAST in any case.",
                           default=False,
                           action='store_true'
                           )

        args = parser.parse_args(sys.argv[2:])

        # start the primer module

        # make sure we can load the sub module
        sys.path.append(os.path.join(os.path.dirname(__file__)))

        import primex.primex
        primex_instance = primex.primex.Primex(args, program_name, version)
        primex_instance.run_module()

    @staticmethod
    def detect():
        parser = argparse.ArgumentParser(
            description="circular RNA detection")
        # NOT prefixing the argument with -- means it"s not optional
        parser.add_argument("-C",
                            "--params",
                            dest="cli_params",
                            help="Defines the input parameters for DCC",
                            default="--help"
                            )
        args = parser.parse_args(sys.argv[2:])

        import os
        os.system("DCC " + args.cli_params)

    @staticmethod
    def circtest():
        parser = argparse.ArgumentParser(
            description="circular RNA statistical testing - Interface to https://github.com/dieterich-lab/CircTest")
        # NOT prefixing the argument with -- means it"s not optional

        ######################################################

        group = parser.add_argument_group("Required")
        group.add_argument("-d",
                           "--DCC",
                           dest="DCC_dir",
                           help="Path to the detect/DCC data directory",
                           required=True
                           )

        group.add_argument("-l",
                           "--condition-list",
                           dest="condition_list",
                           help="Comma-separated list of conditions which should be compared"
                                "E.g. \"RNaseR +\",\"RNaseR -\"",
                           required=True
                           )

        group.add_argument("-c",
                           "--condition-columns",
                           dest="condition_columns",
                           help="Comma-separated list of 1-based column numbers in the detect/DCC output"
                                " which should be compared; e.g. 10,11,12,13,14,15",
                           required=True
                           )

        group.add_argument("-g",
                           "--grouping",
                           dest="grouping",
                           help="Comma-separated list describing the relation of the columns specified via -c to the"
                                " sample names specified via -l; e.g. -g 1,2 and -r 3 would assign sample1 to each "
                                "even column and sample 2 to each odd column",
                           required=True
                           )
        ######################################################

        group = parser.add_argument_group("Processing options")

        group.add_argument("-r",
                           "--replicates",
                           dest="num_replicates",
                           help="Number of replicates used for the circRNA experiment [Default: 3]",
                           type=int,
                           default=3
                           )

        group.add_argument("-f",
                           "--max-fdr",
                           dest="max_fdr",
                           help="Cut-off value for the FDR [Default: 0.05]",
                           type=float,
                           default=0.05
                           )

        group.add_argument("-p",
                           "--percentage",
                           dest="percentage",
                           help="The minimum percentage of circRNAs account for the total "
                                "transcripts in at least one group. [Default: 0.01]",
                           type=float,
                           default=0.01
                           )

        group.add_argument("-s",
                           "--filter-sample",
                           dest="filter_sample",
                           help="Number of samples that need to contain the amount of reads "
                                "specified via -C [Default: 3]",
                           type=int,
                           default=3
                           )

        group.add_argument("-C",
                           "--filter-count",
                           dest="filter_count",
                           help="Number of CircRNA reads that each sample specified via -s has to contain "
                                "[Default: 5]",
                           type=int,
                           default=5
                           )

        ######################################################

        group = parser.add_argument_group("Output options")

        group.add_argument("-o",
                           "--output-directory",
                           dest="output_directory",
                           default="./",
                           help="The output directory for files created by " + program_name + " [Default: .]",
                           )

        group.add_argument("-n",
                           "--output-name",
                           dest="output_name",
                           default="circtest",
                           help="The output name for files created by " + program_name + " [Default: circtest]",
                           )

        group.add_argument("-m",
                           "--max-plots",
                           dest="max_plots",
                           help="How many of candidates should be plotted as bar chart? [Default: 50]",
                           type=int,
                           default=50
                           )

        group.add_argument("-a",
                           "--label",
                           dest="label",
                           help="How should the samples be labeled? [Default: Sample]",
                           default="Sample"
                           )

        group.add_argument("-L",
                           "--limit",
                           dest="range",
                           help="How should the samples be labeled? [Default: Sample]",
                           type=float,
                           default=1.0
                           )

        group.add_argument("-O",
                           "--only-negative-direction",
                           dest="only_negative",
                           help="Only print entries with negative direction indicator [Default: False]",
                           default=False
                           )
        group.add_argument("-H",
                           "--add-header",
                           dest="add_header",
                           help="Add header to CSV output [Default: False]",
                           default=False
                           )

        group.add_argument("-M",
                           "--colour",
                           dest="colour",
                           help="Can be set to bw to create grayscale graphs for manuscripts",
                           choices=("colour", "bw"),
                           default="colour"
                           )
        ######################################################

        args = parser.parse_args(sys.argv[2:])

        # start the primer module

        # make sure we can load the sub module
        sys.path.append(os.path.join(os.path.dirname(__file__)))

        import circtest.circtest
        circtest_instance = circtest.circtest.CircTest(args, program_name, version)
        circtest_instance.run_module()

    @staticmethod
    def quickcheck():
        parser = argparse.ArgumentParser(
            description="circular RNA sequencing library quality assessment")
        # NOT prefixing the argument with -- means it"s not optional

        ######################################################

        group = parser.add_argument_group("Required")
        group.add_argument("-d",
                           "--DCC",
                           dest="DCC_dir",
                           help="Path to the detect/DCC data directory",
                           required=True
                           )

        group.add_argument("-s",
                           "--star",
                           dest="star_dir",
                           help="Path to the base STAR data directory containing sub-folders with per-sample mappings",
                           required=True
                           )

        group.add_argument("-l",
                           "--condition-list",
                           dest="condition_list",
                           help="Comma-separated list of conditions which should be compared"
                                "E.g. \"RNaseR +\",\"RNaseR -\"",
                           required=True
                           )

        group.add_argument("-g",
                           "--grouping",
                           dest="grouping",
                           help="Comma-separated list describing the relation of the columns specified via -c to the"
                                " sample names specified via -l; e.g. -g 1,2 and -r 3 would assign sample1 to each "
                                "even column and sample 2 to each odd column",
                           required=True
                           )
        ######################################################

        group = parser.add_argument_group("Output options")

        group.add_argument("-o",
                           "--output-directory",
                           dest="output_directory",
                           default="./",
                           help="The output directory for files created by " + program_name + " [Default: ./]",
                           )

        group.add_argument("-n",
                           "--output-name",
                           dest="output_name",
                           default="quickcheck",
                           help="The output name for files created by " + program_name + " [Default: quickcheck]",
                           )

        group.add_argument("-c",
                           "--colour",
                           dest="colour",
                           help="Can be set to bw to create grayscale graphs for manuscripts",
                           choices=("colour", "bw"),
                           default="colour"
                           )

        group.add_argument("-C",
                           "--cleanup",
                           dest="cleanup",
                           help="String to be removed from each sample name "
                                "[Default: \"_STARmapping.*Chimeric.out.junction\"]",
                           default="_STARmapping.*Chimeric.out.junction"
                           )

        group.add_argument("-S",
                           "--starfolder",
                           dest="starfolder",
                           help="Suffix string of the STAR folders"
                                "[Default: \"_STARmapping\"]",
                           default="_STARmapping"
                           )

        group.add_argument("-L",
                           "--remove-last",
                           dest="remove_suffix_chars",
                           help="Remove last N characters from each column name of the DCC input data "
                                "[Default: 0]",
                           type=int,
                           default=0
                           )

        group.add_argument("-F",
                           "--remove-first",
                           dest="remove_prefix_chars",
                           help="Remove first N characters from each column name of the DCC input data "
                                "[Default: 0]",
                           type=int,
                           default=0
                           )

        group.add_argument("-R",
                           "--remove-columns",
                           dest="remove_columns",
                           help="Comma-separated list of columns in the DCC data files to not includes in the check",
                           default="0"
                           )
        ######################################################

        args = parser.parse_args(sys.argv[2:])

        # start the primer module

        # make sure we can load the sub module
        sys.path.append(os.path.join(os.path.dirname(__file__)))

        import quickcheck.quickcheck
        quickcheck_instance = quickcheck.quickcheck.QuickCheck(args, program_name, version)
        quickcheck_instance.run_module()

    @staticmethod
    def exon():
        parser = argparse.ArgumentParser(
            description="circular RNA exon usage analysis")
        # NOT prefixing the argument with -- means it"s not optional

        ######################################################

        group = parser.add_argument_group("Required")
        group.add_argument("-d",
                           "--DCC",
                           dest="DCC_dir",
                           help="Path to the detect/DCC data directory",
                           required=True
                           )

        group.add_argument("-l",
                           "--condition-list",
                           dest="condition_list",
                           help="Comma-separated list of conditions which should be compared"
                                "E.g. \"RNaseR +\",\"RNaseR -\"",
                           required=True
                           )

        group.add_argument("-c",
                           "--condition-columns",
                           dest="condition_columns",
                           help="Comma-separated list of 1-based column numbers in the detect/DCC output"
                                " which should be compared; e.g. 10,11,12,13,14,15",
                           required=True
                           )

        group.add_argument("-g",
                           "--grouping",
                           dest="grouping",
                           help="Comma-separated list describing the relation of the columns specified via -c to the"
                                " sample names specified via -l; e.g. -g 1,2 and -r 3 would assign sample1 to each "
                                "even column and sample 2 to each odd column",
                           required=True
                           )

        group.add_argument("-r",
                           "--replicates",
                           dest="replicates",
                           help="Comma-separated list describing the relation of the samples specified via -g to the"
                                " sample names specified via -l; e.g. -g 1,2 and -r 3 would assign sample1 to each "
                                "even column and sample 2 to each odd column",
                           required=True
                           )

        group.add_argument("-b",
                           "--ballgown-data",
                           dest="ballgown_data",
                           help="Path to the ballgown data directory",
                           required=True
                           )

        group.add_argument("-G",
                           "--gtf-file",
                           dest="gtf_file",
                           help="Path to the GTF file containing the employed genome annotation",
                           required=True
                           )

        group.add_argument("-C",
                           "--circtest-output",
                           dest="circtest_file",
                           help="Path to the CircTest CSV file containing the CircTest results",
                           required=True
                           )

        ######################################################

        group = parser.add_argument_group("Additional options")

        group.add_argument("-H",
                           "--has-header",
                           dest="has_header",
                           help="Do the CircTest result files have a header? [Default: No]",
                           type=bool,
                           default=False
                           )

        ######################################################

        group = parser.add_argument_group("Output options")

        group.add_argument("-o",
                           "--output-directory",
                           dest="output_directory",
                           default="./",
                           help="The output directory for files created by " + program_name + " [Default: .]",
                           )

        group.add_argument("-n",
                           "--output-prefix",
                           dest="output_prefix",
                           default="exon_analysis_",
                           help="The output name (prefix) for files created by " + program_name +
                                " [Default: exon_analysis]"
                           )

        # group.add_argument("-p",
        #                    "--max-plots",
        #                    dest="max_plots",
        #                    help="How many of candidates should be plotted as bar chart? [Default: 10]",
        #                    type=int,
        #                    default=10
        #                    )

        ######################################################

        args = parser.parse_args(sys.argv[2:])

        # start the primer module

        # make sure we can load the sub module
        sys.path.append(os.path.join(os.path.dirname(__file__)))

        import exon_usage.exon_usage
        exon_instance = exon_usage.exon_usage.ExonUsage(args, program_name, version)
        exon_instance.run_module()

    @staticmethod
    def reconstruct():
        parser = argparse.ArgumentParser(
            description="circular RNA reconstruction")
        # NOT prefixing the argument with -- means it"s not optional
        parser.add_argument("-C",
                            "--params",
                            dest="cli_params",
                            help="Defines the input parameters for DCC",
                            default="--help"
                            )
        args = parser.parse_args(sys.argv[2:])

        import os
        os.system("FUCHS " + args.cli_params)


if __name__ == "__main__":
    main()
