#!/usr/bin/env python
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Pre-clean pdb files."""

from __future__ import print_function
import sys
import os
import argparse
import csv
from math import isnan

__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2020, Institut Pasteur"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@pasteur.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    parser.add_argument('-sr', dest='structure_quality_ref_file', type=isfile,
                        required=True, help="Structural quality by ref.")
    parser.add_argument('-sn', dest='structure_quality_neg_file', type=isfile,
                        required=True, help="Structural quality by neg.")
    parser.add_argument('-ar', dest='alignment_ref_file', nargs='+', default=[],
                        required=True, help="Structural alignment by ref "
                        "(TMalign + Mammoth).")
    parser.add_argument('-an', dest='alignment_neg_file', nargs='+', default=[],
                        required=True, help="Structural alignment by neg "
                        "(TMalign + Mammoth).")
    parser.add_argument('-n', dest="name", type=str, default="ard",
                        help="Family name")
    parser.add_argument('-o', dest='output_file', type=str,
                        default="pcm_result.tsv", help="Output file.")
    return parser.parse_args()


def load_data(gene_res, data_file, tag_type):
    """Load quality data of the structure"""
    try:
        with open(data_file, "rt") as data:
            data_reader = csv.reader(data, delimiter="\t")
            # Pass header
            #data_reader.next()
            for line in data_reader:
                gene = line[0].split(".")[0]
                pdb_id = line[0].split(".")[1]
                family = line[1]
                #print(gene, family)
                #print(line[2:])
                if (gene, family) in gene_res:
                    gene_res[gene, family].update({tag_type:line[2:], "pdb_id_"+ tag_type:pdb_id})
                else:
                    gene_res[gene,family] = {tag_type:line[2:], "pdb_id_" + tag_type:pdb_id}
    except IOError:
        sys.exit("Error cannot open {0}".format(data_file))
    return gene_res


def differential(ref, neg):
    """
    """
    assert(len(ref) == len(neg))
    res = []
    for i in xrange(len(ref)):
        #ugly
        if neg[i] == "":
            neg[i] = 0
        if ref[i] == "":
            ref[i] = 0
        refval = float(ref[i])
        negval = float(neg[i])
        if isnan(negval) and isnan(refval):
            res += [0]
        elif isnan(negval):
            res += [refval]
        elif isnan(refval):
            res += [negval]
        else:
            res += [refval - negval]
    return res

def write_result(gene_res, name, output_file):
    """
    """
    try:
        with open(output_file, "wt") as output:
            output_writer = csv.writer(output, delimiter="\t")
            output_writer.writerow(["Sequence", "Type", "pdb_id_ref",
                                    "pdb_id_tneg", "d_molpdf",
                                    "d_dope_score", "d_normalized_dope",
                                    "d_GA341_score", "d_zscore", "d_maxsub",
                                    "d_lgscore", "d_mypmfs", "d_Z-score_mammoth",
                                    "d_TM-score_mammoth", "d_RMSD_TMalign",
                                    "d_TM-score_TMalign", "molpdf_ref",
                                    "dope_score_ref", "normalized_dope_ref",
                                    "GA341_score_ref", "zscore_ref",
                                    "maxsub_ref", "lgscore_ref","mypmfs_ref",
                                    "Z-score_mammoth_ref",
                                    "TM-score_mammoth_ref", "RMSD_TMalign_ref",
                                    "TM-score_TMalign_ref", "molpdf_tneg",
                                    "dope_score_tneg", "normalized_dope_tneg",
                                    "GA341_score_tneg", "zscore_tneg",
                                    "maxsub_tneg", "lgscore_tneg", "mypmfs_tneg",
                                    "Z-score_mammoth_tneg",
                                    "TM-score_mammoth_tneg",
                                    "RMSD_TMalign_tneg",
                                    "TM-score_TMalign_tneg"])
            for gene in gene_res:
                # ugly but necessary
                if "pdb_id_quality_ref" in gene_res[gene]:
                    if "pdb_id_alignment_ref_mammoth" in gene_res[gene]:
                        assert(gene_res[gene]["pdb_id_quality_ref"] == gene_res[gene]["pdb_id_alignment_ref_mammoth"])
                    if "pdb_id_alignment_ref_TMalign" in gene_res[gene]:
                        assert(gene_res[gene]["pdb_id_quality_ref"] == gene_res[gene]["pdb_id_alignment_ref_TMalign"])
                if "pdb_id_quality_neg" in gene_res[gene]:
                    if "pdb_id_alignment_neg_mammoth" in gene_res[gene]:
                        #print(gene_res[gene]["pdb_id_quality_neg"])
                        #print(gene_res[gene]["pdb_id_alignment_neg_mammoth"])
                        assert(gene_res[gene]["pdb_id_quality_neg"] == gene_res[gene]["pdb_id_alignment_neg_mammoth"])
                    if "pdb_id_alignment_neg_TMalign" in gene_res[gene]:
                        assert(gene_res[gene]["pdb_id_quality_neg"] == gene_res[gene]["pdb_id_alignment_neg_TMalign"])
                if not "pdb_id_quality_ref" in gene_res[gene]:
                    if "pdb_id_alignment_ref_mammoth" in gene_res[gene]:
                        gene_res[gene]["pdb_id_quality_ref"] = gene_res[gene]["pdb_id_alignment_ref_mammoth"]
                    elif "pdb_id_alignment_ref_TMalign" in gene_res[gene]:
                        gene_res[gene]["pdb_id_quality_ref"] = gene_res[gene]["pdb_id_alignment_ref_TMalign"]
                    else:
                        gene_res[gene]["pdb_id_quality_ref"] = "nan"
                if not "pdb_id_quality_neg" in gene_res[gene]:
                    if "pdb_id_alignment_neg_mammoth" in gene_res[gene]:
                        gene_res[gene]["pdb_id_quality_neg"] = gene_res[gene]["pdb_id_alignment_neg_mammoth"]
                    elif "pdb_id_alignment_neg_TMalign" in gene_res[gene]:
                        gene_res[gene]["pdb_id_quality_neg"] = gene_res[gene]["pdb_id_alignment_neg_TMalign"]
                    else:
                        gene_res[gene]["pdb_id_quality_neg"] = "nan"
                if not "quality_ref" in gene_res[gene]:
                    gene_res[gene]["quality_ref"] = [0] * 7
                if not "quality_neg" in gene_res[gene]:
                    gene_res[gene]["quality_neg"] = [0] * 7
                if not "alignment_ref_mammoth" in gene_res[gene]:
                    gene_res[gene]["alignment_ref_mammoth"] = [0] * 3
                if not "alignment_neg_mammoth" in gene_res[gene]:
                    gene_res[gene]["alignment_neg_mammoth"] = [0] * 3
                if not "alignment_ref_TMalign" in gene_res[gene]:
                    gene_res[gene]["alignment_ref_TMalign"] = [0] * 3
                if not "alignment_neg_TMalign" in gene_res[gene]:
                    gene_res[gene]["alignment_neg_TMalign"] = [0] * 3
                #print(gene_res[gene]["quality_ref"])
                #print(gene_res[gene]["quality_neg"])
                #print("quality")
                print(gene)
                print(gene_res[gene]["quality_ref"])
                print(gene_res[gene]["quality_neg"])
                diff_quality = differential(gene_res[gene]["quality_ref"],
                                            gene_res[gene]["quality_neg"])
                #print("mammoth")
                diff_mammoth = differential(gene_res[gene]["alignment_ref_mammoth"][1:],
                                            gene_res[gene]["alignment_neg_mammoth"][1:])
                #print("tmalign")
                diff_TMalign = differential(gene_res[gene]["alignment_ref_TMalign"][1:],
                                            gene_res[gene]["alignment_neg_TMalign"][1:])
                output_writer.writerow([gene[0], "Candidate_" + gene[1],
                                        gene_res[gene]["pdb_id_quality_ref"],
                                        gene_res[gene]["pdb_id_quality_neg"]] +
                                       diff_quality + diff_mammoth +
                                       diff_TMalign +
                                       gene_res[gene]["quality_ref"] +
                                       gene_res[gene]["alignment_ref_mammoth"][1:] +
                                       gene_res[gene]["alignment_ref_TMalign"][1:] +
                                       gene_res[gene]["quality_neg"] +
                                       gene_res[gene]["alignment_neg_mammoth"][1:] +
                                       gene_res[gene]["alignment_neg_TMalign"][1:])
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))
    #except AssertionError:
    #    sys.exit("Error it is not the pdb_id between quality and alignment for ref or net")


def main():
    """
    """
    args = get_arguments()
    gene_res = {}
    gene_res = load_data(gene_res, args.structure_quality_ref_file,
                         "quality_ref")
    gene_res = load_data(gene_res, args.structure_quality_neg_file,
                         "quality_neg")
    #print(gene_res)
    tag_type = ["ref"] * 2 + ["neg"] *2
    for i, data_file in enumerate(args.alignment_ref_file +
                                  args.alignment_neg_file):
        soft = os.path.basename(data_file).split("_")[0]
        gene_res = load_data(gene_res, data_file, "alignment_{0}_{1}"
                             .format(tag_type[i], soft))
        #print(gene_res)
    write_result(gene_res, args.name, args.output_file)

if __name__ == '__main__':
    main()
