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
import glob
import csv
import urllib2
import string

__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2014, INRA"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@jouy.inra.fr"
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
    parser.add_argument('-i', dest='input_file', type=isfile,
                        help="List of PDB of interest (txt or fasta file).")
    parser.add_argument('-p', dest='pdb', type=str, nargs='+', default=[],
                        help="List of pdb files or codes.")
    parser.add_argument('-d', dest='pdb_dir', type=str, nargs='+',
                        help="Directory that contains PDB files.")
    parser.add_argument('-c', dest='chain_int', type=str, default="A",
                        help="Indicate a chain of interest for extraction "
                        "activity (default chain A).")
    parser.add_argument("-a", dest="activity", type=str,
                        choices=["renum", "sequence", "extract", "clean"],
                        help="Select one operation : renum: to renum the "
                        "pdb file, sequence: get the amino-acid sequence or "
                        "clean: renumerotate and on select ATOM.")
    parser.add_argument('-o', dest='output_dir', type=str, required=True,
                        help="Output directory.")
    return parser.parse_args()


def get_numstring(val, maxval):
    """
    """
    strval = str(val)
    strval_len = len(strval)
    if(strval_len < maxval):
        missing = maxval - strval_len
        strval = " " * missing + strval
    elif(strval_len > maxval):
        sys.exit("There is something not expected, val = {0}".format(val))
    return(strval)


def renum_pdb(pdb_file, activity, seqdict, out, out_str):
    """
    """
    res = 0
    all_possibilities = string.ascii_uppercase + string.digits
    # num = start_position - 1
    num = 0
    num_chain = -1
    aa_prev = ""
    chain_prev = ""
    try:
        with open(pdb_file, "rt") as pdb:
            for line in pdb:
                do_not_print = False
                chain = line[20:22]
                aa = line[17:20]
                newres = line[22:26]
                field = line[0:4]
                if chain != chain_prev and field == "ATOM":
                    num = int(newres) - 1
                    chain_prev = chain
                    num_chain += 1 
                if((newres != res or aa != aa_prev) and field == "ATOM"):
                    aa_prev = aa
                    res = newres
                    num += 1
                    if activity == "sequence":
                        print(seqdict[aa], file=out, end="")
                        #sys.stdout.write(seqdict[aa])
                pdb_line = list(line)
                if(field == "ATOM"):
                    # Check the residue number
                    pdb_line[22:26] = get_numstring(num, 4)
                    if activity != "extract":
                        if pdb_line[16] != " ":
                            if pdb_line[16] == "A":
                                pdb_line[16] = " "
                            else:
                                do_not_print = True
                        pdb_line[21] = all_possibilities[num_chain]
                    # Check the chain 
                    assert(num_chain<36)
                pdb_line = "".join(pdb_line)
                if not do_not_print:
                    if activity == "renum":
                        print(pdb_line, file=out, end="")
                    elif (activity == "clean" or activity == "extract") and field == "ATOM":
                        #print(pdb_line, file=out, end="")
                        out_str += pdb_line
                    #sys.stdout.write(pdb_line)
                # print(line[23:26])
            if activity != "clean" and activity != "extract":
                print(os.linesep, file=out, end="")
            #sys.stdout.write("\n")
    except IOError:
        sys.exit("Error cannot open {0}".format(pdb_file))
    except AssertionError:
        sys.exit("Error in chain renumeration. There is more than 36 chains it's strange")
    return out_str


def check_directory(pdb_dir):
    """Get the list of files in the directory
      Arguments:
          pdb_dir: Path to the directory that contains PDB files
      Returns: List of PDB files.
    """
    return glob.glob('{0}{1}*.pdb'.format(pdb_dir, os.sep))


def read_list(input_file):
    """
    """
    pdb_list = []
    try:
        with open(input_file) as input_f:
            input_f_reader = csv.reader(input_f, delimiter="\t")
            for line in input_f_reader:
                pdb_list.append(line[0].replace("\n", "").replace("\r", "").split(" ")[0])
        assert(len(pdb_list) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(input_file))
    except AssertionError:
        sys.exit("No PDB element read from {0}".format(input_file))
    return pdb_list


def read_fasta(input_file):
    """
    """
    pdb_list = []
    try:
        with open(input_file) as input_f:
            for line in input_f:
                if line[0] == ">":
                    pdb_list.append(line[1:].replace("\n", "").replace("\r", "").split(" ")[0])
        assert(len(pdb_list) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(input_file))
    except AssertionError:
        sys.exit("No PDB element read from {0}".format(input_file))
    return pdb_list


def download_pdb(pdb, results):
    """Download PDB file
     :Parameters:
        pdb: Name of the pdb file
        results: Output directory
    """
    # Open our local file for writing
    outfilename = results + os.sep + pdb + ".pdb"
    try:
        rscb = urllib2.urlopen("http://www.rcsb.org/pdb/files/" + pdb + ".pdb")
        with open(outfilename, "wt") as local_file:
            local_file.write(rscb.read())
    except urllib2.HTTPError, errormes:
        print ("HTTP Error:", errormes.code, pdb, file=sys.stderr)
    except urllib2.URLError, errormes:
        print ("URL Error:", errormes.reason, pdb, file=sys.stderr)
    except IOError:
        print("Something went wrong with {0}".format(outfilename),
              file=sys.stderr)
    return outfilename


def get_unique(seq):
    # Order preserving
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]


def detect_chain_A(out_str):
    """Detect if a chain A is present
    """
    chain_types = []
    for line in out_str.split('\n'):
        chain = line[21:22]
        if chain == "A":
            return True, chain_types
        elif chain != "":
            chain_types += [chain]
    chain_types.sort()
    chain_types = get_unique(chain_types)
    return False, chain_types


def set_chain_A(out_str, chain_types, out):
    """Reset chain A
    """
    for line in out_str.split('\n'):
        chain = line[21:22]
        if chain == chain_types[0]:
            pdb_line = list(line)
            pdb_line[21] = "A"
            pdb_line = "".join(pdb_line)
            print(pdb_line, file=out, end="\n")
        else:
            print(line, file=out, end="\n")


def extract_chain(out_str, chain_int, out):
    """Extract a chain
    """
    extraction_done = False
    try:
        for line in out_str.split('\n'):
            # print(line)
            chain = line[21:22]
            # print(chain)
            if chain == chain_int:
                print(line, file=out, end="\n")
                extraction_done = True
        assert(extraction_done)
    except AssertionError:
        print("Chain {} not found in the PDB".format(chain_int), file=sys.stderr)


def main():
    """
    """
    seqdict = {"ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D",
               "CYS":"C", "GLU":"E", "GLN":"Q", "GLY":"G",
               "HIS":"H", "ILE":"I", "LEU":"L", "LYS":"K",
               "MET":"M", "PHE":"F", "PRO":"P", "SER":"S",
               "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V"}
    args = get_arguments()
    # Get PDB element from a fasta or a list
    if args.input_file:
        if (args.input_file.endswith(".fasta")
            or args.input_file.endswith(".faa")):
            args.pdb += read_fasta(args.input_file)
        else:
            args.pdb += read_list(args.input_file)
    # Get the PDB files from a directory
    if args.pdb_dir:
        args.pdb += check_directory(args.pdb_dir)
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)
    temp_dir = args.output_dir + os.sep + "temp"
    if args.activity == "sequence":
        output_type = ".fasta"
    else:
        output_type = ".pdb"
    # Start working on pdb
    toremove = []
    for pdb in args.pdb:
        if not pdb.endswith(".pdb") and not os.path.isfile(pdb):
            # Create a temporary directory
            if not os.path.isdir(temp_dir):
                os.mkdir(temp_dir)
            # Download pdb file in temp directory
            if args.activity == None:
                pdb = download_pdb(pdb, args.output_dir)
                #print(pdb)
            else:
                pdb = download_pdb(pdb, temp_dir)
                toremove.append(pdb)
        pdb_name = ".".join(os.path.basename(pdb).split(".")[:-1])
        output_file = args.output_dir + os.sep + pdb_name + output_type
        if args.activity == "extract":
            output_file = args.output_dir + os.sep + pdb_name + "_" + args.chain_int + output_type
        if args.activity:
            try:
                with open(output_file, "wt") as out:
                    if args.activity == "sequence":
                        print(">{0}".format(pdb_name), file=out)
                    out_str = renum_pdb(pdb, args.activity, seqdict, out, "")
                    if args.activity == "clean":
                        chain_A_ok, chain_types = detect_chain_A(out_str)
                        if chain_A_ok:
                            print(out_str, file=out)
                        else:
                            set_chain_A(out_str, chain_types, out)
                    elif args.activity == "extract":
                        extract_chain(out_str, args.chain_int, out)
            except IOError:
                sys.exit("Error cannot open {0}".format(output_file))
    #print(toremove)
    for temp_pdb in toremove:
        if os.path.isfile(temp_pdb):
            os.remove(temp_pdb)
    if os.path.isdir(temp_dir):
        os.rmdir(temp_dir)

if __name__ == '__main__':
    main()
