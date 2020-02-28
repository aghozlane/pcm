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

"""Compute GDT_TS of a set of proteins and produce a csv file report."""

from __future__ import print_function
import argparse
import sys
import os
import shutil
import glob
import subprocess
import time
import re
import multiprocessing as mp
import csv
import numpy as np
import ConfigParser
import pandas as pd

__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2014, INRA"
__credits__ = ["Amine Ghozlane"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@jouy.inra.fr"
__status__ = "Developpement"


class Structureconfig:

    def __init__(self, config_file, results):
        """Instantiate Structureconfig object
           Arguments:
            config_file: Configuration file path
            results: Result path
        """
        self.hdict = {}
        self.structconfig_file = '{0}structconfig.cfg'.format(results)
        self.config = ConfigParser.RawConfigParser()
        if(config_file is not None):
            self.structconfig_file = config_file
            self.readconfig()
        elif(os.path.isfile(self.structconfig_file)):
            self.readconfig()
        else:
            self.writeconfig()
            self.readconfig()

    def readconfig(self):
        """Read config data
        """
        # If config parser empty
        if(not os.path.isfile(self.structconfig_file)):
            self.writeconfig()
        # Read config file
        self.config.read(self.structconfig_file)
        # Get parameter value
        self.hdict["iPBA"] = self.config.get('Structure_config', 'iPBA')
        self.hdict["TMalign"] = self.config.get('Structure_config', 'TMalign')
        self.hdict["maxcluster"] = self.config.get('Structure_config',
                                                   'maxcluster')
        self.hdict["mammoth"] = self.config.get('Structure_config', 'mammoth')

    def writeconfig(self):
        """Write Phylogeny config
        """
        self.config.add_section('Structure_config')
        self.config.set('Structure_config', 'iPBA',
                        "python %path_soft{0}MainParserPDB_align.py "
                        "--p1 %pdb1_name --p2 %pdb2_name "
                        "-r %path_pdb1 --type global -o %output > log_iPBA.txt "
                        "2> error_log_iPBA.txt".format(os.sep))
        self.config.set('Structure_config', 'TMalign', "%path_soft{0}TMalign "
                        "%path_pdb1{0}%pdb1_name %path_pdb2{0}%pdb2_name "
                        "> %output{0}output_alignment "
                        "2> error_log_TMalign.txt".format(os.sep))
        self.config.set('Structure_config', 'maxcluster',
                        "%path_soft{0}maxcluster -e %path_pdb1{0}%pdb1_name "
                        "-p %path_pdb2{0}%pdb2_name -gdt 4 "
                        "> %output{0}output_alignment "
                        "2> error_log_maxcluster.txt".format(os.sep))
        self.config.set('Structure_config', 'mammoth',
                        "%path_soft{0}mammoth -p %path_pdb1{0}%pdb1_name "
                        "-e %path_pdb2{0}%pdb2_name "
                        "> %output{0}output_alignment "
                        "2> error_log_mammoth.txt".format(os.sep))
        # Write data
        try:
            # Writing our configuration file to 'example.cfg'
            with open(self.structconfig_file, 'wt') as configfile:
                self.config.write(configfile)
        except IOError:
            sys.exit("Error : cannot open file %s" % self.structconfig_file)


def isfile(path):
    """Check if path is an existing file.
      Arguments:
          path: Path to the file
    """
    # from Jonathan Barnoud
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def isdir(path):
    """Check if path is an existing file.
      Arguments:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def isfileordir(path):
    """Check if path is an existing file or directory.
      Arguments:
          path: Path to the file
    """
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError("{0} does not exist.".format(path))
    if not os.path.isdir(path) and not os.path.isfile(path):
        raise argparse.ArgumentTypeError("{0} is not a file or a directory"
                                         .format(path))
    return path


def getArguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -d input_dir/ or "
                                     "{0} -q file1.pdb "
                                     "-t file2.pdb".format(sys.argv[0]))
    parser.set_defaults(results=".{0}".format(os.sep))
    parser.add_argument('-q', dest='query', type=isfileordir,
                        help='Query PDB file or directory.')
    parser.add_argument('-t', dest='target', type=isfileordir,
                        help='Target PDB file or directory.')
    parser.add_argument('-d', dest='pdb_dir', type=isdir,
                        help='Directory that contains the PDB files.')
    parser.add_argument('-s', dest='soft', type=str, required=True,
                        nargs='+', choices=["TMalign", "mammoth"],
                        # "iPBA", "maxcluster",
                        help="Select one ore several structural alignment "
                        "software.")
    parser.add_argument('-p', dest='path_alignment', type=isdir, nargs='+',
                        help="Indicate the path to structural alignment "
                        "software.")
    parser.add_argument('-b', dest='nbest', type=int, default=0,
                        help="Select only the n-best alignment based on the "
                        "TMscore (only available with -dq and -dt options).")
    parser.add_argument('-i', dest='position_interest_file', type=isfile,
                        help='Check the positions of interest')
    parser.add_argument('-r', dest='results', type=isdir,
                        help='Path to result directory.')
    parser.add_argument('-n', dest='thread', default=mp.cpu_count(),
                        type=int, help="Number of thread (Default = all cpus "
                        "available will be used).")
    parser.add_argument('-l', dest='family_label', type=str, default="PDB",
                        help='Family label.')
    parser.add_argument('-c', dest='config', type=isfile,
                        help="Configure data.")
    args = parser.parse_args()
    return args, parser


def check_directory(pdb_dir):
    """Get the list of files in the directory
      Arguments:
          pdb_dir: Path to the directory that contains PDB files
      Returns: List of PDB files.
    """
    return glob.glob('{0}{1}*.pdb'.format(pdb_dir, os.sep))


def runCommand(cmd):
    """Run command
      Arguments:
          cmd: Command to run
    """
#    output, errors = None, None
    try:
        # execute command line
        retcode = subprocess.call(cmd, shell=True)
        if retcode == None:
            sys.exit("Child was terminated")
#        pipe =  subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
#                                 stderr=subprocess.STDOUT)
#        #No return case
#        if not pipe:
#            sys.exit("Child was terminated")
#        output, errors = pipe.communicate(input=input)
    except OSError, e:
        sys.exit("Execution failed: {0}".format(e))
#    return output


def remove_output_folder(list_folders):
    """Remove temporary directories
    """
    try:
        for i in list_folders:
            shutil.rmtree(i)
    except OSError, e:
        sys.exit("Error during the remove of {0} folder".format(e))


def set_parameter(cmd, path_soft, pdb1, pdb2, path_pdb1, path_pdb2, path_temp):
    """Adjust command parameters
    """
    cmd = cmd.replace('%path_soft', path_soft)
    cmd = cmd.replace('%pdb1_name', pdb1)
    cmd = cmd.replace('%pdb2_name', pdb2)
    cmd = cmd.replace('%output', path_temp)
    cmd = cmd.replace('%path_pdb1', path_pdb1)
    cmd = cmd.replace('%path_pdb2', path_pdb2)
    return cmd


def create_CA_file(pdb_file, out):
    """Get CA atoms of the first chain
    """
    CA_file = out + os.sep + os.path.basename(pdb_file)
    chain = ""
    numatom = 0
    try:
        with open(CA_file, "wt") as CA:
            with open(pdb_file, "rt") as pdb:
                for line in pdb:
                    if(line[0:4] in "ATOM" and line[13:15] in "CA"):
                        newchain = line[21:22]
                        newnumatom = line[23:26]
                        # Get only the first chain and no duplicate atom
                        if ((chain == newchain or chain == "")
                            and newnumatom != numatom):
                            CA.write(line)
                            chain = newchain
                            numatom = newnumatom
    except IOError:
        sys.exit("Error : cannot open file {0}".format(CA_file))


def pdb_set(pdb_element1, pdb_element2, soft, path_soft, regex, cmd_init):
    """
    """
    cm = []
    pdb1 = os.path.basename(pdb_element1)
    pdb2 = os.path.basename(pdb_element2)
    out = "temp_soft_{0}_{1}".format(
            ".".join(pdb1.split(".")[:-1]),
            ".".join(pdb2.split(".")[:-1]))
    os.mkdir(out)
    if soft == "maxcluster":
        create_CA_file(pdb_element1, out)
        create_CA_file(pdb_element2, out)
        cm = [set_parameter(cmd_init, path_soft, pdb1, pdb2,
                               out, out), out, regex, 0]
    elif soft == "iPBA":
        if(os.path.dirname(pdb_element1)
           != os.path.dirname(pdb_element2)):
            shutil.copyfile(pdb_element2,
                            os.path.dirname(pdb_element1)
                            + pdb2)
        cm = [set_parameter(cmd_init, path_soft, pdb1, pdb2,
                           os.path.dirname(pdb_element1),
                           os.path.dirname(pdb_element2),
                           out), out, regex, 0]
    else:
        cm = [set_parameter(cmd_init, path_soft, pdb1, pdb2,
                           os.path.dirname(pdb_element1),
                           os.path.dirname(pdb_element2),
                           out), out, regex, 0]
    return out, cm


def map_run(soft, path_soft, list_pdbfiles, list_query, list_target, cmd_init):
    """Create temp directories and list command
    """
    cmd = []
    list_folders = []
    regex = [soft]
    # First regex detect iPBA result, 2 and 3 detect TMalign,
    # 4 is for maxcluster
    if soft == "iPBA":
        # 1line 2 in total elements "1l2e"
        regex += [re.compile(">4:.*RMSD\s*(\S+):.*GDT_TS\s*(\S+)\s"), [1, 2]]
    elif soft == "TMalign":
        # "4p4e"
        regex += [re.compile(
            "^Length\s+of\s+Chain_[0-9]+:\s+(\S+)\s+residues|"
            "^Aligned\s+length=.+RMSD=\s+(\S+),.*\n|"
            "^TM-score=\s+(\S+)\s+\(if normalized by length of Chain_1\)\s+|"
            "^TM-score=\s+(\S+)\s+\(if normalized by length of Chain_2\)\s+"),
            [4, 4], re.compile("(^[A-Z\-]+)\n")]
    elif soft == "maxcluster":
        # "2l4e"
        regex += [re.compile("^Iter\s+1:.+RMSD=\s+([0-9.]+),\s+MAXSUB=([0-9.]+)"
                            "\..+TM=([0-9.]+)\s+|"
                            "^GDT=([0-9.]+)\s+"), [2, 4]]
    elif soft == "mammoth":
        # "2l2e"
        regex += [re.compile("^Z-score= (\S+)\s+|^TM-score= ([0-9.]+)\s+"),
                  [2, 2], re.compile("^Prediction\s+(\S+)\n|"
                                "Experiment\s+(\S+)\n")]
    try:
        if len(list_pdbfiles) > 0:
            for i in xrange(0, len(list_pdbfiles)):
                for j in xrange(i + 1, len(list_pdbfiles)):
                    out, cm = pdb_set(list_pdbfiles[i], list_pdbfiles[j], soft,
                                      path_soft, regex, cmd_init)
                    cmd += [cm]
                    list_folders += [out]
        else:
            for query in list_query:
                for target in list_target:
                    out, cm = pdb_set(query, target, soft, path_soft, regex,
                                      cmd_init)
                    cmd += [cm]
                    list_folders += [out]
    except OSError as e:
        sys.exit("Error during the creation of {0} folder".format(e))
    return list_folders, cmd


def extract_info(aln_file, regex):
    """Parse result file
    """
    seqlen = []
    result = []
    insert = 1
    line = 1
    alignment = ["", ""]
    count = 0
    try:
        with open(aln_file, "rt") as aln:
            for i in aln:
                # Match score
                match_score = regex[1].match(i)
                if match_score:
                    if regex[0] == "TMalign" and len(seqlen) < 2:
                        # special TMalign
                        seqlen += [int(match_score.group(1))]
                        insert = 2
                    else:
                        for j in xrange(insert, regex[2][1] + 1):
                            try:
                                result += [float(match_score.group(j))]
                            except TypeError:
                                if line < regex[2][0]:
                                    line += 1
                                else:
                                    print("Overflow", file=sys.stderr)
                        insert = j
                # Match sequence
                if len(result) >= 2:
                    match_sequence = regex[3].match(i)
                    if match_sequence:
                        # Take only the first or the third line
                        # UGLY
                        if regex[0] == "mammoth":
                            if count == 0 or count == 3:
                                alignment[count % 2] += match_sequence.group((count % 2) + 1)
                            if count == 3:
                                count = 0
                            else:
                                count += 1
                        elif regex[0] == "TMalign":
                            alignment[count] += match_sequence.group(1)
                            count += 1
            assert(len(alignment) == 2 and
                   len(alignment[0]) == len(alignment[1])
                   and len(alignment[0]) > 0)
    except IOError:
        print("Error : cannot open file {0}".format(aln_file),
              file=sys.stderr)
    except ValueError as e:
        sys.exit("Something went wrong during the conversion of {0}"
                 .format(e))
    except AssertionError:
        print("There is something wrong in the parsing of the different "
                 "sequence :{0}{1}{0}Check the file : {2}"
                 .format(os.linesep, alignment, aln_file), file=sys.stderr)
    # Select minimum length score for tmAlign score
    try:
        if seqlen[0] < seqlen[1]:
            result = result[0:2]
        else:
            result = [result[0], result[2]]
    except IndexError:
        pass
    assert(len(result) >= 2)
    return result + alignment


def adjust_composition(cmd):
    """Compensate iPBA bug by changing pdb order
    """
    cmd = cmd.replace("--p1", "--px")
    cmd = cmd.replace("--p2", "--p1")
    cmd = cmd.replace("--px", "--p2")
    return cmd


def align_structures(cmd):
    """Subprocess job
    """
    result = []
    runCommand(cmd[0])
#    result_data = runCommand(cmd[0])
    try:

        result = extract_info(cmd[1] + os.sep + "output_alignment", cmd[2])
#        result = extract_info(result_data)
    except AssertionError:
        # Retry in the opposite way
        if(cmd[3] is not 1):
            newcmd = adjust_composition(cmd[0])
            if newcmd is not cmd[0]:
                print("Retry with mixed PDB.")
                cmd[3] = 1
                cmd[0] = newcmd
                result = align_structures(cmd)
            else:
                print("1. AssertionError for {0}".format(cmd[0]),
                      file=sys.stderr)
                result = [np.nan] * 4
        else:
            print("2. AssertionError for {0}".format(cmd[0]), file=sys.stderr)
            result = [np.nan] * 4
    if not result:
        print("There is still something wrong here")
        print(cmd[0])
        result = [np.nan] * 4
    return result


def compute_struct_aln(soft, path_soft, list_pdbfiles, list_query, list_target,
                       cmd, thread):
    """
    """
    print("Start extracting sequences")
    # generator not available because multiprocessing pickling
    list_folders, list_run = map_run(soft, path_soft, list_pdbfiles,
                                     list_query, list_target, cmd)
    print("Extraction done. {0} comparisons to do.".format(len(list_run)))

    startTime = time.time()
#     resultList = []
#     for i in list_run:
#         resultList += [align_structures(i)]
    pool = mp.Pool(processes=thread)
    asyncResult = pool.map_async(align_structures, list_run)
    resultList = asyncResult.get()
    print("Comparison done. Time {0}s.".format(time.time() - startTime))
    remove_output_folder(list_folders)
    return resultList


def report(output, list_pdbfiles, pdb_set, struct_matrix, score, nbest,
           header_position, family_label):
    """Write detailed results
    """
    k = 0
    elements = len(score)
    try:
        with open(output, 'wt') as f:
            #if header_position:
            #    writer.writerow(["PDB_1", "PDB_2"] + score + header_position)
            #else:
            #    writer.writerow(["PDB_1", "PDB_2"] + score)
            if len(list_pdbfiles) > 0:
                writer = csv.writer(f, delimiter='\t')
                for i in xrange(0, len(list_pdbfiles)):
                    for j in xrange(i + 1, len(list_pdbfiles)):
                        try:
                            writer.writerow(
                                [".".join(os.path.basename(list_pdbfiles[i]).split(".")[:-1]),
                                 ".".join(os.path.basename(list_pdbfiles[j]).split(".")[:-1]),
                                 family_label]
                                + struct_matrix[k][:elements])
                        except TypeError:
                            writer.writerow(
                                [".".join(os.path.basename(list_pdbfiles[i]).split(".")[:-1]),
                                 ".".join(os.path.basename(list_pdbfiles[j]).split(".")[:-1]),
                                 family_label]
                                + struct_matrix[:elements])
                        except IndexError:
                            writer.writerow(
                                [".".join(os.path.basename(list_pdbfiles[i]).split(".")[:-1]),
                                 ".".join(os.path.basename(list_pdbfiles[j]).split(".")[:-1]),
                                 family_label]
                                + struct_matrix[:elements])
                        k += 1
            else:
                result = pd.DataFrame()
                for query in pdb_set:
                    pdb_set[query].sort(key=lambda x: x[2], reverse=True)
                    df = pd.DataFrame(pdb_set[query])
                    df = df.sort_values(by=2, ascending=False)
                    if nbest > 0:
                        short_set =df[0: nbest]
                    else:
                        short_set = df
                    short_set.insert(0, "query", query)
                    short_set.insert(1, "family tabel", family_label)
                    result = result.append(short_set, ignore_index=True)
                result.to_csv(f, index=False, header=False, sep='\t')
    except IOError:
        sys.exit("Error : cannot open file {0}".format(output))


def mean_table(tab):
    """Compute mean, ignore nan values
    """
    return np.mean(np.ma.masked_array(tab, np.isnan(tab)))
    # return float(sum(tab)) / float(len(tab))


def get_mean_values(struct_matrix, metric):
    """Compute mean values
    """
    try:
        for i in xrange(len(metric)):
            print("Mean {0} : {1:.4}".format(metric[i],
                                             mean_table([values[i]
                                                         for values in
                                                         struct_matrix])))
    except ZeroDivisionError:
        print("No comparison done.")
    except TypeError:
        print("{0[0]} : {1:.4}, Mean {0[1]} : {2:.4}".format(metric,
                                                             struct_matrix[0],
                                                             struct_matrix[1]))

def extract_path(soft, path_software, num):
    """
    """
    result = ""
    try:
        result = path_software[num]
    except IndexError:
        print("List software empty for {0}".format(soft))
    return result


def path_empty(soft, software, num):
    """
    """
    return "."


def read_position(position_interest_file):
    """
    """
    data_position = {}
    try:
        with open(position_interest_file, "rt") as position_interest:
            position_interest_reader = csv.reader(position_interest,
                                                  delimiter="\t")
            # Pass header
            header_position = position_interest_reader.next()[1:]
            for line in position_interest_reader:
                if len(line) > 1:
                    positions = [int(element) for element in line[1:]
                                 if element != "NA"]
                    if sorted(positions) != positions:
                        print("The position should be sorted: {0}"
                              .format(line), file=sys.stderr)
                    data_position[line[0]] = positions
                else:
                    print("The number of column in position_interest_file "
                          "is too short (>=2) : {0}".format(line),
                          file=sys.stderr)
    except IOError:
        sys.exit("Error cannot open {0}".format(position_interest_file))
    except ValueError:
        sys.exit("Error the value in the table should be integers: {0}"
                 .format(line))
    return header_position, data_position


def get_residues(alignment, data_position, starting_position_query,
                 starting_position_ref):
    """Identify the residues of interest in the alignment
     Args :
      alignment: 0 query alignment  1 reference alignment
    """
    data_virtual_position = []
    # Force copy
    temp_position = list(data_position)
    virtual_posit_query = 1
    virtual_posit_ref = 1
    nb_elements_init = len(data_position)
    try:
        for i in xrange(len(alignment[1])):
            if((virtual_posit_ref + starting_position_ref) == temp_position[0]
               and alignment[1][i] != "-"):
                # Remove the position from the set
                temp_position.pop(0)
                if alignment[0][i] == "-":
                    data_virtual_position += ["NA"]
                else:
                    data_virtual_position += ["{0}-{1}".format(
                        alignment[0][i], virtual_posit_query +
                        starting_position_query)]
            # Increase the count in the sequence
            if alignment[0][i] != "-":
                virtual_posit_query += 1
            if alignment[1][i] != "-":
                virtual_posit_ref += 1
            if len(temp_position) == 0:
                break
        # The program has to find the exact same number of position
        assert(nb_elements_init == len(data_virtual_position))
    except AssertionError:
        if len(alignment[0]) > 0 and len(alignment[1]) > 0:
            sys.exit("Error not all the sequence have been found for the "
                     "following sequence :{0}Query :{0}{1[0]}Target :{0}{1[1]}"
                     .format(os.linesep, alignment))
    return data_virtual_position


def extract_position_pdb(pdb_file):
    """Get the position of the first residue in the pdb file
    """
    start_position = 0
    try:
        with open(pdb_file, "rt") as pdb:
            for line in pdb:
                if line[0:4] == "ATOM":
                    start_position = int(line[22:26])
                    break
        assert(start_position > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(pdb_file))
    except ValueError:
        sys.exit("The program failed to get the start position "
                 "in the pdb file : {0}".format(pdb_file))
    except AssertionError:
        sys.exit("Error the start position of the pdb file {0} is 0"
                 .format(pdb_file))
    return start_position


def get_starting_position(list_pdb):
    """Build the dictionary that identify the starting position of 
    the sequence
    """
    starting_position_pdb = {}
    for pdb in list_pdb:
        starting_position_pdb[".".join(
            os.path.basename(pdb).split(".")[:-1])] = extract_position_pdb(pdb)
    return starting_position_pdb


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    header_position = None
    data_position = None
    list_pdbfiles = []
    starting_position_pdb = {}
    # Load arguments
    args, parser = getArguments()
    # Check configuration data
    conf_data = Structureconfig(args.config, args.results)
    # Get pdb
    if args.pdb_dir:
        list_pdbfiles = check_directory(args.pdb_dir)
    elif args.query and args.target:
        if os.path.isdir(args.query):
            list_query = check_directory(args.query)
        else:
            list_query = [args.query]
        if os.path.isdir(args.target):
            list_target = check_directory(args.target)
        else:
            list_target = [args.target]
    else:
        print("Please, indicate a PDB file or a directory "
        "that contains PDB files.")
        sys.exit(parser.print_help())
    if args.position_interest_file:
        (header_position,
         data_position) = read_position(args.position_interest_file)
        # Get starting position in the PDB file
        # For the query
        starting_position_pdb = get_starting_position(list_query)
        # For the reference
        starting_position_pdb.update(get_starting_position(list_target))
#         print(starting_position_pdb)
    # Check software path
    if args.path_alignment:
        get_path_alignment = extract_path
    else:
        get_path_alignment = path_empty
    # Statistics
    numsoft = 0
    for soft in args.soft:
        pdb_set = {}
        if soft == "iPBA":
            metric = ["RMSD_iPBA", "GDT_TS_iPBA"]
        elif soft == "TMalign":
            metric = ["RMSD_TMalign", "TM-score_TMalign"]
        elif soft == "maxcluster":
            metric = ["RMSD_maxcluster", "Maxsub_maxcluster",
                      "TM-score_maxcluster", "GDT_maxcluster"]
        elif soft == "mammoth":
            metric = ["Z-score_mammoth", "TM-score_mammoth"]
        else:
            print("Please, indicate the software to use.")
            sys.exit(parser.print_help())
        path_alignment = get_path_alignment(soft, args.path_alignment, numsoft)
        numsoft += 1
        struct_matrix = compute_struct_aln(soft, path_alignment, list_pdbfiles,
                                           list_query, list_target,
                                           conf_data.hdict[soft], args.thread)
        if len(list_pdbfiles) == 0:
            k = 0
            elements = len(metric)
            aa_interest = []
            for query in list_query:
                for target in list_target:
                    pdb_query = ".".join(os.path.basename(query).split(".")[:-1])
                    pdb_target = ".".join(os.path.basename(target).split(".")[:-1])
                    # Grab position information in the alignment
                    if data_position and pdb_target in data_position:
                        aa_interest = get_residues(
                                        struct_matrix[k][elements:],
                                        data_position[pdb_target],
                                        starting_position_pdb[pdb_query],
                                        starting_position_pdb[pdb_target])
                    if pdb_query in pdb_set:
                        pdb_set[pdb_query] += [[pdb_target] +
                                               struct_matrix[k][:elements] +
                                               aa_interest]
                    else:
                        pdb_set[pdb_query] = [[pdb_target] +
                                              struct_matrix[k][:elements] +
                                              aa_interest]
                    k += 1
        # print(pdb_set)
        report(args.results + "struct_matrix_{0}.tsv".format(soft),
               list_pdbfiles, pdb_set, struct_matrix, metric, args.nbest,
               header_position, args.family_label)
        get_mean_values(struct_matrix, metric)


if __name__ == '__main__':
    main()
