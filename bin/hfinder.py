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

"""Find homologous proteins and extract them in a multifasta."""

from __future__ import print_function
import ConfigParser
import argparse
import sys
import os
import subprocess
import re
import bisect
import multiprocessing as mp
import time
import csv


__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2020, Institut Pasteur"
__credits__ = ["Amine Ghozlane"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@pasteur.fr"
__status__ = "Developpement"


class Inconfig:

    def __init__(self, config_file, results):
        """Instantiate Inconfig object
           Arguments:
            config_file: Configuration file path
            results: Result path
        """
        self.hdict = {}
        self.inconfig_file = '{}{}hconfig.cfg'.format(results, os.sep)
        self.config = ConfigParser.RawConfigParser()
        if(config_file is not None):
            self.inconfig_file = config_file
            self.readconfig()
        elif(os.path.isfile(self.inconfig_file)):
            self.readconfig()
        else:
            self.writeconfig()
            self.readconfig()


    def readconfig(self):
        """Read config data
        """
        # If config parser empty
        if(not os.path.isfile(self.inconfig_file)):
            self.writeconfig()
        # Read config file
        self.config.read(self.inconfig_file)
        # Get parameter value
        self.hdict["ssearch"] = self.config.get('Homology_config','ssearch')
        self.hdict["jackhmmer"] = self.config.get('Homology_config',
                                                  'jackhmmer')
        #self.hdict["hhsearch"] = self.config.get('Homology_config', 'hhsearch')
        self.hdict["psiblast"] = self.config.get('Homology_config', 'psiblast')
        self.hdict["blastp"] = self.config.get('Homology_config', 'blastp')
        self.hdict["blastx"] = self.config.get('Homology_config', 'blastx')
        self.hdict["tblastn"] = self.config.get('Homology_config', 'tblastn')
        self.hdict["hmmsearch"] = self.config.get('Homology_config', 'hmmsearch')
        self.hdict["hmmscan"] = self.config.get('Homology_config', 'hmmscan')
        self.hdict["mmseqs"] = self.config.get('Homology_config', 'mmseqs')
        self.hdict["clustalo"] = self.config.get('Alignment_config', 'clustalo')


    def writeconfig(self):
        """Write Inconfig config
        """
        self.config.add_section('Homology_config')
        # -O option is not working...
        self.config.set('Homology_config', 'ssearch',
                        "%path_softssearch %query %database "
                        "-E %e_value -p -m 8 -T %proc -s BP62 > %output")
        self.config.set('Homology_config', 'jackhmmer',
                        "%path_softjackhmmer --tblout %output --cpu %proc "
                        "-E %e_value %query %database > log_jackhmmer.txt")
        #self.config.set('Homology_config', 'hhsearch', "nothing")
        self.config.set('Homology_config', 'psiblast', "%path_softpsiblast "
                        "-query %query -db %database -out %output "
                        "-evalue %e_value -outfmt 6 -num_threads %proc "
                        "-num_iterations 3")
        self.config.set('Homology_config', 'blastp', "%path_softblastp "
                        "-query %query -db %database -out %output "
                        "-evalue %e_value -outfmt 6 -num_threads %proc "
                        "-max_target_seqs 10")
        self.config.set('Homology_config', 'tblastn', "%path_softtblastn "
                        "-query %query -db %database -out %output "
                        "-evalue %e_value -outfmt 6 -num_threads %proc "
                        "-max_target_seqs 10")
        self.config.set('Homology_config', 'blastx', "%path_softblastx "
                        "-query %query -db %database -out %output "
                        "-evalue %e_value -outfmt 6 -num_threads %proc "
                        "-max_target_seqs 10")
        self.config.set('Homology_config', 'hmmsearch', "%path_softhmmsearch "
                        " --tblout %output --cpu %proc -E %e_value "
                        "%hmm_db %query > log_hmmsearch.txt")
        self.config.set('Homology_config', 'hmmscan', "%path_softhmmscan "
                        " --tblout %output --cpu %proc -E %e_value "
                        "%hmm_db %query > log_hmmscan.txt")
        self.config.set('Homology_config', 'mmseqs', "%path_softmmseqs "
                        "easy-search %query %database %output %temporary --threads %proc "
                        " --search-type 1")
        self.config.add_section('Alignment_config')
        self.config.set('Alignment_config', 'clustalo',
                        "%path_softclustalo -i %multifasta -o %output "
                        "--auto -t Protein --outfmt=fa")
        # Write data
        try:
            # Writing our configuration file to 'example.cfg'
            with open(self.inconfig_file, 'wt') as configfile:
                self.config.write(configfile)
        except IOError:
            sys.exit("Error : cannot open file {0}".format(self.inconfig_file))


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
            msg = "{0} is a file.".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def getArguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "hfinder.py -q query.fasta "
                                     "-d database.fasta")
    parser.add_argument('-q', dest='query', type=isfile, required=True,
                        help='Path to the query file.')
    parser.add_argument('-qm', dest='query_hmm', type=isfile,
                        help='Path to the database file.')
    parser.add_argument('-d', dest='database', type=isfile, required=True,
                        help='Path to the database file (fasta).')
    parser.add_argument('-db', dest='database_blast', type=str,
                        help='Path to the database file (blast).')
    parser.add_argument('-s', dest='software', type=str,
                        choices=["blastp", "tblastn", "blastx", "psiblast",
                                 "jackhmmer", "ssearch", "hmmsearch", "hmmscan",
                                 "mmseqs"],
                        nargs='+', default=["blastp", "jackhmmer", "ssearch"],
                        help='Select protein homology software.')
    parser.add_argument('-p', dest='path_software', type=isdir,
                        nargs='+', help='Path to the software.')
    parser.add_argument('-w', dest='path_clustalo', type=isdir, default=None,
                        help='Path to clustalo software.')
    #parser.add_argument('-rev', dest='reverse', action='store_true',
                        #default=False)
    parser.add_argument('-r', dest='results', type=isdir, default="." + os.sep,
                        help='Path to result directory.')
    parser.add_argument('-n', dest='nbest', type=int, default=0,
                        help='Number of best.')
    parser.add_argument('-t', dest='thread', type=int, default=mp.cpu_count(),
                        help='Number of thread.')
    parser.add_argument('-tmp', dest='tmp', type=isdir, default="." + os.sep + "tmp",
                        help='Temporary folder for mmseqs.')
    parser.add_argument('-b', dest='behavior', default=[],
                        choices=["force_computation", "extract",
                                 "cumulative", "check"],
                        nargs='+', help='Behavior : force computation (if '
                        'previous result), extract sequence, check identity.')
    parser.add_argument('-e', dest='e_value', type=float, default=0.001,
                        help='E-value threshold.')
    parser.add_argument('-idmin', dest='filter_identity_minimum', type=float,
                        default=None,
                        help='Filter on identity minimum threshold.')
    parser.add_argument('-idmax', dest='filter_identity_maximum', type=float,
                        default=None,
                        help='Filter on identity maximum threshold.')
    parser.add_argument('-simmin', dest='filter_similarity_minimum', type=float,
                        default=None,
                        help='Filter on similarity minimum threshold.')
    parser.add_argument('-simmax', dest='filter_similarity_maximum', type=float,
                        default=None,
                        help='Filter on similarity maximum threshold.')
    parser.add_argument('-covmin', dest='filter_coverage_minimum', type=float,
                        default=None,
                        help='Filter on coverage minimum threshold.')
    parser.add_argument('-covmax', dest='filter_coverage_maximum', type=float,
                        default=None,
                        help='Filter on coverage maximum threshold.')
    parser.add_argument('-lmin', dest='filter_length_minimum', type=int,
                        default=None,
                        help='Filter on sequence length minimum threshold.')
    parser.add_argument('-lmax', dest='filter_length_maximum', type=int,
                        default=None,
                        help='Filter on sequence length maximum threshold.')
    parser.add_argument('-c', dest='config', type=isfile,
                        help='Path to configuration file.')
    args = parser.parse_args()
    return args


def run_command(cmd):
    """Run command
      Arguments:
          cmd: Command to run
    """
    try:
        # Execution de la ligne de commande
        retcode = subprocess.call(cmd, shell=True)
        # Cas aucun retour du soft
        if retcode == None:
            sys.exit("Child was terminated")
    except OSError, e:
        sys.exit("Execution failed: {0}".format(e))
    except:
        sys.exit("There is something wrong with the command: {0}".format(cmd))


def replace_motif(build_command, path_soft, query, database, database_hmm,
                  output, thread, e_value, tmp):
    """
    """
    print(build_command, file=sys.stderr)
    if path_soft:
        build_command = build_command.replace('%path_soft', path_soft)
    else:
        build_command = build_command.replace('%path_soft', "")
    build_command = build_command.replace('%proc', str(thread))
    build_command = build_command.replace('%query', query)
    build_command = build_command.replace('%database', database)
    if not database_hmm:
        database_hmm = ""
    build_command = build_command.replace('%hmm_db', database_hmm)
    build_command = build_command.replace('%output', output)
    build_command = build_command.replace('%e_value', str(e_value))
    build_command = build_command.replace('%temporary', tmp)
    print(build_command, file=sys.stderr)
    return build_command


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
    return ""


def extract_homology(output_file, soft): #, reverse):
    """
    """
    homologous = []
    interest = 1
    #if reverse:
        #interest = 0
    #regex = get_regex_by_soft(soft)
    try:
        with open(output_file, "rt") as output:
            if soft == "hmmsearch" or soft == "hmmscan":
                if soft == "hmmsearch":
                    regex_hmm = re.compile(r"(\S+)\s+")
                else:
                    regex_hmm = re.compile(r"\S+\s+\S+\s+(\S+)\s+")
                for line in output:
                    if line[0] != "#":
                        hmm_match = regex_hmm.match(line)
                        if hmm_match:
                            homologous += [hmm_match.group(1)]
            else:
                output_reader = csv.reader(output, delimiter="\t")
                for line in output_reader:
                    #if len(line) == 12:
                    if len(line) > 1:
                        homologous += [line[interest]]
    except IOError:
        sys.exit("Error : cannot read {0}".format(output_file))
    return homologous


def get_unique(data_list):
    """Get unique data
      Arguments:
         data_list: list of data
    """
    return {}.fromkeys(data_list).keys()


def get_element(input_list, name):
    """Search name in input list
      Arguments:
        input_list: List
        name: Search criteria
    """
    # Searching the node with its name
    i = bisect.bisect_left(input_list, name)
    # Object has been found
    if(i != len(input_list) and input_list[i] == name):
        return True
    return False


def register_fasta(results, query_file, listhomology, output):
    """
    """
    reader = False
    try:
        with open(output, "wt") as multifasta_out:
            # Extract homologous
            with open(query_file, "rt") as query:
                for line in query:
                    if line.startswith(">"):
                        reader = False
                        title = line[1:].replace("\n", "").replace("\r", "")
                        if " " in title:
                            title = title.split(" ")[0]
                        # Select homologous
                        if get_element(listhomology, title):
                            multifasta_out.write(line)
                            reader = True
                    elif reader and len(line) > 0:
                        multifasta_out.write(line)
    except IOError as e:
        sys.exit("Something went wrong with the output\n{0}".format(e))


def create_fasta(compair):
    """
    """
    try:
        with open(compair[5], 'wt') as fasta_file:
            fasta_file.write(">{1[0]}{0}{1[1]}{0}"
                             ">{1[2]}{0}{1[3]}{0}".format(os.linesep, compair))
    except IOError:
        sys.exit("Error: cannot open file {0}".format(compair[5]))


def remove_temp_files(temp_files):
    """
    """
    try:
        for temp in temp_files:
            os.remove(temp)
    except OSError, e:
        print("Error: cannot remove file {0}".format(e), file=sys.stderr)


def extract_sequence(multifasta_file, listhomology=None):
    """
    """
    data = []
    seq = False
    try:
        with open(multifasta_file, "rt") as multifasta:
            for line in multifasta:
                if line.startswith(">"):
                    reader = False
                    title = line[1:].replace("\n", "").replace("\r", "")
                    if " " in title:
                        title = title.split(" ")[0]
                    if listhomology:
                        if get_element(listhomology, title):
                            data += [[title, ""]]
                            seq = True
                        else:
                            title = ""
                            seq = False
                    else:
                        data += [[title, ""]]
                        seq = True
                elif seq and len(line) > 0:
                    data[-1][1] += line.replace("\n", "").replace("\r", "")
                else:
                    seq = False
            assert(len(data) > 0)
    except IOError:
        sys.exit("Error : cannot open file {0}".format(multifasta_file))
    except AssertionError:
        sys.exit("Error : no sequence extracted from {0}".format(multifasta_file))
    return data


def replace_sequence(cmd, path_soft, multifasta_file, alignment_file):
    """
    """
    if path_soft:
        cmd = cmd.replace('%path_soft', path_soft)
    else:
        cmd = cmd.replace('%path_soft', "")
    cmd = cmd.replace('%multifasta', multifasta_file)
    cmd = cmd.replace('%output', alignment_file)
    return cmd


def map_sequence(results, cmd, query, database, listhomology, path_clustalo):
    """
    """
    query_data = extract_sequence(query)
    database_data = extract_sequence(database, listhomology)
    idcompair = 0
    comparison = []
    for query in query_data:
        for database in database_data:
            if(query[0] != database[0]):
                # query query_seq database database_seq output_file path_clustalo
                multifasta_file = "{0}{1}{2}.fasta".format(results, os.sep,
                                                      idcompair)
                alignment_file =  "{0}{1}clustalo_{2}.fasta".format(results,
                                                                    os.sep,
                                                                    idcompair)
                comparison += [query + database +
                               [replace_sequence(cmd, path_clustalo,
                                    multifasta_file, alignment_file),
                                multifasta_file, alignment_file]]
                idcompair += 1
    return comparison


def extract_data(alignment_file):
    """
    """
    data_aln = {}
    keys = []
    try:
        with open(alignment_file) as align:
            for line in align:
                if line.startswith(">"):
                    head = line[1:].replace("\n", "").replace("\r", "")
                    data_aln[head] = ""
                    keys += [head]
                elif len(line) > 0 :
                    data_aln[head] += line.replace("\n", "").replace("\r", "")
            assert(data_aln != {} and len(data_aln) == 2)
    except IOError:
        sys.exit("Error cannot open {0}".format(alignment_file))
    except AssertionError:
        sys.exit("Nothing extracted or illegal length (!=2 alignment) "
                 "from {0}".format(alignment_file))
    if len(data_aln[keys[0]]) != len(data_aln[keys[1]]):
        sys.exit("The length of the alignment are different :\n"
                 "{0[0]}:{1} and {0[1]}:{2}".format(keys,
                                                    len(data_aln[keys[0]]),
                                                    len(data_aln[keys[1]])))
    return data_aln


def estimate_parameters(query, database, seq_template, seq_aln):
    """ Compute identity, similarity and coverage
    """
    similarAA = ['AG', 'AP', 'AS', 'AT', 'DE', 'DN', 'DQ', 'ED', 'EN', 'EQ',
                'FW', 'FY', 'GA', 'GP', 'GS', 'GT', 'HK', 'HR', 'IL', 'IM',
                'IV', 'KH', 'KR', 'LI', 'LM', 'LV', 'MI', 'ML', 'MV', 'ND',
                'NE', 'NQ', 'PA', 'PG', 'PS', 'PT', 'QD', 'QE', 'QN', 'RH',
                'RK', 'SA', 'SG', 'SP', 'ST', 'TA', 'TG', 'TP', 'TS', 'VI',
                'VL', 'VM', 'WF', 'WY', 'YF', 'YW']
    aligned = 0.0
    id_aa = 0.0
    similar = 0.0
    s1 = seq_template.replace("-", "")
    s2 = seq_aln.replace("-", "")
    for i in xrange(len(seq_template)):
        if(seq_template[i] == "-" or seq_aln[i] == "-"):
            pass
        elif(seq_template[i] == seq_aln[i] and seq_template[i] is not "X"
             and seq_aln[i] != "X"):
            id_aa += 1.0
            aligned += 1.0
            similar += 1.0
        elif(get_element(similarAA, seq_template[i] + seq_aln[i])):
            similar += 1.0
            aligned += 1.0
        elif(seq_template[i] != "X" and seq_aln[i] != "X"):
            aligned += 1.0
    # identity = 100.0*count/(float(len(aln1.translate(None,"-."))
    # +len(aln2.translate(None,"-.")))/2.0)
    similarity = 100.0 * similar / float(min(len(s1), len(s2)))
    identity = 100.0 * id_aa / float(min(len(s1), len(s2)))
    coverage = 100.0 * aligned / float(max(len(s1), len(s2)))
    return [query, database, len(s1), len(s2), round(identity, 2),
            round(similarity,2), round(coverage, 2)]
    #return {"query":query, "database":database, "identity":identity,
            #"similarity":similarity, "coverage":coverage,"query_len": len(s1),
            #"database_length":len(s2)}


def align_sequences(compair):
    """Thread action
    """
    # Write multifasta file
    create_fasta(compair)
    # Run alignment
    run_command(compair[4])
    # Extract aligment
    data_alignment = extract_data(compair[6])
    # Calculate parameters
    # query_seq database_seq
    result_ali = estimate_parameters(compair[0], compair[2],
                                     data_alignment[compair[0]],
                                     data_alignment[compair[2]])
    # Remove files
    remove_temp_files(compair[5:7])
    return result_ali


def compute_sequence_properties(results, conf_data, path_clustalo, query,
                                database, listhomology, thread):
    """Compute sequence identity
    """
    print("Get comparison")
    # generator not available because multiprocessing pickling
    comparaison = map_sequence(results, conf_data.hdict["clustalo"],
                               query, database, listhomology, path_clustalo)
    print("Extraction done. {0} comparisons to do.".format(len(comparaison)))
    startTime = time.time()
    pool = mp.Pool(processes=thread)
    asyncResult = pool.map_async(align_sequences, comparaison)
    result_list = asyncResult.get()
    #result_list = []
    #for comp in comparaison:
        #result_list += [align_sequences(comp)]
    print("Comparison done. Time {0}s.".format(time.time() - startTime))
    return result_list


def get_unique_element(data_list):
    """Get the unique elements in list
    """
    return {}.fromkeys(data_list).keys()

def difference(lst1, lst2):
    return list(set(lst1)^set(lst2))


def get_sequence_length(fasta_file, listhomology):
    """Get sequence length of homologous
    """
    sequence_length = {}
    interest = False
    header = ""
    try:
        with open(fasta_file, "rt") as fasta:
            for line in fasta:
                if line.startswith(">"):
                    if interest:
                        sequence_length[header] = len(sequence)
                    header = line[1:].replace("\n", "").replace("\r", "").split(" ")[0]
                    if " " in header:
                        header = header.split(" ")[0]
                    # if element in list of homologous
                    interest = get_element(listhomology, header)
                    sequence = ""
                elif interest and len(line) > 0:
                    sequence += line.replace("\n", "").replace("\r", "")
            if interest:
                sequence_length[header] = len(sequence)
                interest = False
            assert(len(sequence_length.keys()) == len(listhomology))
    except IOError:
        sys.exit("Error cannot open {0}".format(fasta_file))
    except AssertionError:
        sys.exit("The length of every blast hit has not "
                 "been found in {0}.\nMissing elements {1}"
                 .format(fasta_file, difference(sequence_length.keys(), listhomology)))
    return sequence_length


def write_hit_length(results, selected, not_selected, sequence_length,
                     output_file=None):
    """Write hit length
    """
    if not output_file:
        output_file =  results + os.sep + soft + "_hit_length.tsv"
    try:
        with open(output_file, "wt") as output:
            output_writer = csv.writer(output,delimiter="\t")
            output_writer.writerow(["Hit", "Length", "Selected"])
            for element in selected:
                output_writer.writerow([element, sequence_length[element],
                                        1])
            for element in not_selected:
                output_writer.writerow([element, sequence_length[element],
                                        0])
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


def filter_sequence_length(results, fasta_file, listhomology, length_min,
                           length_max, soft):
    """
    """
    selected = []
    # Calculate hit length
    sequence_length = get_sequence_length(fasta_file, listhomology)
    # Filter hits
    for element in listhomology:
        select_min = True
        select_max = True
        if length_min:
            select_min = (sequence_length[element] >= length_min)
        if length_max:
            select_max = (sequence_length[element] <= length_max)
        if select_min and select_max:
            selected.append(element)
    # Write hit homology
    not_selected = set(listhomology) - set(selected)
    output_file =  results + os.sep + soft + "_hit_length.tsv"
    print("Write hit sequence length in {0}".format(output_file))
    write_hit_length(results, selected, not_selected, sequence_length,
                     output_file)
    return selected


def write_check_data(results, result_dict, origin, nbest, selection=None):
    """Write the different properties of selected homologues
    """
    output_file = results + os.sep + origin + "_hit_properties.tsv"
    try:
        with open(output_file, "wt") as output:
            output_writer = csv.writer(output, delimiter='\t')
            header = ["Reference", "Query", "Query_length", "Reference_length",
                      "Identity", "Similarity", "Coverage"]
            if selection:
                header += ["Selected"]
            # Write header
            output_writer.writerow(header)
            for query in result_dict:
                if selection:
                    if query in selection:
                        subset_data = [[query] + hit + [1]
                                       for hit in result_dict[query]]
                    else:
                        subset_data = [[query] + hit + [0]
                                       for hit in result_dict[query]]
                else:
                    subset_data = [[query] + hit for hit in result_dict[query]]
                # Sort depending on the identity and coverage
                subset_data.sort(key=lambda x: x[4] + x[6], reverse=True)
                if nbest > 0:
                    short_set = subset_data[0: nbest]
                else:
                    short_set = subset_data
                for hit in short_set:
                    output_writer.writerows(short_set)
    except IOError:
        sys.exit("Error : cannot open file {0}".format(output))


def order_dict(result_stat):
    """Order the information
    """
    result_dict = {}
    for hit in result_stat:
        if hit[1] in result_dict:
            result_dict[hit[1]] += [[hit[0]] + hit[2:]]
        else:
            result_dict[hit[1]] = [[hit[0]] + hit[2:]]
    return result_dict


def check_parameters(result_dict, element, idmin, idmax, simmin, simmax, covmin,
                     covmax):
    for hit in result_dict[element]:
        select_id_min = True
        select_id_max = True
        select_sim_min = True
        select_sim_max = True
        select_cov_min = True
        select_cov_max = True
        if idmin:
            select_id_min = (hit[3] >= idmin)
        if idmax:
            select_id_max = (hit[3] <= idmax)
        if idmin:
            select_cov_min = (hit[5] >= covmin)
        if idmax:
            select_cov_max = (hit[5] <= covmax)
        if(select_id_min and select_id_max
            and select_cov_min and select_cov_max):
            return True
    return False


def filter_sequence_align(listhomology, result_dict, idmin, idmax, simmin,
                           simmax, covmin, covmax):
    selected = []
    for element in listhomology:
        # We could do better with sort, but latter maybe
        if check_parameters(result_dict, element, idmin, idmax, simmin,
                           simmax, covmin, covmax):
            selected.append(element)
    return selected


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    filtered_elements = None
    all_listhomology = []
    all_filtered_homologous = []
    # Load arguments
    args = getArguments()
    # Check software
    if args.path_software:
        get_path = extract_path
    else:
        get_path = path_empty
    # Configure option
    conf_data = Inconfig(args.config, args.results)
    # Analyze the homology for each software
    for num, soft in enumerate(args.software):
        output = args.results + os.sep + soft + "_protein_homology.txt"
        path_soft = get_path(soft, args.path_software, num)
        # Compute homology
        if not os.path.isfile(output) or "force_computation" in args.behavior:
            #UGLY all of this
            try:
                if soft == "hmmsearch" or soft == "hmmscan":
                    if not args.query_hmm:
                        sys.exit("Please provide the hmm db with -dm option")
                    #if args.reverse:
                        ## Reverse query and database
                    run_command(replace_motif(conf_data.hdict[soft], path_soft,
                                        args.database, args.query,
                                        args.query_hmm, output, args.thread,
                                        args.e_value, args.tmp))
                    #else:
                    #run_command(replace_motif(conf_data.hdict[soft],
                                                #path_soft, args.query,
                                                #args.database,
                                                #args.query_hmm, output,
                                                #args.thread, args.e_value))
                else:
                    if((soft == "blastp" or soft == "psiblast"
                        or soft == "blastx" or soft == "tblastn"
                        or soft == "mmseqs")
                       and args.database_blast):
                        run_command(
                            replace_motif(conf_data.hdict[soft], path_soft,
                                          args.query, args.database_blast,
                                          args.query_hmm, output,
                                          args.thread, args.e_value, args.tmp))
                    else:
                        run_command(
                            replace_motif(conf_data.hdict[soft], path_soft,
                                          args.query, args.database,
                                          args.query_hmm, output,
                                          args.thread, args.e_value, args.tmp))
            except KeyError:
                sys.exit("Key {0} is not in the configuration "
                         "file.".format(soft))
        # Get the complete list here
        # Parse result
        print("parse result : " + soft)
        listhomology = extract_homology(output, soft) #, args.reverse)
        # Uniquify elements
        listhomology = get_unique(listhomology)
        # Sort the list for bisect
        listhomology.sort()
        if not "cumulative" in args.behavior and len(listhomology) > 0:
            # Filter candidates with length
            if args.filter_length_minimum or args.filter_length_maximum:
                print("Filter length")
                #if args.reverse:
                    #listhomology = filter_sequence_length(args.results,
                                        #args.query, listhomology,
                                        #args.filter_length_minimum,
                                        #args.filter_length_maximum, soft)
                #else:
                listhomology = filter_sequence_length(args.results,
                                args.database, listhomology,
                                args.filter_length_minimum,
                                args.filter_length_maximum, soft)
            # Check identity, similarity, coverage
            if "check" in args.behavior:
                # No output_fasta
                #if args.reverse:
                #else:
                result_stat = compute_sequence_properties(
                                    args.results, conf_data, args.path_clustalo,
                                    args.query, args.database, listhomology,
                                    args.thread)
                # Order the information
                result_dict = order_dict(result_stat)
                if(args.filter_identity_minimum or args.filter_identity_maximum
                or args.filter_similarity_minimum or args.filter_similarity_maximum
                or args.filter_coverage_minimum or args.filter_coverage_maximum):
                    print("Filter sequence based on identity, similarity or coverage")
                    listhomology = filter_sequence_align(
                                        listhomology, result_dict,
                                        args.filter_identity_minimum,
                                        args.filter_identity_maximum,
                                        args.filter_similarity_minimum,
                                        args.filter_similarity_maximum,
                                        args.filter_coverage_minimum,
                                        args.filter_coverage_maximum)
                    write_check_data(args.results, result_dict, soft, args.nbest,
                                    listhomology)
                else:
                    write_check_data(args.results, result_dict, soft, args.nbest)
                # write statistics
                # write_data(args.results, result_stat, soft)
            # Extract homologous sequence
            if "extract" in args.behavior:
                output_fasta = args.results + os.sep + soft + "_protein_homology.fasta"
                print(listhomology)
                # Get a multifasta that contains only the candidates and
                # homologous proteins
                print("Write sequence of the hits in {0}".format(output_fasta))
                register_fasta(args.results, args.database, listhomology, output_fasta)
        else:
            # Add elements
            all_listhomology += [listhomology]
    # TODO Venn diagram here
    # Venn diagram
    if"cumulative" in args.behavior and all_listhomology:
        group_elements = []
        for list_element in all_listhomology:
            group_elements += list_element
        # Unique elements for bisect
        group_elements = get_unique(group_elements)
        # Sort the list for bisect search
        group_elements.sort()
        # Filter candidates with length
        if((args.filter_length_minimum or args.filter_length_maximum)
           and len(group_elements) > 0):
            print("Filter length")
            #if args.reverse:
                #group_elements = filter_sequence_length(args.results, args.query,
                                        #group_elements,
                                        #args.filter_length_minimum,
                                        #args.filter_length_maximum,
                                        #"all")
            #else:
            group_elements = filter_sequence_length(args.results, args.database,
                                                group_elements,
                                                args.filter_length_minimum,
                                                args.filter_length_maximum,
                                                "all")
        if "check" in args.behavior and len(group_elements) > 0:
            # Analyse candidates
            #if args.reverse:
            result_stat = compute_sequence_properties(args.results, conf_data,
                                                    args.path_clustalo,
                                                    args.query, args.database,
                                                    group_elements,
                                                    args.thread)
            #else:
            #result_stat = compute_sequence_properties(args.results, conf_data,
                                                    #args.path_clustalo,
                                                    #args.database, args.query,
                                                    #group_elements,
                                                    #args.thread)
            # Order the information
            result_dict = order_dict(result_stat)
            if(args.filter_identity_minimum or args.filter_identity_maximum
               or args.filter_similarity_minimum or args.filter_similarity_maximum
               or args.filter_coverage_minimum or args.filter_coverage_maximum):
                print("Filter sequence based on identity, similarity or coverage")
                group_elements = filter_sequence_align(
                                    group_elements, result_dict,
                                    args.filter_identity_minimum,
                                    args.filter_identity_maximum,
                                    args.filter_similarity_minimum,
                                    args.filter_similarity_maximum,
                                    args.filter_coverage_minimum,
                                    args.filter_coverage_maximum)
                if len(group_elements) == 0:
                    print("Warning : No element selected with this constraints "
                          "of identity, similarity and coverage !")
                write_check_data(args.results, result_dict, "all", args.nbest,
                                 group_elements)
            else:
                write_check_data(args.results, result_dict, "all", args.nbest)
        if "extract" in args.behavior and len(group_elements) > 0:
            output_fasta = args.results + os.sep + "all_protein_homology.fasta"
            # Get a multifasta that contains only the homologous proteins
            print("Write sequence of the hits in {0}".format(output_fasta))
            #if args.reverse:
                #register_fasta(args.results, args.query, group_elements,
                               #output_fasta)
            #else:
            register_fasta(args.results, args.database, group_elements,
                            output_fasta)

if __name__ == '__main__':
    main()
