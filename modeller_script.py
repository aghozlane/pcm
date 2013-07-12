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

"""Run modeling with modeller."""


from __future__ import print_function
import argparse
import ConfigParser
import os
import sys
import urllib2
import subprocess
import re
import csv
import math
import textwrap
import multiprocessing as mp
try:
    import requests
    REQUESTS = True
except ImportError:
    REQUESTS = False
    print("Could not import requests.{0}No structure checking with prosa or "
          "ProQ or Verify3D will be available".format(os.linesep),
          file=sys.stderr)
try:
    from modeller import *
    from modeller.parallel import *
    from modeller.automodel import *
    from modeller.scripts import complete_pdb
    MODELLER = True
except ImportError:
    MODELLER = False
    print("Could not import modeller.{0}No modeling "
          "will be available".format(os.linesep), file=sys.stderr)
try:
    from MyModel import RestraintModel
except ImportError:
    sys.exit("modeller_script need to import MyModel.py, leave it in the same "
             "directory as modeller_script.py")
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    MATPLOTLIB = True
except ImportError:
    MATPLOTLIB = False
    print("Could not import matplotlib.{0}The statistic measures "
          "won't be drawn.".format(os.linesep), file=sys.stderr)


__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2013, University of Paris VII"
__credits__ = ["Amine Ghozlane", "Joseph Rebehmed", "Alexandre G. de Brevern"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@inserm.fr"
__status__ = "Developpement"


class ModelingConfig:
    """Configure alignment.
    """

    def __init__(self, config_file, results):
        """Instantiate ModelingConfig object
           :Parameters:
            config_file: Configuration file path
            results: Result path
        """
        self.hdict = {}
        self.phyconfig_file = '{0}modelingconfig.cfg'.format(results)
        self.config = ConfigParser.RawConfigParser()
        if(config_file is not None):
            self.phyconfig_file = config_file
            self.readconfig()
        elif(os.path.isfile(self.phyconfig_file)):
            self.readconfig()
        else:
            self.writeconfig()
            self.readconfig()

    def readconfig(self):
        """Read config data
        """
        # If config parser empty
        if(not os.path.isfile(self.phyconfig_file)):
            self.writeconfig()
        # Read config file
        self.config.read(self.phyconfig_file)
        # Get parameter value
        self.hdict["download_url"] = self.config.get('PDB_config',
                                                     'download_url')
        self.hdict["hhsearch"] = self.config.get('PDB_config', 'hhsearch')
        self.hdict["psiblast"] = self.config.get('PDB_config', 'psiblast')
        self.hdict["jackhmmer"] = self.config.get('PDB_config', 'jackhmmer')
        self.hdict["clustalo"] = self.config.get('Alignment_config',
                                                 'clustalo')
        self.hdict["clustalw2"] = self.config.get('Alignment_config',
                                                  'clustalw2')
        self.hdict["mafft"] = self.config.get('Alignment_config',
                                              'mafft')
        self.hdict["muscle"] = self.config.get('Alignment_config',
                                               'muscle')
        self.hdict["t_coffee"] = self.config.get('Alignment_config',
                                                 't_coffee')
        self.hdict["psipred"] = self.config.get('Secondary_structure_config',
                                                'psipred')
        self.hdict["procheck"] = self.config.get('Checking_config',
                                                 'procheck')
        self.hdict["prosa"] = self.config.get('Checking_config', 'prosa')
        self.hdict["proq"] = self.config.get('Checking_config', 'proq')
        self.hdict["verify3D"] = self.config.get('Checking_config',
                                                 'verify3D')

    def writeconfig(self):
        """Write modeling config
        """
        self.config.add_section('PDB_config')
#         self.config.set("PDB_config", "hhmake",
#                         "%path_softhhsearch -i %multifasta  -d %database "
#                         "-cpu %proc -o %output")
        self.config.set("PDB_config", "download_url",
                        "http://www.rcsb.org/pdb/files/")
        self.config.set("PDB_config", "hhsearch", "%path_softhhsearch "
                        "-i %multifasta  -d %database -cpu %proc -o %output")
        self.config.set("PDB_config", "psiblast", "%path_softpsiblast "
                        "-query %multifasta -db %database -out %output "
                        "-evalue 0.001 -outfmt 6 -num_threads %proc")
        self.config.set("PDB_config", "jackhmmer", "%path_softjackhmmer "
                        "--cpu %proc  -o %output %multifasta %database")
        self.config.add_section('Alignment_config')
        self.config.set('Alignment_config', 'clustalo',
                        "%path_softclustalo -i %multifasta -o %output "
                        "--threads=%proc --auto -t Protein --outfmt=fa")
        self.config.set('Alignment_config', 'clustalw2',
                        "%path_softclustalw2 -INFILE=%multifasta "
                        "-OUTPUT=PIR -OUTFILE=%output")
        self.config.set('Alignment_config', 'mafft',
                        "%path_softmafft --auto --thread %proc "
                        "%multifasta > %output")
        self.config.set('Alignment_config', 'muscle', "%path_softmuscle "
                        "-in %multifasta -out %output")
        self.config.set('Alignment_config', 't_coffee', "%path_softt_coffee "
                        "%multifasta -outfile %output -output pir_aln "
                        "-n_core %proc")
        self.config.add_section('Secondary_structure_config')
        self.config.set('Secondary_structure_config', 'psipred',
                        "%path_softrunpsipred %fasta ")
        # -mode 3dcoffee -template_file %pdb
        self.config.add_section('Checking_config')
        self.config.set('Checking_config', 'procheck',
                        "%path_softprocheck.scr %pdb 1.5")
        self.config.set('Checking_config', 'proq',
                        "http://www.sbc.su.se/~bjornw/ProQ/ProQ.cgi")
        self.config.set('Checking_config', 'prosa',
                        "https://prosa.services.came.sbg.ac.at/")
        self.config.set('Checking_config', 'verify3D',
                        "http://nihserver.mbi.ucla.edu/Verify_3D/")
        # Write configuration data
        try:
            # Writing our configuration file to 'example.cfg'
            with open(self.phyconfig_file, 'wt') as configfile:
                self.config.write(configfile)
        except IOError:
            sys.exit("Error : cannot open "
                     "file {0}".format(self.phyconfig_file))


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


def isdir(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file.".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def islimit(value):
    """Check if value is in confidence limit.
      :Parameters:
          conf_value: Confidence value
    """
    try:
        conf_value = int(value)
        if conf_value < 0 or conf_value > 9:
            raise argparse.ArgumentTypeError("Psipred confidence threshold "
                                             "must be between 0-9")
    except ValueError:
        raise argparse.ArgumentTypeError("Value \"{0}\" is not"
                                         " an integer".format(value))
    return conf_value


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    local_path = ".{0}".format(os.sep)
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    parser.add_argument("-l", dest="list_operations",
                        default=["model", "profile", "check"], type=str,
                        nargs='+', choices=["model", "profile", "check"],
                        help="Select the operations : model and/or profile "
                        "and/or structure checking (default : model, profile, "
                        "and check are done)")
    parser.add_argument('-f', dest='multifasta_file', type=isfile,
                        help='Multifasta file with model and '
                             'template sequences.')
    parser.add_argument('-i', dest='alignment_file', type=isfile,
                        help='Alignment file in PIR format.')
    parser.add_argument('-a', dest='alignment_software', type=str,
                        default="clustalo",
                        choices=["clustalw2", "clustalo", "mafft", "muscle",
                                 "t_coffee"],
                        help="Indicates the software that should be used to "
                        "align sequences (Default = clustalo).")
    parser.add_argument('-g', dest='path_alignment', type=isdir, default=None,
                        help='Path to the alignment software.')
    parser.add_argument('-p', dest='pdb', type=str, nargs='+',
                        help="List of pdb files or codes to use as template.")
    parser.add_argument('-pi', dest='pdb_identification', type=str,
                        nargs='+', default=None,
                        choices=["hhsearch", "psiblast", "jackhmmer"],
                        help="Indicate the software to search homologous.")
    parser.add_argument('-pd', dest='pdb_identification_database',
                        type=isfile, nargs='+', default=[],
                        help="Indicate the support database for template "
                        "research.")
    parser.add_argument('-pp', dest='pdb_identification_path', type=isdir,
                        nargs='+', default=[], help="Path to the software for "
                        "identification.")
    parser.add_argument('-ps', dest='pdb_identification_strategy', type=str,
                        default="best", choices=["best"],  # , "covering", "multiple"
                        help="Select the strategy  to the software for "
                        "identification.")
    parser.add_argument('-e', dest='model_name', type=str,
                        help='Code of the sequence to modelize (required if '
                        'several sequence with no structure associated are '
                        'present in the multifasta).')
    parser.add_argument('-n', dest='number_model', type=int, default=8,
                        help='Number of model to produce (default = 8).')
    parser.add_argument('-q', dest='model_quality', type=str, default="fast",
                        choices=["very_fast", "fast", "normal", "max"],
                        help='Adjust the quality of the modeling (default = '
                        'fast).')
    parser.add_argument('-ht', dest='add_heteroatom', type=int, default=0,
                        help="Indicates the number of hetero-atom residue(s) "
                        "that should be added to the alignment software.")
    parser.add_argument('-hm', dest='heteroatom_models', default=[],
                        type=str, nargs='+',
                        help="Indicate the models for which hetero-atom "
                        "residue(s) should be added to the alignment "
                        "software (choose model + templates).")
    parser.add_argument('-j', dest='path_psipred', type=isdir, default=None,
                        help='Path to Psipred software.')
    parser.add_argument('-d', dest='psipred', type=isfile, default=None,
                        help='Psipred file used for modeling and checking the '
                        'model structures (*.psipass2).')
    parser.add_argument('-b', dest='limit_confidence', type=islimit,
                        default=7, help='Confidence limit for psipred'
                        '(0-9 - default = >7).')
    parser.add_argument('-s', dest='structure_check', type=str, nargs='+',
                        choices=["procheck", "proq", "prosa", "verify3D"],
                        default=["procheck", "proq", "prosa", "verify3D"],
                        help='Select software for structure checking '
                        '(ProQ significance is enhanced with psipred '
                        'results).')
    parser.add_argument('-k', dest='path_check', type=isdir, default=None,
                        help='Indicate the path to procheck software.')
    parser.add_argument('-sb', dest='number_best', type=int, default=5,
                        help='Select number of models for verification'
                        '(from the best model according to dope score - '
                        'default = 5).')
    parser.add_argument('-sm', dest='modeller_summary', type=isfile,
                        default=None, help="Indicate modeller script file "
                        "(instead of alignment file) for checking.")
    parser.add_argument('-r', dest='results', type=isdir, default=local_path,
                        help='Path to result directory. (Default = current '
                        'directory is preferred because of modeller '
                        'function).')
    parser.add_argument('-da', dest='disable_autocorrect',
                        action='store_true', default=False,
                        help='Disable the autocorrect of the multifasta file.')
    parser.add_argument('-m', dest='list_atom', type=isfile, default=None,
                        help='List of interest atom')
    parser.add_argument('-md', dest='default_distances', type=isfile,
                        default=None, help='Distance file')
    parser.add_argument('-t', dest='thread', default=mp.cpu_count(),
                        type=int, help="Number of thread (Default = all cpus "
                        "available will be used).")
    parser.add_argument('-c', dest='config', type=isfile, default=None,
                        help='Path to configuration file.')
    return parser.parse_args(), parser


def extract_elements(template_search_file, regex_text, order):
    """
    """
    elements = []
    regex = re.compile(regex_text)
    nb_elements = 0
    try:
        with open(template_search_file, "rt") as template_search:
            for line in template_search:
                match = regex.match(line)
                if match:
                    # Get name and e-value
                    elements += [[match.group(order[0]),
                                  float(match.group(order[1]))]]
                    # recover
                    # int(match.group(order[2]))
                    nb_elements += 1
                # Consider only the 10 best templates
                if nb_elements > 10:
                    break
        assert(len(elements) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(template_search_file))
    except ValueError:
        sys.exit("Error cannot convert to interger "
                 "{0}".format(match.group(order[1])))
    except AssertionError:
        sys.exit("No template have been identified in "
                 "{0}".format(template_search_file))
    return elements


def identify_template(conf_data, multifasta_file, thread, pdb_identification,
                      pdb_identification_database, pdb_identification_path,
                      strategy, results):
    """
    """
    PDB = []
    elements = []
    if(len(pdb_identification_path) == len(pdb_identification)):
        list_path = pdb_identification_path
    else:
        list_path = [""] * len(pdb_identification)
    if(len(pdb_identification_database) != len(pdb_identification)):
        sys.exit("Please specify explicitly for each software a database "
                 "to explore")
    for i in xrange(len(pdb_identification)):
        output = results + pdb_identification[i] + "_output.txt"
        run_command(replace_motif(conf_data.hdict[pdb_identification[i]],
                                  list_path[i], multifasta_file,
                                  "", output, thread, "", "",
                                  pdb_identification_database[i]))
        # Get elements
        if(pdb_identification[i] == "hhsearch"):
            elements = extract_elements(output, "\s*[0-9]+\s+(\w+)\_[A-Z].+"
                                        "\s+[0-9\.]+\s+(\S+)\s+" + "\S+\s+" * 6
                                        + "\([0-9]+\)",
                                        [1, 2])
            # "\(([0-9]+)\)"
        elif(pdb_identification[i] == "psiblast"):
            elements = extract_elements(output,
                                        "\w+\s+.+\|(\w+)\|[A-Z]\s+.+\s+(\S+)"
                                        "\s+[0-9.]+{0}".format(os.linesep),
                                        [1, 2])
        elif(pdb_identification[i] == "jackhmmer"):
            elements = extract_elements(output,
                                        "^\+\s+(\S+)\s+.+\|(\w+)\|[A-Z]",
                                        [2, 1])
        # Select PDB
        if strategy == "best" and len(elements) > 0:
            PDB += [elements[0][0]]
    nb_pdb = len(PDB)
    if nb_pdb > 1:
        PDB = get_unique_no_order(PDB)
    elif nb_pdb == 0:
        sys.exit("No template structure found !!")
    return PDB


def download_pdb(conf_data, pdb, results):
    """Download PDB file
     :Parameters:
        conf_data: Configuration dictionary
        pdb: Name of the pdb file
        results: Output directory
    """
    # Open our local file for writing
    outfilename = results + pdb + ".pdb"
    try:
        f = urllib2.urlopen(conf_data.hdict["download_url"] + pdb + ".pdb")
        with open(outfilename, "wt") as local_file:
            local_file.write(f.read())
    except urllib2.HTTPError, e:
        print ("HTTP Error:", e.code, pdb, file=sys.stderr)
    except urllib2.URLError, e:
        print ("URL Error:", e.reason, pdb, file=sys.stderr)
    except IOError:
        print("Something went wrong with {0}".format(outfilename),
              file=sys.stderr)
    return outfilename


def get_pdb(conf_data, pdb_list, results):
    """Check pdb extension
     :Parameters:
        conf_data: Configuration dictionary
        pdb_list: List of PDB
        results: Output directory
    """
    pdb_codes = []
    pdb_files = []
    for pdb in pdb_list:
        # Correspond to a PDB file
        if not pdb.endswith('.pdb') or not os.path.isfile(pdb):
            pdb_codes += [os.path.basename(pdb).split(".")[0]]
            pdb_files += [download_pdb(conf_data, pdb, results)]
        # Download corresponding pdb
        else:
            pdb_codes += [".".join(os.path.basename(pdb).split(".")[:-1])]
            pdb_files += [pdb]
    return pdb_codes, pdb_files


def run_command(cmd):
    """Run command
      :Parameters:
          cmd: Command to run
    """
    try:
        # execute the command line
        retcode = subprocess.call(cmd, shell=True)
        # Case of no return
        if retcode is None:
            sys.exit("Child was terminated")
    except OSError as error:
        sys.exit("Execution failed: {0}".format(error))
    except:
        sys.exit("There is something wrong with the command: {0}".format(cmd))


def replace_motif(build_command, path_soft, multifasta_file, pdb_files,
                  output, thread, psipred, fasta, database):
    """Set the software command
     :Parameters:
         build_command: Command
         path_soft: Path to the alignment software
         multifasta_file: Path to the multifasta file
         pdb_files: List of PDB files
         output: Output path
         thread: Number of thread
     Returns:
      Command to execute.
    """
    print(build_command, file=sys.stderr)
    if path_soft:
        build_command = build_command.replace('%path_soft', path_soft)
    else:
        build_command = build_command.replace('%path_soft', "")

    if thread:
        build_command = build_command.replace('%proc', str(thread))
    else:
        build_command = build_command.replace('%proc', str(mp.cpu_count()))
    build_command = build_command.replace('%multifasta', multifasta_file)
    try:
        if(len(pdb_files) > 0):
            build_command = build_command.replace('%pdb', " ".join(pdb_files))
    except TypeError:
        pass
    build_command = build_command.replace('%output', output)
    build_command = build_command.replace('%psipred', psipred)
    build_command = build_command.replace('%fasta', fasta)
    build_command = build_command.replace('%database', database)
    print(build_command, file=sys.stderr)
    return build_command


def check_format(aln_data, pdb_codes):
    """Check the pir file format
     :Parameters:
      aln_data: Alignment read file
      pdb_codes: List of PDB
     Returns:
      Status of the format and extraction result
    """
    ide = ""
    status = True
    count = 0
    ens = ["P1", "F1", "DL", "DC", "RL", "RC", "XX"]
    data = {}
    pdb_elements = pdb_codes[:]
    regex_suite = (re.compile(r"^>" + "|".join([i + ";([\w-]+)"for i in ens])),
                   re.compile(r"^([sequence|structureX]:[\w-]+" + ":\S+" * 8
                              + ")"),
                   re.compile(r"([\w*-]+{0})".format(os.linesep)))
    for line in aln_data:
        # Match title
        if count in (0, 2):
            match = regex_suite[0].match(line)
            if match:
                ide = match.group(1)
                data[ide] = [None, ""]
                if ide in pdb_elements:
                    pdb_elements.pop(pdb_elements.index(match.group(1)))
                count = 1
        # Match subtitle
        elif count == 1:
            match = regex_suite[count].match(line)
            if match:
                data[ide][0] = match.group(1)
                count += 1
            else:
                print("Warning: next line does not correspond to subtitle "
                      "format required for modeller"
                      " :{0}{1}".format(os.linesep, line), file=sys.stderr)
                count += 1
                status = False
        # Match Sequence
        if count == 2:
            match = regex_suite[count].match(line)
            if not match:
                print("Warning: next line does not correspond to a correct "
                      "protein sequence :{0}{1}".format(os.linesep, line),
                      file=sys.stderr)
                status = False
            else:
                data[ide][1] += match.group(1)
    assert(len(pdb_elements) == 0)
    return status, data


def adjust_PIR_format(aln_pir_file, data_dict, pdb_codes, pdb_files,
                      add_heteroatom, heteroatom_models):
    """Correct pir format and add heteroatom in the alignment
     :Parameters:
      aln_pir_file: Path to the PIR file
      data_dict: Dictionary containing alignment data
      pdb_codes: List of PDB
      pdb_files: List of PDB files
      add_heteroatom: Number of heteroatom
      heteroatom_models: List of model for which heteroatom should be added
    """
    output_file = (os.path.dirname(aln_pir_file) + os.sep +
                   os.path.basename(aln_pir_file).split(".")[0]
                   + "_corrected.pir")
    het_atm_flag = True
    het_atm = ""
    structure_present = False
    if add_heteroatom > 0:
        het_atm = "." * add_heteroatom
    try:
        with open(output_file, "wt") as aln_pir:
            for element in data_dict:
#                 pdb_start_posit = 1
                if(element not in heteroatom_models
                   and len(heteroatom_models) > 0):
                    het_atm_flag = False
                aln_pir.write(">P1;{0}{1}".format(element, os.linesep))
                if data_dict[element][0] and not add_heteroatom:
                    aln_pir.write(data_dict[element][0])
                else:
                    # Write midline characteristic
                    type_data = "sequence"
                    if element in pdb_codes:
                        type_data = "structureX"
                        structure_present = True
                    sequence = data_dict[element][1].replace("-", "")
                    sequence = sequence.replace(os.linesep, "")
                    sequence = sequence.replace("*", "")
                    if het_atm_flag:
                        aln_pir.write("{0}:{1}: :A: :A"
                                      .format(type_data, element)
                                      + ": " * 4 + os.linesep)
                    else:
                        aln_pir.write("{0}:{1}:  :A: :A"
                                      .format(type_data, element)
                                      + ": " * 4 + os.linesep)
                end_aln_posit = data_dict[element][1].index("*")
                if het_atm_flag:
                    aln_pir.write("{0}{1}{2}{3}".format(
                        data_dict[element][1][0:end_aln_posit],
                        het_atm, data_dict[element][1][end_aln_posit:],
                        os.linesep))
                else:
                    aln_pir.write("{0}{1}{2}".format(
                        data_dict[element][1][0:end_aln_posit],
                        data_dict[element][1][end_aln_posit:], os.linesep))
                het_atm_flag = True
    except IOError:
        sys.exit("Error cannot open {0}".format())
    assert(structure_present)
    return output_file


def check_pir(aln_pir_file, pdb_codes, pdb_files, add_heteroatom,
              heteroatom_models):
    """Run checking of PIR alignment file
     :Parameters:
      aln_pir_file: Path to the PIR file
      pdb_codes: List of PDB
      add_heteroatom: Number of heteroatom
      heteroatom_models: List of model for which heteroatom should be added
    """
    try:
        if not aln_pir_file.endswith('.pir'):
            print("Warning : The alignment file is expected "
                  "to be in PIR format", file=sys.stderr)
        with open(aln_pir_file, "rt") as aln_pir:
            aln_data = aln_pir.readlines()
        status, data_dict = check_format(aln_data, pdb_codes)
        if not status or add_heteroatom > 0:
            print("Try to correct PIR alignment format.", file=sys.stderr)
            aln_pir_file = adjust_PIR_format(aln_pir_file, data_dict,
                                             pdb_codes, pdb_files,
                                             add_heteroatom, heteroatom_models)
    except AssertionError:
        sys.exit("All PDB structures must be referenced in the alignment.")
    except IOError:
        sys.exit("Error cannot open {0}".format(aln_pir_file))
    return aln_pir_file


def get_fasta_data(aln_fasta_file):
    """Extract fasta data
     :Parameters:
       aln_fasta_file: Path to the fasta file
    """
    head = ""
    data_fasta = {}
    regex_head = re.compile(r"^>([\w-]+)")
    regex_sequence = re.compile(r"([\w*-]+{0})".format(os.linesep))
    try:
        with open(aln_fasta_file, "rt") as aln_fasta:
            for i in aln_fasta:
                match_head = regex_head.match(i)
                match_sequence = regex_sequence.match(i)
                if match_head:
                    head = match_head.group(1)
                    data_fasta[head] = ""
                elif match_sequence:
                    data_fasta[head] += match_sequence.group(1)
    except IOError:
        sys.exit("Error cannot open {0}".format(aln_fasta_file))
    return data_fasta


def write_pir_file(aln_pir_file, data_fasta, pdb_codes, pdb_files,
                   add_heteroatom, heteroatom_models):
    """Write new pir alignment
     :Parameters:
      aln_pir_file: Path to the PIR file
      data_dict: Dictionary containing alignment data
      pdb_codes: List of PDB
      add_heteroatom: Number of heteroatom
      heteroatom_models: List of model for which heteroatom should be added
    """
    hetatm_flag = True
    hetatm = ""
    if add_heteroatom > 0:
        hetatm = "." * add_heteroatom
    try:
        with open(aln_pir_file, "wt") as aln_pir:
            for element in data_fasta:
#                 pdb_start_posit = 1
                aln_pir.write(">P1;{0}{1}".format(element, os.linesep))
                type_data = "sequence"
                # Identify sequences with a pdb structure associated
                if element in pdb_codes:
                    type_data = "structureX"
                # Check if we have to add heteroatom in the alignment for the
                # current PDB
                if (element not in heteroatom_models and
                        len(heteroatom_models) > 0):
                    hetatm_flag = False
                sequence = data_fasta[element].replace("-", "")
                sequence = sequence.replace(os.linesep, "")
                if hetatm_flag:
                    aln_pir.write("{0}:{1}: :A: :A{2}{3}".format(
                                type_data, element, ": " * 4, os.linesep))
                    aln_pir.write("{0}{1}*{2}".format(data_fasta[element],
                                                     hetatm, os.linesep))
                else:
                    aln_pir.write("{0}:{1}: :A: :A{2}{3}"
                                  .format(type_data, element,
                                          ": " * 4, os.linesep))
                    aln_pir.write("{0}*{1}".format(data_fasta[element],
                                                  os.linesep))
                hetatm_flag = True
    except IOError:
        sys.exit("Error cannot open {0}".format(aln_pir_file))
    except AssertionError:
        sys.exit("Error cannot find start residue position in"
                 " the pdb : {0}".format(pdb_files[pdb_codes.index(element)]))


def get_multifasta_data(multifasta_file):
    """
    """
    regex_head = re.compile(r"^>([\w-]+)")
    regex_protein = re.compile(r"^([A-Za-z]+)")
    multifasta_data = {}
    try:
        with open(multifasta_file, "rt") as multifasta:
            for line in multifasta:
                match_head = regex_head.match(line)
                match_protein = regex_protein.match(line)
                if match_head:
                    head = match_head.group(1)
                    multifasta_data[head] = ""
                elif match_protein:
                    multifasta_data[head] += match_protein.group(1)
        for seq_head in multifasta_data:
            assert(len(multifasta_data[seq_head]) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(multifasta_file))
    except AssertionError:
        sys.exit("The program has faild to extract the sequence of "
                 "\"{0}\".".format(seq_head))
    return multifasta_data


def get_pdb_sequence(pdb_file, seqdict):
    """
    """
    pdb_seq = ""
    res = 0
    chain = None
    try:
        with open(pdb_file, "rt") as pdb:
            for line in pdb:
                aa = line[17:20]
                newres = line[22:26]
                field = line[0:4]
                if newres != res and field == "ATOM":
                    res = newres
                    newchain = line[21:22]
                    if not chain:
                        chain = line[21:22]
                    if newchain != chain:
                        break
                    pdb_seq += seqdict[aa]
    except IOError:
        sys.exit("Error cannot open {0}".format(pdb_file))
    except TypeError:
        print("Unknown residue \"{0}\" :{1}{2}".format(aa, os.linesep, line))
    return pdb_seq


def adjust_multifasta_format(multifasta_file, multifasta_data, pdb_seq,
                             wrong_pdb, pdb_codes, pdb_files, seqdict, results):
    """Create a new multifasta
    """
    corrected_fasta_file = (results +
                            ".".join(os.path.basename(
                                        multifasta_file).split(".")[:-1])
                            + "_corrected.fasta")
    try:
        with open(corrected_fasta_file, "wt") as corrected_fasta:
            for head in multifasta_data:
                # PDB is missing in the multifasta
                if head in wrong_pdb:
                    corrected_fasta.write(
                        ">{0}{1}{2}{1}".format(
                            head, os.linesep,
                            "{0}".format(os.linesep).join(
                                textwrap.wrap(pdb_seq[head], 80))))
                    # Resolved problems
                    wrong_pdb.pop(wrong_pdb.index(head))
                # The line was correct
                else:
                    corrected_fasta.write(
                        ">{0}{1}{2}{1}".format(
                            head, os.linesep,
                            "{0}".format(os.linesep).join(
                                textwrap.wrap(multifasta_data[head], 80))))
            # Download missing sequences
            for pdb in wrong_pdb:
                pdb_index = pdb_codes.index(pdb)
                corrected_fasta.write(">{0}{1}{2}{1}".format(
                            pdb, os.linesep, "{0}".format(os.linesep).join(
                                textwrap.wrap(get_pdb_sequence(
                                    pdb_files[pdb_index], seqdict), 80))))
    except IOError:
        sys.exit("Error cannot open {0}".format(corrected_fasta_file))
    return corrected_fasta_file


def check_multifasta(multifasta_file, pdb_codes, pdb_files, seqdict,
                     disable_autocorrect, results):
    """Check multifasta
    """
    check_wrong = False
    wrong_pdb = pdb_codes[:]
    pdb_seq = {}
    # Check extension
    if not multifasta_file.endswith("fasta"):
        print("Warning :  the sequence are supposed to be "
              "in fasta format", file=sys.stderr)
    # Get multifasta data
    multifasta_data = get_multifasta_data(multifasta_file)
    # Check multifasta data and PDB
    for head in multifasta_data:
        try:
            pdb_seq[head] = get_pdb_sequence(pdb_files[pdb_codes.index(head)],
                                             seqdict)
            # IF PDB is OK
            if pdb_seq[head] == multifasta_data[head]:
                wrong_pdb.pop(wrong_pdb.index(head))
            else:
                print("Warning : Amino-acid sequence in the PDB and in the "
                      "multifasta file for {0} is different.".format(head),
                      file=sys.stderr)
                check_wrong = True
        except ValueError:
            pass
    if len(wrong_pdb) != 0 and not check_wrong:
        print("Warning : The multifasta file does not contain the sequence "
              "of every PDB files.", file=sys.stderr)
        check_wrong = True
    # Autocorrect
    if not disable_autocorrect and check_wrong:
        print("Try to correct multifasta file.", file=sys.stderr)
        multifasta_file = adjust_multifasta_format(multifasta_file,
                                                   multifasta_data, pdb_seq,
                                                   wrong_pdb, pdb_codes,
                                                   pdb_files, seqdict, results)
    return multifasta_file


def run_alignment(conf_data, multifasta_file, pdb_codes, pdb_files,
                  alignment_software, path_soft, add_heteroatom,
                  heteroatom_models, thread, results):
    """Compute alignment and adjust PIR information
     :Parameters:
      conf_data: Configuration dictionary
      multifasta_file: Path to the multifasta file
      pdb_codes: List of PDB
      pdb_files: List of PDB files
      alignment_software: Name of the alignment softwaire
      path_soft: Path to alignment software
      add_heteroatom: Number of heteroatom
      heteroatom_models: List of model for which heteroatom should be added
      thread: Number of process
      results: Output path
     Returns:
    """
    aln_pir_file = (results + alignment_software + "_" +
                    str(os.getpid()) + "_aln.pir")
    # compute
    if alignment_software in ("t_coffee", "clustalw2"):
        run_command(replace_motif(conf_data.hdict[alignment_software],
                                  path_soft, multifasta_file, pdb_files,
                                  aln_pir_file, thread, "", "", ""))
    else:
        aln_fasta_file = (results + alignment_software + "_"
                          + str(os.getpid()) + "_aln.fasta")
        run_command(replace_motif(conf_data.hdict[alignment_software],
                                  path_soft, multifasta_file, pdb_files,
                                  aln_fasta_file, thread, "", "", ""))
        data_fasta = get_fasta_data(aln_fasta_file)
        write_pir_file(aln_pir_file, data_fasta, pdb_codes, pdb_files,
                       add_heteroatom, heteroatom_models)
    if alignment_software in ("t_coffee", "clustalw2"):
        aln_pir_file = check_pir(aln_pir_file, pdb_codes, pdb_files,
                                 add_heteroatom, heteroatom_models)
    return aln_pir_file


def get_model(aln_file, pdb_codes):
    """Get the model
     :Parameters:
      aln_data: Alignment read file
      pdb_codes: List of PDB
     Returns:
    """
    regex = re.compile(r"^>\w+;([\w-]+)")
    model = None
    try:
        with open(aln_file) as aln:
            for line in aln:
                match = regex.match(line)
                if match:
                    # Get first element if not specified
                    if match.group(1) not in pdb_codes:
                        model = match.group(1)
                        break
        assert(model)
    except IOError:
        sys.exit("Error cannot open {0}".format(aln_file))
    except AssertionError:
        sys.exit("Model name not found")
    return model


def get_environment(pdb_files):
    """Set Modeller environment parameters
     :Parameters:
      pdb_files: List of PDB files
    """
    env = environ(restyp_lib_file="${LIB}/restyp.lib")
    # Set PDB directory
    env.io.atom_files_directory = [os.path.dirname(pdb) for pdb in pdb_files]
    # Read in HETATM records from template PDBs
    env.io.hetatm = True
    # Set topology parameters
    env.libs.topology.read(file="$(LIB)/top_heav.lib")
    env.libs.parameters.read(file="$(LIB)/par.lib")
    return env


def get_parallel(process):
    """Start modeling slave
     :Parameters:
      process: Number of process to launch
     Returns:
      List of worker
    """
    job_worker = job()
    for i in xrange(0, process):
        job_worker.append(local_slave())
    return job_worker


def write_extract_fasta(multifasta_data, model):
    """
    """
    if len(multifasta_data) == 1 and not model:
        model = multifasta_data.keys()[0]
        print("Warning model name is not indicated.{0}\"{1}\" "
              "sequence is considered for psipred"
              " prediction".format(os.linesep, model))
    elif(len(multifasta_data) > 1 and not model):
        sys.exit("The program has failed to extract the fasta "
                 "sequence for psipred.{0}Please, indicate the "
                 "name of the model".format(os.linesep))
    out_file = model + '_psipred.fasta'
    try:
        with open(out_file, "wt") as out:
            out.write(">{0}{1}".format(model, os.linesep))
            out.write("{0}".format(os.linesep).join(
                textwrap.wrap(multifasta_data[model], 80)))
    except IOError:
        sys.exit("Error cannot open {0}".format(out_file))
    except KeyError:
        sys.exit("The program has failed to extract the fasta sequence of "
                 ".{1}Check the multifasta file.".format(os.linesep))
    return out_file


def run_secondary_structure_pred(conf_data, multifasta_file, model,
                                 path_psipred):
    """
    """
    multifasta_data = get_multifasta_data(multifasta_file)
    fasta_file = write_extract_fasta(multifasta_data, model)
    run_command(replace_motif(conf_data.hdict['psipred'], path_psipred, "", "",
                              "", "", "", fasta_file, ""))
    return fasta_file.replace(".fasta", ".horiz")


def load_psipred(psipred_file):
    """Load psipred data
     :Parameters:
      psipred_file: Path to psipred file
     Returns:
      conf: List of confidence level associated to the conformation
      pred: List of residue conformations (helix, strand)
    """
    conf = []
    pred = []
    regex_conf = re.compile(r"^Conf:\s+([0-9]+)")
    regex_pred = re.compile(r"^Pred:\s+([ECH]+)")
    if (not psipred_file.endswith("psipass2")
        and not psipred_file.endswith("horiz")):
        print("Warning :  the psipred file is supposed to be "
              "in PSIPRED HFORMAT format", file=sys.stderr)
    try:
        with open(psipred_file) as psipred:
            for line in psipred:
                match_conf = regex_conf.match(line)
                match_pred = regex_pred.match(line)
                if match_conf:
                    conf += map(int, match_conf.group(1))
                elif match_pred:
                    pred += list(match_pred.group(1))
    except IOError:
        sys.exit("Error cannot open {0}".format(psipred_file))
    except ValueError:
        sys.exit("There is a text value in confidence "
                 "line : {0}{1}".format(os.linesep, line))
    return conf, pred


def cluster_psipred(conf, pred, limit_confidence):
    """Clusterize psipred data according to the threshold
     :Parameters:
      conf: List of confidence level associated to the conformation
      pred: List of residue conformations (helix, strand)
      limit_confidence: Confidence threshold
     Returns:
      psipred_clusters:
    """
    psipred_clusters = []
    struct = None
    for i in xrange(len(pred)):
        if(conf[i] >= limit_confidence and pred[i] in "HE" and not struct):
            struct = ["", 0, 0]
            struct[0] = pred[i]
            struct[1] = i + 1
        if(struct):
            if(struct[0] != pred[i] or conf[i] < limit_confidence):
                struct[2] = i
                psipred_clusters += [struct]
                struct = None
    return psipred_clusters


def compute_models(env, job_worker, alignment_file, pdb_codes, pdb_files,
                   model_name, model_quality, number_model, psipred_result):
    """Define modeling parameters and start modeling
    """
    # Load psipred
    if(psipred_result):
        atm = RestraintModel(env, alnfile=alignment_file,
                             knowns=pdb_codes, sequence=model_name,
                             assess_methods=[assess.GA341, assess.DOPE,
                                             assess.normalized_dope])
        atm.psipred_result = psipred_result
    else:
#         # classical models
        atm = automodel(env, alnfile=alignment_file,
                        knowns=pdb_codes, sequence=model_name,
                        assess_methods=[assess.GA341, assess.DOPE,
                                        assess.normalized_dope])
    # Indicate number of template and model
    atm.starting_model = 1
    atm.ending_model = int(number_model)
    # Start slave
    atm.use_parallel_job(job_worker)
    if model_quality == "very_fast":
        print("Start very fast modeling")
        atm.very_fast()
    elif model_quality == "fast":
        print("Start fast modeling")
    elif model_quality == "normal":
        print("Start normal speed modeling")
        atm.md_level = refine.slow
    else:
        print("Start max quality modeling")
        # Thorough MD optimization
        atm.md_level = refine.very_slow
        atm.repeat_optimization = 2
        atm.max_molpdf = 1e6
        # Very thorough VTFM optimization
        atm.library_schedule = autosched.slow
        atm.max_var_iterations = 300
    # Start modeling
    atm.make()
    return atm


def histplot(data, label, result_data):
    """Barplot general method
    """
    plt.xlabel(label[0])
    plt.ylabel(label[1])
    plt.title(label[2])
    plt.hist(data)
    plt.savefig(result_data)
    plt.clf()


def summary_data(data, results, sessionid):
    """Write summary results
     :Parameters:
      data:
      results: Output path
      sessionid: ID number
    """
    summary_file = (results + "modeller_summary_" + str(sessionid) + ".csv")
    try:
        with open(summary_file, "wt") as out:
            writer = csv.writer(out, delimiter='\t')
            writer.writerow(["Model", "molpdf", "DOPE score",
                             "Normalized DOPE score", "GA341 score"])
            for i in xrange(len(data)):
                writer.writerow([data[i]["name"], data[i]["molpdf"],
                                 data[i]["DOPE score"],
                                 data[i]["Normalized DOPE score"],
                                 data[i]["GA341 score"][0]])
    except IOError:
        sys.exit("Error cannot open {0}".format(summary_file))


def save_general_data(atm, results, sessionid):
    """Write general parameters of the produced models
     :Parameters:
       atm:
       results: Output path
       sessionid: ID number
    """
    # Get converged models
    ok_models = filter(lambda x: x['failure'] is None, atm.outputs)
    # Rank the models by molpdf score
    ok_models.sort(lambda a, b: cmp(a["DOPE score"], b["DOPE score"]))
    # Barplot Dope score of each ok models
    if MATPLOTLIB:
        histplot([i["DOPE score"] for i in ok_models],
                 ["Dope score", "Frequency", "Dope score histogram"],
                 results + "dope_per_model_{0}.svg".format(sessionid))
    # Write data summary
    summary_data(ok_models, results, sessionid)


def insert_gaps(vals, seq):
    """
    """
    # Insert gaps into the profile corresponding to those in seq:
    for n, res in enumerate(seq.residues):
        for gap in xrange(res.get_leading_gaps()):
            vals.insert(n, None)
    # Add a gap at position '0', so that we effectively count from 1:
    vals.insert(0, None)
    return vals


def get_profile(profile_file, seq):
    """Read `profile_file` into a Python array, and add gaps corresponding to
       the alignment sequence `seq`.
     :Parameters:
     Returns:
    """
    # Read all non-comment and non-blank lines from the file:
    vals = []
    try:
        with open(profile_file, "rt") as f:
            for line in f:
                if not line.startswith('#') and len(line) > 10:
                    spl = line.split()
                    vals.append(float(spl[-1]))
    except IOError:
        sys.exit("Error cannot open {0}".format(profile_file))

    return insert_gaps(vals, seq)


def compute_profile(data):
    """Compute the profile of each model
    """
    env = get_environment([data[0]])
    aln = alignment(env, file=data[2])
    aln_template = aln[os.path.basename(data[0]).split(".")[0]]
    # env and pdbfile
    pdb = complete_pdb(env, data[0])
    # all atom selection
    s = selection(pdb)
    # Energy profile
    profile = s.get_dope_profile()
    profile.write_to_file(".".join(os.path.basename(data[0]).split(".")[:-1])
                          + ".profile")
    # Get smooth curve
    profile_smooth = profile.get_smoothed(window=15)
    profile_smooth.write_to_file(
        ".".join(os.path.basename(data[0]).split(".")[:-1])
        + "_smooth.profile")
    return [insert_gaps([i.energy for i in profile], aln_template),
            insert_gaps([i.energy for i in profile_smooth], aln_template)]


def get_session_id(idfile):
    """Get the id of previous calculation
    """
    session_id = os.getpid()
    regex_aln = re.compile(r"[\w-]+_([0-9]+)_aln.*\.pir")
    regex_modeller = re.compile(r"[\w-]+_([0-9]+)\.csv")
    match_aln = regex_aln.match(os.path.basename(idfile))
    match_modeller = regex_modeller.match(os.path.basename(idfile))
    if match_aln:
        session_id = match_aln.group(1)
    elif match_modeller:
        session_id = match_modeller.group(1)
    return session_id


def load_profiles(pdb_files, pdb_codes, alignment_file, process, results):
    """Compute and load the profile of each model
    """
    list_run = [[pdb_files[i], (results + pdb_codes[i] + '.profile'),
                 alignment_file]
                for i in xrange(len(pdb_files))]
    pool = mp.Pool(processes=process)
    asyncResult = pool.map_async(compute_profile, list_run)
    profile = asyncResult.get()
    return [i[0] for i in profile], [i[1] for i in profile]


def plot_DOPE_profile(list_template, list_model, list_model_files, sessionid,
                      pdb_codes, results, note=None):
    """Plot the dope result for each residue of each model
    """
    # Get color map
    color_map = cm.get_cmap('gist_rainbow')
    colors = [color_map(1. * i / len(list_model))
              for i in xrange(len(list_model))]
    for t in xrange(len(list_template)):
        # Plot the template and model profiles in the same plot for comparison:
        fig = plt.figure(figsize=(20, 7))
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.set_xlabel('Alignment position')
        ax1.set_ylabel('DOPE score per-residue')
        if len(list_model) > 0:
            ax1.set_xlim(0, len(list_model[0]))
        # Plot templates
        temp_col = ax1.plot(list_template[t], color="green", linewidth=1,
                            label='Template')
        # Plot models
        for i in xrange(len(list_model)):
            model_col = ax1.plot(list_model[i], color=colors[i], linewidth=1,
                                 label='Model')
        if(len(list_model_files) <= 10):
            ax1.legend([pdb_codes[t]] +
                       [(os.path.basename(list_model_files[i]).split(".")[0]
                         + "_{0}".format(i + 1))
                        for i in xrange(len(list_model_files))],
                       loc="upper center", numpoints=1,
                       bbox_to_anchor=(0.5, 1.12),
                       ncol=3, fancybox=True, shadow=True)
        else:
            ax1.legend([temp_col, model_col[-1]], [pdb_codes[t]] +
                       [os.path.basename(list_model_files[0]).split(".")[0]
                        + "_*"],
                       loc="upper center", numpoints=1,
                       bbox_to_anchor=(0.5, 1.12),
                       ncol=3, fancybox=True, shadow=True)
        if note:
            plt.savefig(results + os.sep +
                        'dope_profile_{0}_{1}_{2}.svg'
                        .format(sessionid, pdb_codes[t], note))
        else:
            plt.savefig(results + os.sep +
                        'dope_profile_{0}_{1}.svg'.format(sessionid,
                                                          pdb_codes[t]))
        plt.clf()


def plot_DOPE_profile_all(list_template, list_model, list_model_files,
                          sessionid, pdb_codes, results, note=None):
    """Plot the dope result for each residue of each model
    """
    # Get color map
    color_map = cm.get_cmap('gist_rainbow')
    colors = [color_map(1. * i / (len(list_model) + len(list_template)))
              for i in xrange(len(list_model) + len(list_template))]
    # Plot the template and model profiles in the same plot for comparison:
    fig = plt.figure(figsize=(20, 7))
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.set_xlabel('Alignment position')
    ax1.set_ylabel('DOPE score per-residue')
    if len(list_model) > 0:
        ax1.set_xlim(0, len(list_model[0]))
    for t in xrange(len(list_template)):
        # Plot templates
        temp_col = ax1.plot(list_template[t], color=colors[t], linewidth=1,
                            label='Template')
        # Plot models
    for i in xrange(len(list_model)):
        model_col = ax1.plot(list_model[i], color=colors[t + i + 1],
                             linewidth=1, label='Model')
    if(len(list_model_files) <= 10):
        ax1.legend(pdb_codes +
                   [(os.path.basename(list_model_files[i]).split(".")[0]
                     + "_{0}".format(i + 1))
                    for i in xrange(len(list_model_files))],
                   loc="upper center", numpoints=1,
                   bbox_to_anchor=(0.5, 1.12),
                   ncol=3, fancybox=True, shadow=True)
    else:
        ax1.legend([temp_col, model_col[-1]], pdb_codes +
                   [os.path.basename(list_model_files[0]).split(".")[0]
                    + "_*"],
                   loc="upper center", numpoints=1,
                   bbox_to_anchor=(0.5, 1.12),
                   ncol=3, fancybox=True, shadow=True)
    if note:
        plt.savefig(results + os.sep +
                    'all_dope_profile_{0}_{1}_{2}.svg'
                    .format(sessionid, "_".join(pdb_codes), note))
    else:
        plt.savefig(results + os.sep +
                    'all_dope_profile_{0}_{1}.svg'.format(sessionid,
                                                      "_".join(pdb_codes)))
    plt.clf()


def plot_partial_DOPE_profile(list_template, list_model, list_model_files,
                              sessionid, pdb_codes, results):
    """Divide the protein in two side and plot dope
       :Parameters:
        list_template: List of template
        list_model: List of model
        list_model_files: List of model files
        sessionid: Id number
        pdb_codes: List of PDB codes
        results: Output path
    """
    # Plot templates
    for i in xrange(len(list_template)):
        middle = math.ceil(len(list_template[i]) / 2.0)
        part = [[0, int(middle)]]
        part += [[int(middle) + 1, len(list_template[i])]]
        for e in part:
            # Plot models
            for j in xrange(len(list_model)):
                dist = abs(len(list_model[j]) - len(list_template[i]))
                if(len(list_model[j]) < len(list_template[i])):
                    list_model[j] += list_template[i][-dist:]
                elif(len(list_model[j]) > len(list_template[i])):
                    list_template[i] += list_model[j][-dist:]
                assert(len(list_model[j]) == len(list_template[i]))
                # Plot the template and model profiles
                # in the same plot for comparison:
                fig = plt.figure(figsize=(20, 7))
                ax1 = fig.add_subplot(1, 1, 1)
                ax1.set_xlabel('Alignment position')
                ax1.set_ylabel('DOPE score per-residue')
                ax1.set_xlim(e[0], e[1])
                ax1.plot(list_template[i], color="green", linewidth=1,
                         label=pdb_codes[i])
                ax1.plot(list_model[j], color="black", linewidth=1,
                         label='Model')
                model_name = ".".join(os.path.basename(list_model_files[j])
                                      .split(".")[:-1])
                ax1.legend(
                    [pdb_codes[i]] +
                    [model_name.split(".")[0] + "_{0}".format(i + 1)],
                    loc="upper center", numpoints=1,
                    bbox_to_anchor=(0.5, 1.12),
                    ncol=2, fancybox=True, shadow=True)
                plt.savefig(results + os.sep +
                            'partial_dope_profile_{0}.svg'
                            .format("_".join([str(sessionid), pdb_codes[i],
                                              model_name, str(e[0] + 1),
                                              str(e[1] + 1)])))
                plt.clf()


def compute_delta_DOPE(template_profile, list_model_profile):
    """Compute delta DOPE
     :Parameters:
      template_profile:
      list_model_profile:
    """
    delta = []
    for model_profile in list_model_profile:
        dist = abs(len(model_profile) - len(template_profile))
        if(len(model_profile) < len(template_profile)):
            model_profile += template_profile[-dist:]
        elif(len(model_profile) > len(template_profile)):
            template_profile += model_profile[-dist:]
        assert(len(model_profile) == len(template_profile))
        profile = []
        for i in xrange(len(model_profile)):
            if(model_profile[i] and template_profile[i]):
                profile.append(model_profile[i] - template_profile[i])
            else:
                profile.append(None)
        delta += [profile]
    return delta


def plot_delta_DOPE_profile(list_delta_dope, list_model_files,
                            sessionid, results, template_pdb):
    """Plot delta dope profile
      :Parameters:
        list_delta_dope: List of delta dope
        list_model_files: List of model files
        sessionid: Id number
        results: Output path
        template_pdb: List of template PDB file
    """
    # Plot models
    for i in xrange(len(list_delta_dope)):
        # Plot the template and model profiles in the same plot for comparison:
        fig = plt.figure(figsize=(20, 7))
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.set_xlabel('Alignment position')
        ax1.set_ylabel('delta DOPE per-residue score')
        ax1.set_xlim(0, len(list_delta_dope[0]))
        model_name = ".".join(
                        os.path.basename(list_model_files[i]).split(".")[:-1])
        ax1.plot(list_delta_dope[i], color="black", linewidth=1,
                 label=template_pdb + model_name)
        plt.savefig(results + os.sep + "delta_dope_profile_{0}.svg"
                    .format("_".join([str(sessionid), template_pdb,
                                      model_name])))
        plt.clf()


def get_unique(seq):
    """Get unique elements with order preserving
      :Parameters:
       seq: List of elements
    """
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]


def get_unique_no_order(seq):
    """Get unique elements with no order preserving
     :Parameters:
           seq: List of elements
    """
    # Not order preserving
    return {}.fromkeys(seq).keys()


def get_sequence(pdb_file):
    """Extract PDB file sequence
     :Parameters:
      pdb_file: Path to the PDB file
    """
    seqpdb = []
    try:
        with open(pdb_file, "rt") as pdb:
            for line in pdb:
                if line[0:4] == "ATOM":
                    seqpdb.append(line[17:20] + "_" + str(int(line[23:26])))
            assert(len(seqpdb) > 0)
            seqpdb = get_unique(seqpdb)
    except IOError:
        sys.exit("Error cannot open {0}".format(pdb_file))
    except ValueError:
        sys.exit("Nothing read from {0}".format(pdb_file))
    return seqpdb


def save_delta_DOPE_profile(list_model_files, list_model_profile,
                            list_delta_dope, results, template_pdb):
    """Write delta DOPE profile
     :Parameters:
        list_model_files: List of model files
        list_model_profile: List of profile
        list_delta_dope: List of delta dope
        results: Output path
        template_pdb: List of template PDB file
    """
    try:
        for mod in list_model_files:
            seq_pdb = get_sequence(mod)
            for i in xrange(len(list_delta_dope)):
                num = 0
                model_name = ".".join(
                    os.path.basename(list_model_files[i]).split(".")[:-1])
                output_file = (results + os.sep + template_pdb + "_" +
                               model_name + ".delta_dope_profile")
                with open(output_file, "wt") as output:
                    for j in xrange(len(list_delta_dope[i])):
                        if(list_model_profile[i][j] and num < len(seq_pdb)):
                            if list_delta_dope[i][j]:
                                output.write("{0}\t{1}{2}".format(
                                    seq_pdb[num].split("_")[0],
                                    list_delta_dope[i][j], os.linesep))
                            else:
                                output.write("{0}\t0{1}".format(
                                    seq_pdb[num].split("_")[0], os.linesep))
                            num += 1
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


def load_summary(summary_file):
    """Load summary data
     :Parameters:
      summary_file: File describing Dope energy
    """
    summary_data = {}
    try:
        with open(summary_file, "rt") as summary:
            summary_reader = csv.reader(summary, delimiter='\t')
            # Pass head
            summary_reader.next()
            for line in summary_reader:
                if(len(line) == 5):
                    summary_data[line[0]] = [float(val) for val in line[1:]]
    except ValueError:
        sys.exit("Error cannot convert the element of {0}"
                 " into float".format(line))
    except IOError:
        sys.exit("Error cannot open {0}".format(summary_file))
    return summary_data


def save_picture(urlfile, output_file):
    """
    """
    try:
        request = requests.get(urlfile)
        with open(output_file, "wb") as output_image:
            for block in request.iter_content(1024):
                if not block:
                    break
                output_image.write(block)
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


def run_proq(website_path, pred, pdb):
    """Run ProQ
    """
    # Set PDB file
    try:
        with open(pdb, 'rt') as pdb_sent:
            files = {'pdbfile': pdb_sent}
            # Send PDB file
            if pred:
                req = requests.post(website_path, files=files,
                                    data={'ss': "".join(pred)})
            else:
                req = requests.post(website_path, files=files)
        assert(req.status_code == requests.codes.ok)
    except IOError:
        sys.exit("Error cannot open {0}".format(pdb_sent))
    except AssertionError:
        sys.exit("Something went wrong with {0}".format(0))
    try:
        if(req.text):
            maxsub = float(req.text.split("Predicted MaxSub              "
                                          ": <b>")[1].split("</b><BR>")[0])
            lgscore = float(req.text.split("Predicted LGscore             "
                                           ": <b>")[1].split("</b><BR>")[0])
        else:
            sys.exit("No data received from ProQ")
    except ValueError:
        sys.exit("Something went wrong with {0}".format(0))
    return [pdb, maxsub, lgscore]


def run_prosa(website_path, pdb, results):
    """Run prosa
    """
    zscore = None
    path_hrplot = None
    path_eplot = None
    # Set cache properties
    payload = {'max_file_size': '26214400', "maxlength": '26214400'}
    # Set PDB file
    try:
        # "action":"/prosa.phb"
        with open(pdb, 'rt') as pdb_sent:
            files = {'userfile': pdb_sent}
            # Send PDB file
            req = requests.post(website_path + "prosa.php", files=files,
                                data=payload)
        assert(req.status_code == requests.codes.ok)
    except IOError:
        sys.exit("Error cannot open {0}".format(pdb_sent))
    except AssertionError:
        sys.exit("Something went wrong with {0}".format(0))
    try:
        if (req.text):
            zscore = float(req.text.split("<span class='zscore'>")[1]
                           .split("</span>")[0])
            path_hrplot = req.text.split("<a href='upload/")[1].split(
                                                                "' alt='")[0]
            path_eplot = req.text.split("<img src='upload/")[2].split(
                                                                "' alt='")[0]
        else:
            sys.exit("No data received from verify3D")
    except ValueError:
        sys.exit("Something went wrong with {0}".format(0))
    # Download profile picture
    if(path_hrplot):
        save_picture(website_path + "upload/" + path_hrplot,
                     results + "prosa_hrplot_" +
                     ".".join(os.path.basename(pdb).split(".")[:-1]) + ".png")
    if(path_eplot):
        save_picture(website_path + "upload/" + path_eplot,
                     results + "prosa_eplot_" +
                     ".".join(os.path.basename(pdb).split(".")[:-1]) + ".png")
    return [pdb, zscore]


def load_verify3D(verify3D_result, type_data=0):
    """
    """
    request = r"^[a-z]+\s+([A-Z])\s+([0-9]+)\s+[0-9.-]+\s+([0-9.-]+)"
    if(type_data):
        request = r"^[a-z]+\s+([A-Z])\s+([0-9]+)\s+([0-9.-]+)"
    regex_verify = re.compile(request)
    data_verify3D = []
    try:
        for line in verify3D_result.split(os.linesep):
            match_verify = regex_verify.match(line)
            if match_verify:
                # Residue Position First_score Second_score
                data_verify3D += [[match_verify.group(1),
                                   int(match_verify.group(2)),
                                   float(match_verify.group(3))]]
        assert(len(data_verify3D) > 0)
    except ValueError:
        sys.exit("There is something wrong with :" + os.linesep
                 + "\"{0}\"".format(line))
    except AssertionError:
        sys.exit("Failed to load verify3D raw scores :{0}\"{1}\""
                 .format(os.linesep, line))
    return data_verify3D


def run_verify3D(website_path, pdb, results):
    """Run Verify3D
    """
    req = None
    req2 = None
    # Set cache properties
    payload = {'MAX_FILE_SIZE': '1073741824', 'pagestatus': 'sendit'}
    # Set PDB file
    try:
        with open(pdb, 'rt') as pdb_sent:
            files = {'pdbfile': pdb_sent}
            # Indicate original page
            headers = {'referer': website_path}
            # Send PDB file
            req = requests.post(website_path, data=payload, files=files,
                                headers=headers)
        assert(req.status_code == requests.codes.ok)
    except IOError:
        sys.exit("Error cannot open {0}".format(pdb_sent))
    except AssertionError:
        sys.exit("Something went wrong with {0}".format(0))
    if(req):
        # Get session id
        ver = req.text.split('action="/Verify_3D/temp/raw')[1].split(".dat")[0]
        # Download profile picture
        save_picture(website_path + "temp/vplot" + ver + ".png",
                     results + "verify3D_"
                     + ".".join(os.path.basename(pdb).split(".")[:-1])
                     + ".png")
        # Download raw score
        req2 = requests.get(website_path + "temp/raw" + ver + ".dat?",
                            headers=headers)
        # Download average score
        req3 = requests.get(website_path + "temp/avg" + ver + ".dat?",
                            headers=headers)
    else:
        sys.exit("No data received from verify3D")

    if(req2):
        # Load verify3D data
        data_verify3D = load_verify3D(req2.text)
        # Write verify3D data
        write_checking(data_verify3D, ["Residue", "Position", "Score"],
                       results + "verify3D_raw_"
                       + ".".join(os.path.basename(pdb).split(".")[:-1])
                       + ".txt")
    else:
        sys.exit("No data received from verify3D")
    if(req3):
        # Load verify3D data
        data_verify3D_smooth = load_verify3D(req3.text, 1)
        # Write verify3D data
        write_checking(data_verify3D_smooth,
                       ["Residue", "Position", "Average score"],
                       results + "verify3D_avg_"
                       + ".".join(os.path.basename(pdb).split(".")[:-1])
                       + ".txt")
    else:
        sys.exit("No data received from verify3D")
    return data_verify3D, data_verify3D_smooth


def run_checking(conf_data, summary_data, structure_check, path_check,
                 number_best, pred, REQUESTS, results):
    """
     :Parameters:
      - conf_data: Configuration dictionary
      - structure_check:
    """
    data_proq = []
    data_prosa = []
    data_verify3D = {}
    data_verify3D_smooth = {}
    num_struct = 0
    for pdb in sorted(summary_data.iteritems(), key=lambda x: x[1][1]):
        if num_struct >= number_best:
            break
        if('procheck' in structure_check):
            print("Run procheck for " + pdb[0])
            run_command(replace_motif(conf_data.hdict['procheck'],
                                      path_check, "", [pdb[0]], "", "", "",
                                      "", ""))
        if('proq' in structure_check and REQUESTS):
            print("Run ProQ for " + pdb[0])
            status = True
            try:
                result_proq = run_proq(conf_data.hdict['proq'], pred, pdb[0])
            except requests.exceptions.ConnectionError:
                print("Error cannot connect to ProQ for {0}".format(pdb[0]),
                      file=sys.stderr)
                status = False
            if status:
                data_proq += [result_proq]
        if('prosa' in structure_check and REQUESTS):
            print("Run prosa for " + pdb[0])
            status = True
            try:
                result_prosa = run_prosa(conf_data.hdict['prosa'], pdb[0],
                                         results)
            except requests.exceptions.ConnectionError:
                print("Error cannot connect to prosa for {0}".format(pdb[0]),
                      file=sys.stderr)
                status = False
            if status:
                data_prosa += [result_prosa]
        if('verify3D' in structure_check and REQUESTS):
            print("Run verify3D for " + pdb[0])
            try:
                (data_verify3D[pdb[0]],
                 data_verify3D_smooth[pdb[0]]) = run_verify3D(
                                            conf_data.hdict['verify3D'],
                                            pdb[0], results)
            except requests.exceptions.ConnectionError:
                print("Error cannot connect to Verify3D", file=sys.stderr)
                data_verify3D = {}
        num_struct += 1
    return data_proq, data_prosa, data_verify3D, data_verify3D_smooth


def write_checking(data, title, filename):
    """Write checking value
    """
    try:
        with open(filename, "wt") as output:
            outwriter = csv.writer(output, delimiter="\t")
            outwriter.writerow(title)
            outwriter.writerows(data)
    except IOError:
        sys.exit("Error cannot open {0}".format(filename))


def plot_verify3D_profile(summary_data, data_verify3D, number_best, results,
                          sessionid, ylabel, type_data):
    """Plot verify3D profile for several PDB
    """
    pdb_files = []
    num_struct = 0
    # Get color map
    color_map = cm.get_cmap('gist_rainbow')
    colors = [color_map(1. * i / number_best)
              for i in xrange(number_best)]
    # Plot the template and model profiles in the same plot for comparison:
    fig = plt.figure(figsize=(30, 7))
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.set_xlabel('Residue number')
    ax1.set_ylabel(ylabel)
    for pdb in sorted(summary_data.iteritems(), key=lambda x: x[1][1]):
        if num_struct >= number_best:
            break
        # Plot templates
        ax1.plot(
            [element[1] for element in data_verify3D[pdb[0]]],
            [element[2] for element in data_verify3D[pdb[0]]],
            color=colors[num_struct], linewidth=1, label=pdb[0])
        num_struct += 1
        pdb_files += [pdb[0]]
    ax1.legend(
         pdb_files,
         loc="upper center", numpoints=1,
         bbox_to_anchor=(0.5, 1.12),
         ncol=3, fancybox=True, shadow=True)
    plt.savefig(results + os.sep + "verify3D_profile_{0}_{1}.svg"
                .format(type_data, sessionid))
    plt.clf()


def load_interest_atom(interest_atom_file):
    """
    """
    atom_dict = {}
    try:
        with open(interest_atom_file, "rt") as interest_atom:
            reader = csv.reader(interest_atom, delimiter='\t')
            for line in reader:
                if len(line) == 3:
                    atom_dict[(int(line[0]), line[1].upper(),
                               line[2].upper())] = ""
        assert(len(atom_dict) != 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(interest_atom_file))
    except ValueError:
        sys.exit("Error an integer is expected instead "
                 "of : \"{0}\"".format(line[0]))
    except AssertionError:
        sys.exit("The load of the interest atom has failed,"
                 " atom_dict is empty")
    return atom_dict


def load_PDB_interest_atom(pdb_file, interest_atom_dict):
    """
    """
    atom_list = []
    try:
        with open(pdb_file, "rt") as pdb:
            for line in pdb:
                if line[0:6] == "ATOM  " or line[0:6] == "HETATM":
                    residue_number = int(line[22:26])
                    residue = line[17:20].replace(" ", "").upper()
                    atom = line[11:16].replace(" ", "").upper()
                    if(interest_atom_dict.has_key((residue_number,
                                                  residue, atom))):
                        atom_list += [["_".join([str(residue_number),
                                                residue, atom]),
                                      float(line[30:38]),
                                      float(line[38:46]),
                                      float(line[46:54])]]
        assert(len(atom_list) != 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(pdb_file))
    except ValueError:
        sys.exit("Error an numeric value is expected instead "
                 "of :{0}\"{1}\"".format(os.linesep, line))
    except AssertionError:
        sys.exit("The load of the interest atom has failed,"
                 " atom_list is empty")
    return atom_list


def load_default_distances(default_distances_file):
    """Load distance file
    """
    default_distances_dict = {}
    try:
        with open(default_distances_file, "rt") as default_distances:
            distances_reader = csv.reader(default_distances, delimiter='\t')
            for line in distances_reader:
                if len(line) == 3:
                    default_distances_dict[(line[0].upper(),
                                            line[1].upper())] = float(line[2])
        assert(len(default_distances_dict) != 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(default_distances_file))
    except ValueError:
        sys.exit("Error an numeric value is expected instead "
                 "of :{0}\"{1}\"".format(os.linesep, line))
    except AssertionError:
        sys.exit("The load of default distances has failed, dict is empty")
    return default_distances_dict


def get_distances(pdb_atom_list):
    """Compute distance between a list of atoms
    """
    list_dist = []
    for i in xrange(0, len(pdb_atom_list)):
        for j in xrange(i + 1, len(pdb_atom_list)):
            list_dist += [[pdb_atom_list[i][0], pdb_atom_list[j][0],
                           math.sqrt(pow(pdb_atom_list[i][1]
                                         - pdb_atom_list[j][1], 2.0) +
                                     pow(pdb_atom_list[i][2]
                                         - pdb_atom_list[j][2], 2.0) +
                                     pow(pdb_atom_list[i][3]
                                         - pdb_atom_list[j][3], 2.0))]]
    return list_dist


def write_distances(results, pdb_file, pdb_distances_list):
    """Write ATOM distances
    """
    distances_out_file = (results +
                          ".".join(os.path.basename(pdb_file).split(".")[:-1])
                          + "_distances.csv")
    try:
        with open(distances_out_file, "wt") as distances_out:
            distance_writer = csv.writer(distances_out, delimiter='\t')
            distance_writer.writerow(["ATOM_1", "ATOM_2",
                                      "Distance (angstrom)"])
            distance_writer.writerows(pdb_distances_list)
    except IOError:
        sys.exit("Error cannot open {0}".format(distances_out_file))


def compute_distance_variation(pdb_distances_list, default_distances_dict):
    """Compute distance 
    """
    distance_variation = []
    try:
        for element in pdb_distances_list:
            if(default_distances_dict.has_key((element[0], element[1]))):
                distance_variation += [abs(element[2] -
                                           default_distances_dict[
                                            (element[0], element[1])])]
        assert(len(distance_variation) != 0)
    except AssertionError:
        sys.exit("The calculation of distance variation has failed."
                 + os.linesep + "No element has been found in the atom list")
    return distance_variation


def run_check_distances(interest_atom_dict, default_distances, results,
                        sessionid):
    """Run distance checking
    """
    distance_variation = []
    # Load default distances
    if default_distances:
        default_distances_dict = load_default_distances(default_distances)
    # Check for each PDB the distances
    for pdb in sorted(summary_data.iteritems(), key=lambda x: x[1][1]):
        # Load atom position
        pdb_atom_list = load_PDB_interest_atom(pdb[0], interest_atom_dict)
        # Compute position
        pdb_distances_list = get_distances(pdb_atom_list)
        # Write distances
        write_distances(results, pdb[0], pdb_distances_list)
        if default_distances_dict:
            distance_variation += [pdb[0],
                                   compute_distance_variation(
                                    pdb_distances_list,
                                    default_distances_dict)]
    # Save distance variation
    write_checking(distance_variation, ["PDB", "Distance variation"],
                   results + os.sep + "result_prosa_{0}.txt"
                           .format(sessionid))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    pdb_codes = None
    pdb_files = None
    conf = None
    pred = None
    psipred_result = None
    # List of Amino-acid and their three-letter and one-letter code
    seqdict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
               "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G",
               "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
               "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
               "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
    # Load parameters
    args, parser = get_arguments()
    if (("profile" in args.list_operations
        or "model" in args.list_operations) and not args.pdb
        and not args.pdb_identification):
        print("You must provide at least one structure template "
              "for modeling and/or profiling", file=sys.stderr)
        sys.exit(parser.print_help())
    # Configure option
    conf_data = ModelingConfig(args.config, args.results)
    # Identify templates
    if(args.pdb_identification):
        args.pdb = identify_template(conf_data, args.multifasta_file,
                                     args.thread, args.pdb_identification,
                                     args.pdb_identification_database,
                                     args.pdb_identification_path,
                                     args.pdb_identification_strategy,
                                     args.results)
    # Prepare modeling
    if ("profile" in args.list_operations or "model" in args.list_operations
        and args.pdb):
        pdb_codes, pdb_files = get_pdb(conf_data, args.pdb, args.results)
    # Compute alignment
    if args.multifasta_file and not args.alignment_file:
        if pdb_codes and pdb_files:
            args.multifasta_file = check_multifasta(args.multifasta_file,
                                                    pdb_codes, pdb_files,
                                                    seqdict,
                                                    args.disable_autocorrect,
                                                    args.results)
        print("Run alignment")
        if ("profile" in args.list_operations
            or "model" in args.list_operations):
            args.alignment_file = run_alignment(conf_data,
                                                args.multifasta_file,
                                                pdb_codes, pdb_files,
                                                args.alignment_software,
                                                args.path_alignment,
                                                args.add_heteroatom,
                                                args.heteroatom_models,
                                                args.thread, args.results)
    elif args.multifasta_file and args.alignment_file:
        sys.exit("You have provided an alignment and the multifasta "
                 "which indicates to compute the alignment. "
                 "Please choose only one.")
    elif (not args.multifasta_file and not args.alignment_file
          and ("profile" in args.list_operations
               or "model" in args.list_operations)):
        print("You must provide an alignment or a multifasta sequence "
              "to continue...", file=sys.stderr)
        sys.exit(parser.print_help())
    if (("profile" in args.list_operations or "model" in args.list_operations)
        and not args.model_name):
        args.model_name = get_model(args.alignment_file, pdb_codes)
    # Run psipred
    if args.path_psipred and args.multifasta_file:
        print("Run psipred")
        args.psipred = run_secondary_structure_pred(conf_data,
                                                    args.multifasta_file,
                                                    args.model_name,
                                                    args.path_psipred)
    elif args.path_psipred and not args.multifasta_file:
        sys.exit("Please indicate the path to multifasta/fasta file containing"
                 " the sequence of the protein of interest")
    # Load psipred file
    if args.psipred:
        conf, pred = load_psipred(args.psipred)
        psipred_result = cluster_psipred(conf, pred, args.limit_confidence)
    # Get session id
    if args.alignment_file:
        sessionid = get_session_id(args.alignment_file)
    elif args.modeller_summary:
        sessionid = get_session_id(args.modeller_summary)
    else:
        sys.exit("One alignment file or a modeller summary file is required.")
    # Summary file reference
    summary_file = (args.results + "modeller_summary_"
                    + str(sessionid) + ".csv")
    if (MODELLER and ("model" in args.list_operations
                      or "profile" in args.list_operations)):
        # request verbose output
        log.level(output=1, notes=1, warnings=1, errors=1, memory=0)
        log.verbose()
        # Load environment
        env = get_environment(pdb_files)
    if MODELLER and "model" in args.list_operations:
        # Use several CPUs in a parallel job on this machine
        job_worker = get_parallel(int(args.thread))
        # Start Modelling
        atm = compute_models(env, job_worker, args.alignment_file, pdb_codes,
                             pdb_files, args.model_name, args.model_quality,
                             args.number_model, psipred_result)
        save_general_data(atm, args.results, sessionid)
    if MODELLER and MATPLOTLIB and "profile" in args.list_operations:
        # List of models
        list_model = [args.model_name + ".B" + str(99990000 + i)
                      for i in xrange(1, args.number_model + 1)
                      if os.path.isfile(args.model_name + ".B"
                                        + str(99990000 + i) + ".pdb")]
        list_model_files = [args.results + i + ".pdb" for i in list_model]
        # Load template profile
        (list_template_profile,
         list_template_profile_smooth) = load_profiles(pdb_files,
                                                       pdb_codes,
                                                       args.alignment_file,
                                                       args.thread,
                                                       args.results)
        # Load model profile
        (list_model_profile,
         list_model_profile_smooth) = load_profiles(list_model_files,
                                                    list_model,
                                                    args.alignment_file,
                                                    args.thread,
                                                    args.results)
        # Plot DOPE profile
        plot_DOPE_profile(list_template_profile, list_model_profile,
                          list_model_files, sessionid, pdb_codes, args.results)
        # Plot DOPE smooth profile
        plot_DOPE_profile(list_template_profile_smooth,
                          list_model_profile_smooth, list_model_files,
                          sessionid, pdb_codes, args.results, "smooth")
        # Plot DOPE profile
        plot_DOPE_profile_all(list_template_profile, list_model_profile,
                              list_model_files, sessionid, pdb_codes,
                              args.results)
        # Plot DOPE smooth profile
        plot_DOPE_profile_all(list_template_profile_smooth,
                              list_model_profile_smooth, list_model_files,
                              sessionid, pdb_codes, args.results, "smooth")
        # Plot partial DOPE
        plot_partial_DOPE_profile(list_template_profile, list_model_profile,
                                  list_model_files, sessionid, pdb_codes,
                                  args.results)
        # Delta DOPE
        for i in xrange(len(list_template_profile)):
            list_delta_dope = compute_delta_DOPE(list_template_profile[i],
                                            list_model_profile)
            plot_delta_DOPE_profile(list_delta_dope, list_model_files,
                                    sessionid, args.results, pdb_codes[i])
            save_delta_DOPE_profile(list_model_files, list_model_profile,
                                    list_delta_dope, args.results,
                                    pdb_codes[i])
        # Histogram of DOPE
        if (os.path.isfile(summary_file)):
            summary_data = load_summary(summary_file)
            histplot([summary_data[model][1]
                      for model in summary_data],
                     ["Dope score", "Frequency", "Dope score histogram"],
                     args.results + "dope_per_model_{0}.svg".format(sessionid))
        else:
            sys.exit("{0} does not exist".format(summary_file))
    if args.structure_check and os.path.isfile(summary_file):
        if not args.psipred and "proq" in args.structure_check:
            print("It is recommanded to provide the secondary structure "
                  "prediction using psipred to use ProQ.")
        # Histogram of DOPE
        summary_data = load_summary(summary_file)
        # Check number of model
        if args.number_best > args.number_model:
            args.number_best = args.number_model
        # Start structure checking
        (data_proq, data_prosa,
         data_verify3D, data_verify3D_smooth) = run_checking(
                                        conf_data, summary_data,
                                        args.structure_check, args.path_check,
                                        args.number_best, pred,
                                        REQUESTS, args.results)
        if(args.list_atom):
            interest_atom_dict = load_interest_atom(args.list_atom)
            run_check_distances(interest_atom_dict, args.default_distances,
                                args.results, sessionid)
        if data_proq:
            write_checking(data_proq, ["PDB", "maxsub", "lgscore"],
                           args.results + os.sep + "result_proq_{0}.txt"
                           .format(sessionid))
        if data_prosa:
            write_checking(data_prosa, ["PDB", "zscore"],
                           args.results + os.sep + "result_prosa_{0}.txt"
                           .format(sessionid))
        if data_verify3D:
            plot_verify3D_profile(summary_data, data_verify3D,
                                  args.number_best, args.results,
                                  sessionid, "Raw score per-residue",
                                  "raw")
        if data_verify3D_smooth:
            plot_verify3D_profile(summary_data, data_verify3D_smooth,
                                  args.number_best, args.results,
                                  sessionid, "Average score per-residue",
                                  "avg")
    elif args.structure_check and not os.path.isfile(summary_file):
        sys.exit("Summary file is required for structure checking")


if __name__ == '__main__':
    main()
