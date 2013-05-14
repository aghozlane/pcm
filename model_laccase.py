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
import multiprocessing as mp
try:
    from modeller import *
    from modeller.parallel import *
    from modeller.automodel import *
    from modeller.scripts import complete_pdb
    MODELLER = True
except ImportError:
    MODELLER = False
    print("Could not import modeller\nNo modeling will be available",
          file=sys.stderr)
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    MATPLOTLIB = True
except ImportError:
    MATPLOTLIB = False
    print("Could not import matplotlib.\n"
      "The statistic measures won't be drawn.", file=sys.stderr)

class ModelingConfig:
    """Configure alignment.
    """

    def __init__(self, config_file, results):
        """Instantiate ModelingConfig object
           Arguments:
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


    def writeconfig(self):
        """Write modeling config
        """
        self.config.add_section('PDB_config')
        self.config.set("PDB_config", "download_url",
                        "http://www.rcsb.org/pdb/files/")
        self.config.add_section('Alignment_config')
        self.config.set('Alignment_config', 'clustalo',
                        "%path_soft{0}clustalo -i %multifasta -o %output "
                        "--threads=%proc --auto -t Protein "
                        "--outfmt=fa".format(os.sep))
        self.config.set('Alignment_config', 'clustalw2',
                        "%path_soft{0}clustalw2 -INFILE=%multifasta "
                        "-OUTPUT=PIR -OUTFILE=%output".format(os.sep))
        self.config.set('Alignment_config', 'mafft',
                        "%path_soft{0}mafft --auto --thread %proc "
                        "%multifasta > %output".format(os.sep))
        self.config.set('Alignment_config', 'muscle',
                        "%path_soft{0}muscle -in %multifasta "
                        "-out %output".format(os.sep))
        self.config.set('Alignment_config', 't_coffee',
                        "%path_soft{0}t_coffee %multifasta "
                        "-outfile %output -output pir_aln  "
                        "-n_core %proc".format(os.sep))
        # -mode 3dcoffee -template_file %pdb
        # Write data
        try:
            # Writing our configuration file to 'example.cfg'
            with open(self.phyconfig_file, 'wt') as configfile:
                self.config.write(configfile)
        except IOError:
            sys.exit("Error : cannot open "
                     "file {0}".format(self.phyconfig_file))


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


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    local_path = ".{0}".format(os.sep)
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "model_laccase.py ")
    parser.set_defaults(pdb_dir=local_path, results=local_path,
                        number_model=8, thread=detect_cpus(),
                        path_t_coffee=local_path, model_quality="normal",
                        alignment_software="t_coffee",
                        path_alignment=local_path,
                        list_operations=["modeling", "profile"],
                        add_heteroatom=0, heteroatom_models=[])
    parser.add_argument("-l", "--list_operations", type=str, nargs='+',
                        choices=["modeling", "profile", "verification"],
                        help='Select operations.')
    parser.add_argument('-f', '--multifasta_file', type=isfile,
                        help='Fasta sequence.')
    parser.add_argument('-i', '--alignment_file', type=isfile,
                        help='Alignment file in PIR format.')
    parser.add_argument('-a', '--alignment_software', type=str,
                        choices=["clustalw2", "clustalo", "mafft", "muscle",
                                 "t_coffee"],
                        help="Indicate the software that should be used to "
                        "align sequences.")
    parser.add_argument('-ht', '--add_heteroatom', type=int,
                        help="Indicate the number of hetero-atom residue(s) that "
                        "should be added to the alignment software.")
    parser.add_argument('-hm', '--heteroatom_models',
                        type=str, nargs='+',
                        help="Indicate the models for which hetero-atom "
                        "residue(s) should be added to the alignment "
                        "software.")
    parser.add_argument('-p', '--pdb', type=str, required=True, nargs='+',
                        help="List of pdb files or codes.")
    parser.add_argument('-e', '--model_name', type=str,
                        help='Code of the sequence to modelize.')
    parser.add_argument('-g', '--path_alignment', type=isdir,
                        help='Path to t_coffee software.')
    parser.add_argument('-n', '--number_model', type=int,
                        help='Number of model to produce.')
    parser.add_argument('-q', '--model_quality', type=str,
                        choices=["very_fast", "fast", "normal", "max"],
                        help='Adjust the quality of the modeling.')
    parser.add_argument('-s', '--structure_check', type=str,
                        nargs='+', choices=["proq", "procheck", "verify3d"],
                        help='Select phylogeny software.')
    parser.add_argument('-r', '--results', type=isdir,
                        help='Path to result directory.')
    parser.add_argument('-k', '--path_check', type=isdir,
                        nargs='+', help='Path to alignment software.')
    parser.add_argument('-t', '--thread', type=int, help='Number of thread.')
    parser.add_argument('-c', '--config', type=isfile,
                        help='Path to configuration file.')
    return parser.parse_args(), parser


def detect_cpus():
    """
    Detects the number of CPUs on a system. Cribbed from pp.
    """
    # Linux, Unix and MacOS: # Linux, Unix and MacOS:
    if hasattr(os, "sysconf"):
        if "SC_NPROCESSORS_ONLN" in os.sysconf_names:
            # Linux & Unix: # Linux and Unix:
            ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
            if isinstance (ncpus, int) and ncpus > 0 :
                return ncpus
            else :  # OSX:
                return int(os.popen2("sysctl -n hw.ncpu")[ 1 ].read())
        # Windows:
        if "NUMBER_OF_PROCESSORS" in os.environ:
            ncpus = int(os.environ["NUMBER_OF_PROCESSORS"]);
            if ncpus > 0 :
                return ncpus
    return 1  # Default 1


def download_pdb(conf_data, pdb, results):
    """
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
        print("Something went wrong with %s" % outfilename,
              file=sys.stderr)
    return outfilename


def get_pdb(conf_data, pdb_list, results):
    """
    """
    pdb_codes = []
    pdb_files = []
    for pdb in pdb_list:
        if pdb.endswith('.pdb') and os.path.isfile(pdb):
            pdb_files += [pdb]
            pdb_codes += [os.path.basename(pdb).split(".")[0]]
        else:
            pdb_codes += [os.path.basename(pdb).split(".")[0]]
            pdb_files += [download_pdb(conf_data,
                                       os.path.basename(pdb).split(".")[0],
                                       results)]
    return pdb_codes, pdb_files


def run_command(cmd):
    """Run command
      Arguments:
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
                  output, thread):
    """Set the software command
    """
    print(build_command, file=sys.stderr)
    build_command = build_command.replace('%path_soft', path_soft)
    if thread:
        build_command = build_command.replace('%proc', str(thread))
    else:
        build_command = build_command.replace('%proc', str(detect_cpus()))
    build_command = build_command.replace('%multifasta', multifasta_file)
    build_command = build_command.replace('%pdb', " ".join(pdb_files))
    build_command = build_command.replace('%output', output)
    print(build_command, file=sys.stderr)
    return build_command


def check_format(aln_data, pdb_codes):
    """
    """
    ide = ""
    status = True
    count = 0
    ens = ["P1", "F1", "DL", "DC", "RL", "RC", "XX"]
    data = {}
    pdb_elements = pdb_codes[:]
    regex_suite = (re.compile(r"^>" + "|".join([i + ";([\w-]+)"for i in ens])),
                   re.compile(r"^([sequence|structureX]:[\w-]+" + ":\S+"*8
                              + ")"),
                   re.compile(r"([\w*-]+\n)"))
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
                      "format required for modeller :\n{0}".format(line),
                      file=sys.stderr)
                count += 1
                status = False
        # Match Sequence
        if count == 2:
            match = regex_suite[count].match(line)
            if not match:
                print("Warning: next line does not correspond to a correct "
                      "protein sequence :\n{0}".format(line), file=sys.stderr)
                status = False
            else:
                data[ide][1] += match.group(1)
    assert(len(pdb_elements) == 0)
    return status, data


def adjust_format(aln_pir_file, data_dict, pdb_codes,
                  add_heteroatom, heteroatom_models):
    """Correct pir format and add heteroatom in the alignment
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
                if(element not in heteroatom_models
                   and len(heteroatom_models) > 0):
                    het_atm_flag = False
                aln_pir.write(">P1;{0}\n".format(element))
                if data_dict[element][0] and not add_heteroatom:
                    aln_pir.write(data_dict[element][0])
                else:
                    # Write midline characteristic
                    type_data = "sequence"
                    if element in pdb_codes:
                        type_data = "structureX"
                        structure_present = True
                    sequence = data_dict[element][1].replace("-", "")
                    sequence = sequence.replace("\n", "")
                    sequence = sequence.replace("*", "")
                    if het_atm_flag:
                        aln_pir.write("{0}:{1}:1 :A:{2}:A"
                                      .format(type_data, element,
                                              len(sequence) + add_heteroatom)
                                      + ": "*4 + "\n")
                    else:
                        aln_pir.write("{0}:{1}:1 :A:{2}:A"
                                      .format(type_data, element,
                                              len(sequence))
                                      + ": "*4 + "\n")
                end_aln_posit = data_dict[element][1].index("*")
                if het_atm_flag:
                    aln_pir.write("{0}{1}{2}\n".format(
                        data_dict[element][1][0:end_aln_posit],
                        het_atm, data_dict[element][1][end_aln_posit:]))
                else:
                    aln_pir.write("{0}{1}\n".format(
                        data_dict[element][1][0:end_aln_posit],
                        data_dict[element][1][end_aln_posit:]))
                het_atm_flag = True
    except IOError:
        sys.exit("Error cannot open {0}".format())
    assert(structure_present)
    return output_file


def check_pir(aln_pir_file, pdb_codes, add_heteroatom, heteroatom_models):
    """Run checking of pir alignment file
    """
    try:
        if not aln_pir_file.endswith('.pir'):
            print("Warning : The alignment file is expected "
                  "to be in pir format", file=sys.stderr)
        with open(aln_pir_file, "rt") as aln_pir:
            aln_data = aln_pir.readlines()
        status, data_dict = check_format(aln_data, pdb_codes)
        if not status or add_heteroatom > 0:
            print("Try to correct pir alignment format.", file=sys.stderr)
            aln_pir_file = adjust_format(aln_pir_file, data_dict, pdb_codes,
                                         add_heteroatom, heteroatom_models)
    except AssertionError:
        sys.exit("All PDB structures must be referenced in the alignment.")
    except IOError:
        sys.exit("Error cannot open {0}".format(aln_pir_file))
    return aln_pir_file


def get_fasta_data(aln_fasta_file):
    """Extract fasta data
    """
    head = ""
    data_fasta = {}
    regex_head = re.compile(r"^>([\w-]+)")
    regex_sequence = re.compile(r"([\w*-]+\n)")
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


def write_pir_file(aln_pir_file, data_fasta, pdb_codes,
                   add_heteroatom, heteroatom_models):
    """Write new pir alignment
    """
    hetatm_flag = True
    hetatm = ""
    if add_heteroatom > 0:
        hetatm = "." * add_heteroatom
    try:
        with open(aln_pir_file, "wt") as aln_pir:
            for element in data_fasta:
                aln_pir.write(">P1;{0}\n".format(element))
                type_data = "sequence"
                # Identify sequences with a pdb structure associated
                if element in pdb_codes:
                    type_data = "structureX"
                #
                if (element not in heteroatom_models
                    and len(heteroatom_models) > 0):
                    hetatm_flag = False
                sequence = data_fasta[element].replace("-", "")
                sequence = sequence.replace("\n", "")
                if hetatm_flag:
                    aln_pir.write("{0}:{1}:1 :A:{2}:A".format(
                                type_data, element,
                                len(sequence) + add_heteroatom)
                                + ": "*4 + "\n")
                    aln_pir.write("{0}{1}*\n".format(data_fasta[element],
                                                     hetatm))
                else:
                    aln_pir.write("{0}:{1}:1 :A:{2}:A"
                                  .format(type_data, element, len(sequence))
                                  + ": "*4 + "\n")
                    aln_pir.write("{0}*\n".format(data_fasta[element]))
                hetatm_flag = True
    except IOError:
        sys.exit("Error cannot open {0}".format(aln_pir_file))


def run_alignment(conf_data, multifasta_file, pdb_codes, pdb_files,
                  alignment_software, path_soft, add_heteroatom,
                  heteroatom_models, thread, results):
    """Compute alignment and adjust pir information
    """
    aln_pir_file = (results + alignment_software + "_" +
                    str(os.getpid()) + "_aln.pir")
    if not multifasta_file.endswith("fasta"):
        print("Warning :  the sequence are supposed to be "
              "in fasta format", file=sys.stderr)
    # compute
    if alignment_software in ("t_coffee", "clustalw2"):
        run_command(replace_motif(conf_data.hdict[alignment_software],
                                  path_soft, multifasta_file, pdb_files,
                                  aln_pir_file, thread))
    else:
        aln_fasta_file = (results + alignment_software + "_"
                          + str(os.getpid()) + "_aln.fasta")
        run_command(replace_motif(conf_data.hdict[alignment_software],
                                  path_soft, multifasta_file, pdb_files,
                                  aln_fasta_file, thread))
        data_fasta = get_fasta_data(aln_fasta_file)
        write_pir_file(aln_pir_file, data_fasta, pdb_codes,
                       add_heteroatom, heteroatom_models)
    if alignment_software in ("t_coffee", "clustalw2"):
        aln_pir_file = check_pir(aln_pir_file, pdb_codes,
                                 add_heteroatom, heteroatom_models)
    return aln_pir_file


def get_model(aln_file, pdb_codes):
    """
    """
    regex = re.compile("^>\w+;([\w-]+)")
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


def get_parallel(thread):
    """Start modeling slave
    """
    job_worker = job()
    for i in xrange(0, thread):
        job_worker.append(local_slave())
    return job_worker


def compute_models(env, job_worker, alignment_file, pdb_codes, pdb_files,
                   model_name, model_quality, number_model):
    """Define modeling parameters and start modeling 
    """
    # Compute models
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
    """
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
       the alignment sequence `seq`."""
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
    # profile result
#     s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=data[1],
#                   normalize_profile=False, smoothing_window=15)
#     return get_profile(data[1], aln[os.path.basename(data[0]).split(".")[0]])
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


def get_session_id(alignment_file):
    """Get the id of previous calculation
    """
    session_id = os.getpid()
    regex = re.compile("[\w-]+_([0-9]+)_aln.*\.pir")
    match = regex.match(os.path.basename(alignment_file))
    if match:
        session_id = match.group(1)
    return session_id


def load_profiles(pdb_files, pdb_codes, alignment_file, thread, results):
    """Compute and load the profile of each model
    """
    list_run = [[pdb_files[i], (results + pdb_codes[i] + '.profile'),
                 alignment_file]
                for i in xrange(len(pdb_files))]
    pool = mp.Pool(processes=thread)
    asyncResult = pool.map_async(compute_profile, list_run)
    profile = asyncResult.get()
    return [i[0] for i in profile], [i[1] for i in profile]


def plot_DOPE_profile(list_template, list_model, list_model_files, sessionid,
                      pdb_codes, results, note=None):
    """Plot the dope result for each residue of each model
    """
    # Get color map
    color_map = cm.get_cmap('gist_rainbow')
    colors = [color_map(1.*i / len(list_model))
              for i in xrange(len(list_model))]
    for t in xrange(len(list_template)):
        # Plot the template and model profiles in the same plot for comparison:
        fig = plt.figure(figsize=(10, 7))
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.set_xlabel('Alignment position')
        ax1.set_ylabel('DOPE score per-residue')
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
            ax1.legend([temp_col, model_col[-1]], pdb_codes[t] +
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


def plot_partial_DOPE_profile(list_template, list_model, list_model_files,
                              sessionid, pdb_codes, results):
    """Divide the protein in two side and plot dope
    """
    # Plot templates
    for i in xrange(len(list_template)):
        middle = math.ceil(len(list_template[i]) / 2.0)
        part = [[0, int(middle)]]
        part += [[int(middle) + 1 , len(list_template[i])]]
        for e in part:
            # Plot models
            for j in xrange(len(list_model)):
                dist = abs(len(list_model[j]) - len(list_template[i]))
                if(len(list_model[j]) < len(list_template[i])):
                    list_model[j] += list_template[i][-dist :]
                elif(len(list_model[j]) > len(list_template[i])):
                    list_template[i] += list_model[j][-dist :]
                assert(len(list_model[j]) == len(list_template[i]))
                # Plot the template and model profiles
                # in the same plot for comparison:
                fig = plt.figure(figsize=(10, 7))
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
    """
    """
    delta = []
    for model_profile in list_model_profile:
        dist = abs(len(model_profile) - len(template_profile))
        if(len(model_profile) < len(template_profile)):
            model_profile += template_profile[-dist :]
        elif(len(model_profile) > len(template_profile)):
            template_profile += model_profile[-dist :]
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
    """
    # Plot models
    for i in xrange(len(list_delta_dope)):
        # Plot the template and model profiles in the same plot for comparison:
        fig = plt.figure(figsize=(10, 7))
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
    """
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]


def get_sequence(pdb_file):
    """
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
    """
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
                                output.write("{0}\t{1}\n".format(
                                    seq_pdb[num].split("_")[0],
                                    list_delta_dope[i][j]))
                            else:
                                output.write("{0}\t0\n".format(
                                    seq_pdb[num].split("_")[0]))
                            num += 1
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


def load_summary(summary_file):
    """Load summary data
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


def run_verification(conf_data, structure_check):
    """
    """
    raise NotImplemented
#    for soft in structure_check:
#        conf_data[soft]


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    args, parser = get_arguments()
    # Configure option
    conf_data = ModelingConfig(args.config, args.results)
    # Prepare modeling
    pdb_codes, pdb_files = get_pdb(conf_data, args.pdb, args.results)
    # Compute alignment
    if args.multifasta_file and not args.alignment_file:
        args.alignment_file = run_alignment(conf_data, args.multifasta_file,
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
    elif not args.multifasta_file and not args.alignment_file:
        print("You must provide an alignment or a multifasta sequence "
              "to continue...", file=sys.stderr)
        sys.exit(parser.print_help())
    if not args.model_name:
        args.model_name = get_model(args.alignment_file, pdb_codes)
    # Get session id
    sessionid = get_session_id(args.alignment_file)
    if MODELLER:
        # request verbose output
        log.level(output=1, notes=1, warnings=1, errors=1, memory=0)
        log.verbose()
        # Load environment
        env = get_environment(pdb_files)
    if MODELLER and "modeling" in args.list_operations:
        # Use several CPUs in a parallel job on this machine
        job_worker = get_parallel(int(args.thread))
        atm = compute_models(env, job_worker, args.alignment_file, pdb_codes,
                             pdb_files, args.model_name, args.model_quality,
                             args.number_model)
        save_general_data(atm, args.results, sessionid)
    if MODELLER and MATPLOTLIB and "profile" in args.list_operations:
        summary_file = (args.results + "modeller_summary_"
                        + str(sessionid) + ".csv")
        list_model = [args.model_name + ".B" + str(99990000 + i)
                      for i in xrange(1, args.number_model + 1)
                      if os.path.isfile(args.model_name + ".B"
                                        + str(99990000 + i) + ".pdb")]
        list_model_files = [args.results + i + ".pdb" for i in list_model]
        # Load alignment
        (list_template_profile,
         list_template_profile_smooth) = load_profiles(pdb_files,
                                                       pdb_codes,
                                                       args.alignment_file,
                                                       args.thread,
                                                       args.results)
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
    if args.structure_check:
        run_verification(None, args.structure_check)


if __name__ == '__main__':
    main()
