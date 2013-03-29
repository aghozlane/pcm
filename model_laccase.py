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

"""Run modelling."""

from __future__ import print_function
import argparse
import ConfigParser
import os
import sys
import urllib2
import subprocess
import re
import csv
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
    """Configure alignment and phylogeny.
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
        self.hdict["download_url"] = self.config.get("PDB_config",
                                                     "download_url")
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
        self.config.add_section("PDB_config")
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
                        "-outfile %output -output pir_aln -mode 3dcoffee "
                        "-n_core %proc -template_file %pdb".format(os.sep))
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
                        alignment_software="t_coffee", path_alignment="./",
                        list_operations=["modeling", "profile"],
                        add_heteroatom=0)
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
                        help="Indicate the number of hetero-atom residue that "
                        "should be added to the alignment software.")
    parser.add_argument('-p', '--pdb', type=str, required=True, nargs='+',
                        help="List of pdb files or codes.")
    parser.add_argument('-e', '--model_name', type=str,
                        help='Code of the sequence to modelize.')
    parser.add_argument('-g', '--path_alignment', type=isdir,
                        help='Path to t_coffee software.')
    parser.add_argument('-n', '--number_model', type=int,
                        help='Number of model to produce.')
    parser.add_argument('-q', '--model_quality', type=str,
                        choices=["fast", "normal", "max"],
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


# def get_pdb_code(pdb_file):
#     """Extract the PDB code from the header
#     """
#     pdb_code = ""
#     try:
#         with open(pdb_file, "rt") as pdb:
#             pdb_code = pdb.next().split()[4]
#     except IOError:
#         sys.exit("Error cannot open {0}".format(pdb_file))
#     assert(len(pdb_code) == 4)
#     return pdb_code


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
    print(path_soft)
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


def adjust_format(aln_pir_file, data_dict, pdb_codes, add_hetatm):
    """Correct pir format and add heteroatom in the alignment
    """
    output_file = (os.path.dirname(aln_pir_file) + os.sep +
                   os.path.basename(aln_pir_file).split(".")[0]
                   + "_corrected.pir")
    het_atm = ""
    structure_present = False
    if add_hetatm > 0:
        het_atm = "." * add_hetatm
    try:
        with open(output_file, "wt") as aln_pir:
            for element in data_dict:
                aln_pir.write(">P1;{0}\n".format(element))
                if data_dict[element][0] and not add_hetatm:
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
                    aln_pir.write("{0}:{1}:1 :A:{2}:A".format(
                                    type_data, element,
                                    len(sequence) + add_hetatm)
                                  + ": "*4 + "\n")
                end_aln_posit = data_dict[element][1].index("*")
                aln_pir.write("{0}{1}{2}\n".format(
                    data_dict[element][1][0, end_aln_posit],
                    het_atm, data_dict[element][1][end_aln_posit:]))
    except IOError:
        sys.exit("Error cannot open {0}".format())
    assert(structure_present)
    return output_file


def check_pir(aln_pir_file, pdb_codes, add_hetam):
    """Run checking of pir alignment file
    """
    try:
        with open(aln_pir_file, "rt") as aln_pir:
            aln_data = aln_pir.readlines()
        status, data_dict = check_format(aln_data, pdb_codes)
        if not status or add_hetam > 0:
            print("Try to correct pir alignment format.", file=sys.stderr)
            aln_pir_file = adjust_format(aln_pir_file, data_dict, pdb_codes,
                                         add_hetam)
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


def write_pir_file(aln_pir_file, data_fasta, pdb_codes, add_hetatm):
    """Write new pir alignment
    """
    hetatm = ""
    if add_hetatm > 0:
        hetatm = "." * add_hetatm
    try:
        with open(aln_pir_file, "wt") as aln_pir:
            for element in data_fasta:
                aln_pir.write(">P1;{0}\n".format(element))
                type_data = "sequence"
                if element in pdb_codes:
                    type_data = "structureX"
                sequence = data_fasta[element].replace("-", "")
                sequence = sequence.replace("\n", "")
                aln_pir.write("{0}:{1}:1 :A:{2}:A".format(
                                type_data, element,
                                len(sequence) + add_hetatm)
                                + ": "*4 + "\n")
                aln_pir.write("{0}{1}*\n".format(data_fasta[element], hetatm))
    except IOError:
        sys.exit("Error cannot open {0}".format(aln_pir_file))


def run_alignment(conf_data, multifasta_file, pdb_codes, pdb_files,
                  alignment_software, path_soft, add_hetatm, thread, results):
    """Compute alignment and adjust pir information
    """
    aln_pir_file = results + alignment_software + "_aln.pir"
    # compute
    if alignment_software in ("t_coffee", "clustalw2"):
        run_command(replace_motif(conf_data.hdict[alignment_software],
                                  path_soft, multifasta_file, pdb_files,
                                  aln_pir_file, thread))
    else:
        aln_fasta_file = results + alignment_software + "_aln.fasta"
        run_command(replace_motif(conf_data.hdict[alignment_software],
                                  path_soft, multifasta_file, pdb_files,
                                  aln_fasta_file, thread))
        data_fasta = get_fasta_data(aln_fasta_file)
        write_pir_file(aln_pir_file, data_fasta, pdb_codes, add_hetatm)
    if alignment_software in ("t_coffee", "clustalw2"):
        aln_pir_file = check_pir(aln_pir_file, pdb_codes, add_hetatm)
    return aln_pir_file


def get_model(aln_file, pdb_codes):
    """
    """
    regex = re.compile("^>\w+;([\w-]+)")
    try:
        with open(aln_file) as aln:
            for line in aln:
                match = regex.match(line)
                if match:
                    if match.group(1) not in pdb_codes:
                        model = match.group(1)
                        break
    except IOError:
        sys.exit("Error cannot open {0}".format(aln_file))
    return model


def get_environment(pdb_files):
    """Set Modeller environment parameters
    """
    env = environ()
    # Set PDB directory
    env.io.atom_files_directory = [os.path.dirname(pdb) for pdb in pdb_files]
    # Read in HETATM records from template PDBs
    env.io.hetatm = True
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
    atm.starting_model = len(pdb_codes)
    atm.ending_model = int(number_model)
    # Start slave
    atm.use_parallel_job(job_worker)
    if model_quality == "fast":
        print("Start fast modeling")
        atm.very_fast()
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


def barplot(data, label, result_data):
    """Barplot general method
    """
    plt.xlabel(label[0])
    plt.ylabel(label[1])
    plt.title(label[2])
    plt.bar(range(0, len(data)), data)
    plt.savefig(result_data)
    plt.clf()


def summary_data(data, results):
    """
    """
    output_file = results + "modeller_summary.csv"
    try:
        with open(output_file, "wt") as out:
            writer = csv.writer(out, delimiter='\t')
            writer.writerow(["molpdf", "DOPE score", "Normalized DOPE score",
                             "GA341 score"])
            for i in xrange(len(data)):
                writer.writerow([(data[i]["name"].split(".")[0]
                                  + "_{0}".format(i + 1)), data[i]["molpdf"],
                                 data[i]["DOPE score"],
                                 data[i]["Normalized DOPE score"],
                                 data[i]["GA341 score"][0]])
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


def save_general_data(atm, results):
    """Write general parameters of the produced models
    """
    # Get converged models
    ok_models = filter(lambda x: x['failure'] is None, atm.outputs)
    # Rank the models by molpdf score
    ok_models.sort(lambda a, b: cmp(a["DOPE score"], b["DOPE score"]))
    # Barplot Dope score of each ok models
    if MATPLOTLIB:
        barplot([i["DOPE score"] for i in ok_models],
                ["Model", "Dope score", "Dope score per model"],
                results + "dope_per_model.svg")
    # Write data summary
    summary_data(ok_models, results)


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
    # Insert gaps into the profile corresponding to those in seq:
    for n, res in enumerate(seq.residues):
        for gap in xrange(res.get_leading_gaps()):
            vals.insert(n, None)
    # Add a gap at position '0', so that we effectively count from 1:
    vals.insert(0, None)
    return vals


def compute_profile(data):
    """Compute the profile of each model
    """
    env = get_environment([data[0]])
    aln = alignment(env, file=data[2])
    # env and pdbfile
    pdb = complete_pdb(env, data[0])
    # all atom selection
    s = selection(pdb)
    # profile result
    s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=data[1],
                  normalize_profile=False, smoothing_window=15)
    return get_profile(data[1], aln[os.path.basename(data[0]).split(".")[0]])



def load_profiles(env, pdb_files, pdb_codes, alignment_file, thread, results):
    """Compute and load the profile of each model
    """
    list_run = [[ pdb_files[i], (results + pdb_codes[i] + '.profile'),
                 alignment_file]
                for i in xrange(len(pdb_files))]
    pool = mp.Pool(processes=thread)
    asyncResult = pool.map_async(compute_profile, list_run)
    profile = asyncResult.get()
    return profile


def plot_DOPE_profile(list_template, list_model, list_model_files, results):
    """Plot the dope result for each residue of each model
    """
    # Get color map
    color_map = cm.get_cmap('gist_rainbow')
    colors = [color_map(1.*i / len(list_model))
              for i in xrange(len(list_model))]
    # Plot the template and model profiles in the same plot for comparison:
    fig = plt.figure(figsize=(10, 7))
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.set_xlabel('Alignment position')
    ax1.set_ylabel('DOPE per-residue score')
    ax1.set_xlim(0, len(list_model[0]))
    # Plot templates
    for template in list_template:
        ax1.plot(template, color="green", linewidth=1, label='Template')
    # Plot models
    for i in xrange(len(list_model)):
        ax1.plot(list_model[i], color=colors[i], linewidth=1, label='Model')
    ax1.legend(["template"] + [(list_model_files[i].split(".")[0]
                                + "_{0}".format(i + 1))
                               for i in xrange(len(list_model_files))],
               loc="upper center", numpoints=1, bbox_to_anchor=(0.5, 1.12),
               ncol=3, fancybox=True, shadow=True)
    plt.savefig(results + os.sep + 'dope_profile.svg')


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
        save_general_data(atm, args.results)
    if MODELLER and MATPLOTLIB and "profile" in args.list_operations:
        list_model = [args.model_name + ".B" + str(99990000 + i)
                      for i in xrange(1, args.number_model + 1)]
        list_model_files = [i + ".pdb" for i in list_model]
        # Load alignment
        list_template = load_profiles(env, pdb_files, pdb_codes,
                                      args.alignment_file, args.thread,
                                      args.results)
        list_model = load_profiles(env, list_model_files, list_model,
                                   args.alignment_file, args.thread,
                                   args.results)
        # Plot DOPE profile
        plot_DOPE_profile(list_template, list_model, list_model_files,
                          args.results)
    if args.structure_check:
        run_verification(None, args.structure_check)


if __name__ == '__main__':
    main()
