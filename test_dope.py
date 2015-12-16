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
import glob
import os
import sys


try:
    from modeller import *
    from modeller.automodel import *
    from modeller.scripts import complete_pdb
    MODELLER = True
except ImportError:
    MODELLER = False
    print("Could not import modeller\nNo modeling will be available",
          file=sys.stderr)


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
    parser.set_defaults(results=local_path)
    parser.add_argument('-p', '--pdb_file', type=isfile,
                        help='PDB files.')
    parser.add_argument('-d', '--pdb_dir', type=isdir,
                        help='Directory that contains pdb files.')
#     parser.add_argument('-i', '--alignment_file', type=isfile, required=True,
#                         help='Alignment in pir format.')
    parser.add_argument('-r', '--results', type=isdir,
                        help='Path to the result directory.')
    return parser.parse_args(), parser


def check_directory(pdb_dir):
    """Get the list of files in the directory
      Arguments:
          pdb_dir: Path to the directory that contains PDB files
      Returns: List of PDB files.
    """
    return glob.glob('{0}{1}*.pdb'.format(pdb_dir, os.sep))


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


def get_profile(profile_file):
    """Read `profile_file` into a Python array, and add gaps corresponding to
       the alignment sequence `seq`."""
    # Read all non-comment and non-blank lines from the file:
    vals = []
    try:
        with open(profile_file, "rt") as f:
            for line in f:
                if not line.startswith('#') and len(line) > 10:
                    spl = line.split()
#                     print(spl[-1])
#                     print(float(spl[-1]))
                    vals.append(float(spl[-1]))
    except IOError:
        sys.exit("Error cannot open {0}".format(profile_file))
    except ValueError:
        sys.exit("Error with the conversion to float of {0}".format(line))
    return vals


def compute_profile(pdb_file, profile_file):
    """Compute the profile of each model
    """
    env = get_environment([pdb_file])
    # env and pdbfile
    pdb = complete_pdb(env, pdb_file)
    # all atom selection
    s = selection(pdb)
    # profile result
    normal_result = s.assess_dope()
    # output='ENERGY_PROFILE NO_REPORT',
    # normalize_profile=False
    # get_profile(profile_file)
    profile = s.get_dope_profile()
#     profile = s.get_dopehr_profile()
#     profile.write_to_file(profile_file)
#     val = []
#     with open(profile_file, "rt") as test:
#         for i in test:
# #             print(i.split())
#             val += [float(i.split()[1])]
#     for i in profile:
#         print(i)
    return [i.energy for i in profile], normal_result


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    list_pdbfiles = []
    profile = []
    args, parser = get_arguments()
    log.level(output=0, notes=0, warnings=0, errors=0, memory=0)
    # Get pdb
    try:
        if args.pdb_dir:
            list_pdbfiles += check_directory(args.pdb_dir)
        if args.pdb_file:
            list_pdbfiles += [args.pdb_file]
        assert(list_pdbfiles != [])
    except ValueError:
        print("Please, indicate a PDB file or a directory "
        "that contains PDB files.")
        sys.exit(parser.print_help())
    for i in xrange(len(list_pdbfiles)):
        profile_file = (args.results +
                        ".".join(os.path.basename(
                            list_pdbfiles[i]).split(".")[:-1]) +
                        '.profile')
        profile_temp, result = compute_profile(list_pdbfiles[i], profile_file)
        profile = [profile_temp]
    if(len(profile) > 0):
        for i in xrange(len(list_pdbfiles)):
            print("{0}: {1} {2}".format(os.path.basename(list_pdbfiles[i]),
                                    sum(profile[i]), result))


if __name__ == '__main__':
    main()
