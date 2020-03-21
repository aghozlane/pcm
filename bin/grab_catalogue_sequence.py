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

"""Extract a list of sequence from the catalogue."""


from __future__ import print_function
import argparse
import sys
import os
import csv
import bisect
import textwrap
import collections

__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2014, INRA"
__credits__ = ["Amine Ghozlane"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@jouy.inra.fr"
__status__ = "Developpement"


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


def getArguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    parser.set_defaults(results=".{0}".format(os.sep))
    parser.add_argument('-i', dest='list_sequences_file', type=isfile,
                        help='List of sequence to extract (in a file).')
    parser.add_argument('-l', dest='list_sequences', type=str,
                        help='List of sequence to extract (in a file).')
    parser.add_argument('-s', dest='set_sequences', type=str, nargs='+',
                        help="Directory that contains PDB files.")
    parser.add_argument('-d', dest='catalogue_file', type=isfile,
                        required=True, help='Database query.')
    parser.add_argument('-n', dest='not_in_database', action='store_true',
                        help='Select instead elements which are not in the'
                        ' list.')
    parser.add_argument('-o', dest='output_file', type=str,
                        help='Output file.')
    parser.add_argument('-r', dest='results', type=isdir,
                        help='Path to result directory.')
    return parser.parse_args()


def extract_interest_elements(list_sequences_file):
    """Get a list of the element of interest
    """
    list_sequences = []
    try:
        with open(list_sequences_file, "rt") as list_seq:
            list_sequences_reader = csv.reader(list_seq, delimiter="\t")
            for line in list_sequences_reader:
                 list_sequences.append(line[0])
    except IOError:
        sys.exit("Error cannot the file : {0}".format(list_sequences_file))
    return list_sequences


# def is_selected(header, list_sequences):
#     """
#     """
#     for element in list_sequences:
#         if element in header:
#             return element
#     return None


def get_element(name, input_list):
    """Search name in input list
      Arguments:
        input_list: List
        name: Search criteria
    """
    # Searching the node with its name
    i = bisect.bisect_left(input_list, name)
    # Object has been found
    if(i != len(input_list) and input_list[i] == name):
        return True #input_list[i]
    return False #None

def get_unique(seq):
    """Get unique elements with order preserving
      :Parameters:
       seq: List of elements
    """
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]

def extract_catalogue_sequence(list_sequences, catalogue_file, not_in_database):
    """
    """
    grab_sequence = False
    interest_sequence = {}
    title = ""
    try:
        with open(catalogue_file, "rt") as catalogue:
            for line in catalogue:
                if line.startswith(">"):
                    grab_sequence = False
                    title = line[1:].replace("\n", "").replace("\r", "")
                    if " " in title:
                        title = title.split(" ")[0]
                    selection = get_element(title, list_sequences)
                    if selection and not not_in_database:
                        interest_sequence[title] = ""
                        grab_sequence = True
                    elif not selection and not_in_database:
                        interest_sequence[title] = ""
                        grab_sequence = True
                elif grab_sequence and len(line) > 0:
                    interest_sequence[title] += line.replace("\n", "").replace("\r", "")
            assert(len(interest_sequence) > 0)
    except IOError:
        sys.exit("Error cannot the file : {0}".format(catalogue_file))
    except AssertionError:
        sys.exit("Error no element detected in the file : {0}"
                 .format(catalogue_file))
    return interest_sequence


def fill(text, width=80):
    """Split text"""
    return os.linesep.join(text[i:i+width] for i in xrange(0, len(text), width))


def write_interest_sequence(interest_sequence, output_file):
    """
    """
    try:
        with open(output_file, "wt") as output:
            for key in interest_sequence:
                output.write(">{1}{0}{2}{0}".format(
                                os.linesep, key, fill(interest_sequence[key])))
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get the arguments
    args = getArguments()
    # Get List of sequence of interest
    print("Load the list of sequence of interest ...")
    try:
        if args.list_sequences_file:
            list_sequences = extract_interest_elements(args.list_sequences_file)
        elif args.set_sequences:
            list_sequences = args.set_sequences.sort()
        elif args.list_sequences:
            if args.list_sequences.startswith("[") and args.list_sequences.endswith("]"):
                args.list_sequences = args.list_sequences.replace("[", "").replace("]", "")
                list_sequences = args.list_sequences.split(",")
            else:
                sys.exit("The list should be in the format [name1,name2]")
        else:
            sys.exit("Please provide a file containing the list of query sequences or a set of sequences")
        assert(list_sequences > 0)
        list_sequences.sort()
        list_sequences = get_unique(list_sequences)
    except AssertionError:
        sys.exit("Error no element detected in the query list")
    
    print("{0} (unique) sequences to search".format(len(list_sequences)))
    # Extract catalogue sequence
    print("Extract sequences from the catalogue...")
    interest_sequence = extract_catalogue_sequence(list_sequences,
                                                   args.catalogue_file,
                                                   args.not_in_database)
    print("{0} extracted sequences".format(len(interest_sequence)))
    # Write sequences
    if not args.output_file:
       args.output_file = "extracted_sequence.fasta"
    if not args.not_in_database:
        if len(list_sequences) != len(interest_sequence.keys()):          
            missing_elements = list(set(list_sequences) - set(interest_sequence.keys()))

            if len(missing_elements) > 0:
                print("There is {0} missing elements.".format(len(missing_elements)))
                for miss in missing_elements:
                    print(miss, file=sys.stderr)
            else:
                duplicated_elements = [item for item, count in collections.Counter(list_sequences).items() if count > 1]
                if len(duplicated_elements) > 0:
                    print("There is duplicated elements in the input list:")
                    print(duplicated_elements)
    print("Write sequences to {0}".format(args.output_file))
    write_interest_sequence(interest_sequence, args.output_file)
    print("Done.")


if __name__ == '__main__':
    main()
