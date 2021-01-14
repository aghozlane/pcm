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

"""Convert names"""

import os
import argparse
import sys

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
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='input_file', type=isfile, required=True,
                        help="Input file")
    parser.add_argument('-o', dest='output_file', 
                        default=os.curdir + os.sep + "output.fasta", 
                        help="Output file")
    parser.add_argument('-a', dest='association_file',
                        default=os.curdir + os.sep + "association.tsv", 
                        help="Association file")
    return parser.parse_args()


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    args = get_arguments()
    count_prot = 1
    seq = ""
    header = ""
    with open(args.input_file, "rt") as input_file, open(args.output_file, "wt") as output_file, open(args.association_file, "wt") as association_file:
        for line in input_file:
            if line.startswith(">"):
                if len(seq) > 0:
                    output_file.write(">protein{}\n{}\n".format(count_prot, fill(seq)))
                    association_file.write("protein{}\t{}\n".format(count_prot, header))
                    seq = ""
                header = line[1:].strip().replace("\n", "")
                count_prot += 1
            else:
                seq += line.strip().replace("\n", "")
        if len(seq) > 0:
            output_file.write(">protein{}\n{}\n".format(count_prot, fill(seq)))
            association_file.write("protein{}\t{}\n".format(count_prot, header))



if __name__ == '__main__':
    main()