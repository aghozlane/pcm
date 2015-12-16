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

"""Extract fasta sequence from multifasta."""

from __future__ import print_function
import sys
import os

def extract_data(output_dir, input_file):
  """
  """
  output = None
  file_open = False
  try:
    with open(input_file, "rt") as input:
      for line in input:
	if line[0] == ">":
	  # Close old file
	  if output:
	    output.close()
	    output = None
	    file_open = False
	  # Open new file
	  output = open(output_dir + os.sep + line[1:].strip() + ".fasta", "wt")
	  output.write(line)
	  file_open = True
	elif file_open and len(line.strip()) > 0:
	  output.write(line)
      
  except IOError:
    sys.exit("Error cannot open")

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    if len(sys.argv) != 2:
      sys.exit("Error, the program take only one argument : ./{0} file.fasta".format(sys.argv[0]))
    output_dir = ".".join(os.path.basename(sys.argv[1]).split(".")[:-1])
    try:
      os.mkdir(output_dir)
    except:
      print("Failed to create : {0}".format(output_dir), file=sys.stderr)
    extract_data(output_dir, sys.argv[1])
    

if __name__ == '__main__':
    main()

