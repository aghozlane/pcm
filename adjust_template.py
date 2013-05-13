#!/usr/bin/env python
import sys
import os
import argparse

seqdict = {"ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D",
           "CYS":"C", "GLU":"E", "GLN":"Q", "GLY":"G",
           "HIS":"H", "ILE":"I", "LEU":"L", "LYS":"K",
           "MET":"M", "PHE":"F", "PRO":"P", "SER":"S",
           "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V"}


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} file.pdb".format(sys.argv[0]))
    parser.set_defaults(activity="renum")
    parser.add_argument("-a", "--activity", type=str,
                        choices=["renum", "sequence"],
                        help='Select operations.')
    parser.add_argument('-p', '--pdb', type=str, required=True, nargs='+',
                        help="List of pdb files or codes.")
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


def renum_pdb(pdb_file, activity):
    """
    """
    num = 0
    try:
        with open(pdb_file, "rt") as pdb:
            for line in pdb:
                aa = line[17:20]
                atom = line[13:16]
                field = line[0:4]
                if(atom == "N  "):
                    num += 1
                    if activity == "sequence" and field == "ATOM":
                        sys.stdout.write(seqdict[aa])
                pdb_line = list(line)
                if(field == "ATOM"):
                    pdb_line[23:26] = get_numstring(num, 3)
                pdb_line = "".join(pdb_line)
                if activity == "renum":
                    sys.stdout.write(pdb_line)
                # print(line[23:26])
            sys.stdout.write("\n")
    except IOError:
        sys.exit("Error cannot open {0}".format(pdb_file))


def main():
    """
    """
    args = get_arguments()
    try:
        for pdb in args.pdb:
            if pdb.endswith(".pdb") and os.path.isfile(pdb):
                if args.activity == "sequence":
                    print(">{0}"
                    .format(".".join(os.path.basename(pdb).split(".")[:-1])))
                renum_pdb(pdb, args.activity)
            else:
                raise ValueError
    except ValueError:
        sys.exit("Please indicate pdb file as entry (*.pdb)")

if __name__ == '__main__':
    main()
