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
    parser.add_argument("-a", dest="activity", type=str,
                        choices=["renum", "sequence"],
                        help='Select operations.')
    parser.add_argument('-p', dest='pdb', type=str, required=True, nargs='+',
                        help="List of pdb files or codes.")
    # parser.add_argument('-s', dest='start_position', type=int, default=1,
    #                    help="List of pdb files or codes (default : 1).")
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


# def renum_pdb(pdb_file, activity, start_position):
def renum_pdb(pdb_file, activity):
    """
    """
    res = 0
    # num = start_position - 1
    num = 0
    aa_prev = ""
    chain_prev = ""
    try:
        with open(pdb_file, "rt") as pdb:
            for line in pdb:
                chain = line[20:22]
                aa = line[17:20]
                newres = line[22:26]
                field = line[0:4]
                if chain != chain_prev and field == "ATOM":
                    num = int(newres) - 1
                    chain_prev = chain
                if((newres != res or aa != aa_prev) and field == "ATOM"):
                    aa_prev = aa
                    res = newres
                    num += 1
                    if activity == "sequence":
                        sys.stdout.write(seqdict[aa])
                pdb_line = list(line)
                if(field == "ATOM"):
                    pdb_line[22:26] = get_numstring(num, 4)
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
                # renum_pdb(pdb, args.activity, args.start_position)
                renum_pdb(pdb, args.activity)
            else:
                raise ValueError
    except ValueError:
        sys.exit("Please indicate pdb file as entry (*.pdb)")

if __name__ == '__main__':
    main()
