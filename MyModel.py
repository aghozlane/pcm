from modeller import *
from modeller.automodel import *
import sys, re

class RestraintModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        for struct in self.psipred_result:
            # Constrain all residues to be alpha-helical
            if(struct[0] == "H"):
                rsr.add(secondary_structure.alpha(
                    self.residue_range("{0}:".format(struct[1]),
                                       "{0}:".format(struct[2]))))
            # Constrain all residues to be beta-strand
            elif(struct[0] == "E"):
                rsr.add(secondary_structure.strand(
                    self.residue_range("{0}:".format(struct[1]),
                                       "{0}:".format(struct[2]))))

