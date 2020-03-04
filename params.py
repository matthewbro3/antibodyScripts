#!/usr/bin/env python3
import sys

"""
import re
"""

threeToOne = dict(zip(['ALA','ARG','ASN','ASP','CYS','GLN','GLY','GLU','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'],
                      ['A',  'R',  'N',  'D',  'C',  'Q',  'G',  'E',  'H',  'I',  'L',  'K',  'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V']))

oneToThree = dict(zip(['A',  'R',  'N',  'D',  'C',  'Q',  'G',  'E',  'H',  'I',  'L',  'K',  'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V'],
                      ['ALA','ARG','ASN','ASP','CYS','GLN','GLY','GLU','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']))

gCharge=dict(zip(['A','R','N','D','C','Q','G','E','H','I','L','K','M','F','P','S','T','W','Y','V'],
                 [ 0,  1,  0, -1,  0,  0,  0,  -1, 0.5,0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0]))

def readH3Sequence(file):
    """ Reads a .seq file containing resid/1-letter aa codes (space separated)
        and identifies the residues between H95 and H102 assembling their
        1-letter codes into a string
    """
    with open(file) as f:
        h3 = ''
        inRange = 0
        protein = f.readlines()
        for line in protein:
            idres = line.split()
            if idres[0] == 'H95':
                inRange = 1
            elif idres[0] == 'H103':
                inRange = 0

            if inRange == 1:
                h3 = h3 + idres[1]

        f.close()
        return h3

def chargeCalc(chain):
    """ 
    Takes a sequence as a string using 1-letter code and
    calculates the net charge. His is given a charge of +0.5

    H3 = "ACDEFGHIKLMNPQRSTVWY"
    charge = chargeCalc(H3)
    print(charge)
    
    Test code should return 0.5
    """
    netCharge = 0
    for residue in chain:
        netCharge = netCharge + gCharge[residue]
    return netCharge

def scdCalc(chain):
    """ 
    Takes a sequences as a string using 1-letter code
    and calculates the sequence charge decoration.
    Histidine given a charge of +0.5.
    Beta value for SCD equation set to 1/2; change as needed.

    H3 = "ACDEFGHIKLMNPQRSTVWY"
    scd = scdCalc(H3)
    print(scd)

    Test code should return 2.45
    """
    scd = 0
    beta = 0.5
    for i in range(len(chain)):
        for j in range(i):
            qi = gCharge[chain[i]]
            qj = gCharge[chain[j]]
            scd = scd+qj*qi*(i-j)**beta
    return scd


def hpCalc(chain):
    """
    Takes a sequence as a string using 1-letter code and 
    calculates the net hydrophobicity.
    Requires a table of hydrophobicity values, given as consensus.hpb

    H3 = "ACDEFGHIKLMNPQRSTVWY"
    hp = hpCalc(H3)
    print(hp)

    Test code should return 2.97
    """

    with open('consensus.hpb') as f:
        hpTable = {}
        for line in f.readlines():
            entry = line.split()
            ID = entry[0]
            value = float(entry[1])
            hpTable[ID] = value
        netHP = 0
        for residue in chain:
            code = oneToThree[residue]
            netHP = netHP + hpTable[code]
        f.close()
        return netHP

def shdCalc(chain):
    """
    Takes a sequence as a string using 1-letter code and
    calculates the sequence hydrophobicity pattern (SCD).
    Requires a table of hydrophobicity values, given as consensus.hpb
    Beta value for SCD equation set to 1/2; change as needed.
    
    
    H3 = "ACDEFGHIKLMNPQRSTVWY"
    shd = shdCalc(H3)
    print(shd)

    Test code should return 0.0074
    """
    with open('consensus.hpb') as f:
        hpTable = {}
        for line in f.readlines():
            entry = line.split()
            ID = entry[0]
            value = float(entry[1])
            hpTable[ID] = value
        shd = 0
        beta = 0.5
        for i in range(len(chain)):
            for j in range(i):
                qi = hpTable[oneToThree[chain[i]]]
                qj = hpTable[oneToThree[chain[j]]]
                scd = shd+qj*qi*(i-j)**beta
        f.close()
        return scd
    
            
            
# Use following line to test program

#h3 = "ACDEFGHIKLMNPQRSTVWY"

# Main program
h3     = readH3Sequence(sys.argv[1])
charge = chargeCalc(h3)
scd    = scdCalc(h3)
hp     = hpCalc(h3)
shd    = shdCalc(h3)
print(h3, " charge: ", charge, " scd: ", scd, " hp: ", hp, " shd: ", shd)

