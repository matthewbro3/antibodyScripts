#!/usr/bin/env python3
import sys, math

threeToOne = dict(zip(['ALA','ARG','ASN','ASP','CYS','GLN','GLY','GLU','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'],
                      ['A',  'R',  'N',  'D',  'C',  'Q',  'G',  'E',  'H',  'I',  'L',  'K',  'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V']))

oneToThree = dict(zip(['A',  'R',  'N',  'D',  'C',  'Q',  'G',  'E',  'H',  'I',  'L',  'K',  'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V'],
                      ['ALA','ARG','ASN','ASP','CYS','GLN','GLY','GLU','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']))

gCharge = dict(zip(['ALA','ARG','ASN','ASP','CYS','GLN','GLY','GLU','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'],
                   [ 0,  1,  0, -1,  0,  0,  0,  -1, 0.5,0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0]))

resID = []
xCoord = []
yCoord = []
zCoord = []

def pdbH3Ca(file):
    """
    Strips out the identities and alpha carbon coordinates
    of the residues in the CDR-H3 from a .pdb file
    and stores them as global lists.

    To check the extracted data for correctness, use following test code:
    
    for x in range(len(resID)):
        print((resID[x], xCoord[x], yCoord[x], zCoord[x]))  
    """
    with open(file) as f:
        inRange = 0
        global resID
        global xCoord
        global yCoord
        global zCoord
        for line in f.readlines():
            atom = line.split()
            if atom[0] == 'ATOM':
                if atom[5] == '104':
                    inRange =  0
                elif atom[5] == '95':
                    inRange = 1
                if(atom[2] == 'CA' and atom[4] == 'H' and inRange == 1):
                    resID.append(atom[3])
                    xCoord.append(float(atom[6]))
                    yCoord.append(float(atom[7]))
                    zCoord.append(float(atom[8]))
        f.close()

def ccdCalc(file):
    ccd = 0
    beta = 0.5
    for i in range(len(resID)):
        for j in range(i):
            qi = gCharge[resID[i]]
            qj = gCharge[resID[j]]
            distance = math.sqrt((xCoord[i]-xCoord[j])**2+(yCoord[i]-yCoord[j])**2+(zCoord[i]-zCoord[j])**2)
            ccd = ccd+qj*qi*(distance)**beta
    return ccd    

def chdCalc(file):
    with open('consensus.hpb') as c:
        hpTable = {}
        for line in c.readlines():
            entry = line.split()
            ID = entry[0]
            value = float(entry[1])
            hpTable[ID] = value
        chd = 0
        beta = 0.5
        for i in range(len(resID)):
            for j in range(i):
                qi = hpTable[resID[i]]
                qj = hpTable[resID[j]]
                distance = math.sqrt((xCoord[i]-xCoord[j])**2+(yCoord[i]-yCoord[j])**2+(zCoord[i]-zCoord[j])**2)
                chd = chd+qj*qi*(distance)**beta
        c.close()
    return chd

def cccdCalc(file):
    """
    Strips out the coordinates of charged side chains
    and calculates a cccd value based on those coordinates.

    To check the extracted data for correctness, use following test code:
    
    for x in range(len(resID)):
        print((resID[x], xCoord[x], yCoord[x], zCoord[x]))  
    """

    with open(file) as f:
        inRange = 0
        cresID = []
        cxCoord = []
        cyCoord = []
        czCoord = []
        for line in f.readlines():
            atom = line.split()
            if atom[0] == 'ATOM':
                if atom[5] == '104':
                    inRange =  0
                elif atom[5] == '95':
                    inRange = 1
                if(inRange == 1 and atom[4] == 'H'):
                    if((atom[3] == 'ASP' and atom[2] == 'CG')or
                       (atom[3] == 'ARG' and atom[2] == 'NH1') or
                       (atom[3] == 'GLU' and atom[2] == 'CD') or
                       (atom[3] == 'LYS' and atom[2] == 'NZ') or
                       (atom[3] == 'HIS' and atom[2] == 'NE2')):
                        cresID.append(atom[3])
                        cxCoord.append(float(atom[6]))
                        cyCoord.append(float(atom[7]))
                        czCoord.append(float(atom[8]))
        cccd = 0
        beta = 0.5
        for i in range(len(cresID)):
            for j in range(i):
                qi = gCharge[cresID[i]]
                qj = gCharge[cresID[j]]
                distance = math.sqrt((cxCoord[i]-cxCoord[j])**2+(cyCoord[i]-cyCoord[j])**2+(czCoord[i]-czCoord[j])**2)
                cccd = cccd+qj*qi*(distance)**beta
        return cccd  
        f.close()


target = sys.argv[1]
pdbH3Ca(target) 
ccd = ccdCalc(target)
chd = chdCalc(target)
cccd = cccdCalc(target)
print((ccd,chd,cccd))
for x in range(len(resID)):
    print((resID[x], xCoord[x], yCoord[x], zCoord[x]))  
