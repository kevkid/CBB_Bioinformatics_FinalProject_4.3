#!/usr/bin/python

__author__ = "Peter Williams"
__copyright__ = "Copyright 2016"
__credits__ = ["Peter Williams"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Peter Williams"
__email__ = "peter.williams@yale.edu"

### Usage:      python3 dihedrals.py -i <input file> -a <Atom ID> -b <Atom ID> -c <Atom ID> -d <Atom ID> -o <output file>
### Examples:   python3 dihedrals.py -i 1g8p.pdb -a 1 -b 12 -c 3 -d 2 -o output.txt
###				python3 dihedrals.py -i 1g8p.pdb -a 1 -b 12 -c 3 -d 2 
### Note:       	Calculates the dihedral angle of a set of 4 atoms in a pdb
###						with corresponding Atom IDs
###
### Input Formats:	-i 			pdb file
###					-a,-b,-c	integers for atom ids
### Output Format:	txt file with statement indicating the dihedral angle

### Import libraries
import argparse
import numpy as np

### This is one way to read in arguments in Python. We need to read input file and score file.
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file name', required=True)
parser.add_argument('-a', '--AtomID1', help='Atom ID of one of 4 atoms whose dihedral angle we wish to calculate', required=True)
parser.add_argument('-b', '--AtomID2', help='Atom ID of one of 4 atoms whose dihedral angle we wish to calculate', required=True)
parser.add_argument('-c', '--AtomID3', help='Atom ID of one of 4 atoms whose dihedral angle we wish to calculate', required=True)
parser.add_argument('-d', '--AtomID4', help='Atom ID of one of 4 atoms whose dihedral angle we wish to calculate', required=True)
parser.add_argument('-o', '--output', nargs='?',default='output.txt', help='Output file name')
args = parser.parse_args()

### Implementation
def dihedral(filename,id1,id2,id3,id4,outfile):
	#######################################################################
	### Initialize variables
	ids=[int(id1),int(id2),int(id3),int(id4)]
	ids=np.sort(ids)
	Rs=np.zeros([4,3])
	#######################################################################
	### Collect coordinates of four atoms of interest
	f=open(filename,'r')
	for x in f:
		# If the line in the pdb file contains coordinates it starts with Atom and has the following form
		# x[0:4] ="Atom"
		# x[6:11] = Atom serial number
		# x[12:16] = Atom Name
		# x[16] = Alternate location indicator
		# x[17:20] = Residue Name
		# x[21] = Chain Identifier
		# x[22:26] = Residue Sequence number
		# x[26] = Code for insertions of residues
		# x[30:38] = x orthogonal angstrom coordinate
		# x[38:46] = y orthogonal angstrom coordinate
		# x[46:54] = z orthogonal angstrom coordinate
		
		if x[0:4]=="ATOM":
			if int(x[6:12])==ids[0]:
				Rs[0,0]=float(x[30:38])
				Rs[0,1]=float(x[38:46])
				Rs[0,2]=float(x[46:54])
			if int(x[6:12])==ids[1]:
				Rs[1,0]=float(x[30:38])
				Rs[1,1]=float(x[38:46])
				Rs[1,2]=float(x[46:54])
			if int(x[6:12])==ids[2]:
				Rs[2,0]=float(x[30:38])
				Rs[2,1]=float(x[38:46])
				Rs[2,2]=float(x[46:54])
			if int(x[6:12])==ids[3]:
				Rs[3,0]=float(x[30:38])
				Rs[3,1]=float(x[38:46])
				Rs[3,2]=float(x[46:54])
				break;
	f.close()
	
	#####################################################################
	### Calculate phi angle
	
	# Creat matrix of displacement vectors of consecutive atoms
	Rij=Rs[1:,:]-Rs[0:-1,:]
	# Use formulae for cos and sin of dihedral angles
	cosphi=np.dot(np.cross(Rij[0,:],-Rij[1,:]),np.cross(Rij[1,:],-Rij[2,:]))/(np.linalg.norm(np.cross(Rij[0,:],Rij[1,:]))*np.linalg.norm(np.cross(Rij[1,:],Rij[2,:])))
	sinphi=np.dot(Rij[1,:],np.cross(np.cross(Rij[0,:],-Rij[1,:]),np.cross(Rij[1,:],-Rij[2,:])))/(np.linalg.norm(Rij[1,:])*np.linalg.norm(np.cross(Rij[0,:],Rij[1,:]))*np.linalg.norm(np.cross(Rij[1,:],Rij[2,:])))
	# Determine quadrant of angle based on signs of cos and sin and 
	#    place angle accordingly and convert to degrees
	if cosphi>0:
		phi=np.arcsin(sinphi)*180/np.pi
	elif sinphi<0:
		phi=-np.arccos(cosphi)*180/np.pi
	else:
		phi=np.arccos(cosphi)*180/np.pi
	
	#####################################################################
	### Output result into output text file
	
	# Open outfile to set up for writing
	f=open(outfile,'w')
	f.write("The dihedral angle formed by atoms with atomIDs [%d, %d, %d, %d] = %.3f" %(ids[0],ids[1],ids[2],ids[3],float(phi)))
	f.write(u"\u00b0")
	f.write('\n')
	f.close()
	
### Run
dihedral(args.input, args.AtomID1, args.AtomID2, args.AtomID3, args.AtomID4, args.output)

