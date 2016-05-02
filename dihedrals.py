#!/usr/bin/python

__author__ = "Peter Williams"
__copyright__ = "Copyright 2016"
__credits__ = ["Peter Williams"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Peter Williams"
__email__ = "peter.williams@yale.edu"

### Usage:      python3 dihedrals.py -i <input file> -a <Atom ID> -b <Atom ID> -c <Atom ID> -d <Atom ID> -o <output file>
### Examples:   python3 dihedrals.py -i sample-input.pdb -n sample-ids.txt -o sample-output.txt
###				python3 dihedrals.py -i sample-input.pdb -n sample-ids.txt
### Note:       	Calculates the dihedral angle of a set of 4 atoms in a pdb
###						with corresponding Atom IDs
###					The -o flag is optional
###
### Input Formats:	-i 	string of corresponding pdb file
###					-n	string corresponding to text file containing
###							integers for atom ids. Each line is a set 
###							of 4 which defines the dihedral angle.
###					-o	string containing the name of the file to which the angles are output
### Output Format:	txt file that is tab deliminated indicating the atom names, 
###						the atom ids, and the angle in degrees

### Import libraries
import argparse
import numpy as np

### This is one way to read in arguments in Python. We need to read input file and score file.
parser = argparse.ArgumentParser(description='Dihedral Angle Calculator')
parser.add_argument('-i', '--input', help='input file name', required=True)
parser.add_argument('-n', '--inputids', help='id file name', required=True)
parser.add_argument('-o', '--output', nargs='?',default='output.txt', help='Output file name')
args = parser.parse_args()

### Implementation
def dihedral(filename,idfile,outfile):
	###########################################################################
	### Load textfile of ids
	idarray=np.loadtxt(idfile)
	
	###########################################################################
	### Initialize outputfile
		
	# Open outfile to set up for writing
	g=open(outfile,'w')
	g.write('Atom Names')
	g.write('\t\t')
	g.write('Atom IDs')
	g.write('\t\t')
	g.write('Angle (')
	g.write(u"\u00b0")
	g.write(')')
	g.write('\n')
	g.write('A1')
	g.write('\t')
	g.write('A2')
	g.write('\t')
	g.write('A3')
	g.write('\t')
	g.write('A4')
	g.write('\t')
	g.write('A1')
	g.write('\t')
	g.write('A2')
	g.write('\t')
	g.write('A3')
	g.write('\t')
	g.write('A4')
	g.write('\n')
	
	###########################################################################
	### For each row in idfile, loop through calculating the angle and outputting to outfile
	for i in range(0,len(idarray[:,0])):
		
		### Initialize variables
		ids=idarray[i,:]
		ids=np.sort(ids)
		Rs=np.zeros([4,3])
		names=[]
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
					names.append(x[12:16].strip())
				if int(x[6:12])==ids[1]:
					Rs[1,0]=float(x[30:38])
					Rs[1,1]=float(x[38:46])
					Rs[1,2]=float(x[46:54])
					names.append(x[12:16].strip())
				if int(x[6:12])==ids[2]:
					Rs[2,0]=float(x[30:38])
					Rs[2,1]=float(x[38:46])
					Rs[2,2]=float(x[46:54])
					names.append(x[12:16].strip())
				if int(x[6:12])==ids[3]:
					Rs[3,0]=float(x[30:38])
					Rs[3,1]=float(x[38:46])
					Rs[3,2]=float(x[46:54])
					names.append(x[12:16].strip())
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
		#g=open(outfile,'w')
		for j in range(0,4):
			g.write(names[j])
			g.write('\t')
		for j in range(0,4):
			g.write(str(int(ids[j])))
			g.write('\t')
		g.write('%.3f' %(phi))
		g.write('\n')
	g.close()
		
### Run
dihedral(args.input, args.inputids, args.output)

