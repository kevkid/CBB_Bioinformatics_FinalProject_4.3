**README for CBB_Bioinformatics_FinalProject_4.3**
---------------------------------------------------------------
Tool that calculates the dihedral angle based on a pdb file and an input-id file containing the corresponing atom ids of four atoms in said file.

# The python tool that accomplishes this task is named dihedrals.py
dihedrals.py takes 2 required inputs (input, inputids) and has one optional (output)
The tool is called from the command line as exemplified below.

Usage:      python3 distcalc.py -i <input file name> -n <input id file name> -o <.txt output filename>

Examples:
```{r NCBI_python, engine="python", highlight=TRUE}
# Usage from terminal:
	     python3 dihedrals.py -i sample-input.pdb -n sample-ids.txt -o sample-output.txt
   		 python3 dihedrals.py -i sample-input.pdb -n sample-ids.txt
```

Input Formats:
```{r NCBI_python, engine="python", highlight=TRUE}
	-i 	string of corresponding pdb file
	-n	string corresponding to text file containing
        integers for atom ids. Each line is a set 
  			of 4 which defines the dihedral angle.
    -o	string containing the name of the file to which the angles are output
```

Format of Input ID file:
  * Text file containing integers in the format of four atom ids per row. The atom ids of each row define each dihedral angle.

| Atom ID1  | Atom ID2 | Atom ID3 | Atom ID4 |
|---|---|---|---|
| 1 | 2 | 3 | 12 |
| 1 | 2 | 5 | 6 |
| 2 | 5 | 6 | 7 |
| 5 | 6 | 7 | 8 |
| ... | ... | ... | ... |
| ... | ... | ... | ... |

Order of Output file formatting:

Atom Names			Atom IDs

| A1 | A2 | A3 | A4 | A1 | A2 | A3 | A4 |	Angle (°) |
|---|---|---|---|---|---|---|---|---|
| N	| CA | C | N | 1 | 2 | 3 | 12 | 99.285 |
| N | CA | CB | CG | 1 | 2 | 5 | 6 | -69.859 |
| CA | CB | CG | CD |	2 | 5 | 6 | 7 | -176.012 |
| CB | CG | CD | NE | 5 | 6 | 7 | 8 | 176.505 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |
