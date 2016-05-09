**README for Calculating Dihedral Angles from a PDB-file**
---------------------------------------------------------------
R tool that calculates the dihedral angle based on a pdb file and an input-id file containing the corresponing atom ids of four atoms in said file. A tool that accomplishes this task in the language python can be found [here] (https://github.com/peter-mm-williams/Dihedral_Angle_Calc.git). This tool is part of a set of bioinformatic and biological structure tools created for CBB752 at Yale University in the Spring 2016. The website containing the set of tools can be found [here] (http://cbb752spring2016.github.io/Structure).


# The python tool that accomplishes this task is named dihedrals.py
# The R tool that accomplishes this task is named 4_3.r
## General

The tool is named 4_3.r. It takes 2 required inputs (input, inputids) and has one optional (output)
The tool is called from the command line as exemplified below. It takes in a pdb file and a text file containing sets of 4 atom ids which define the dihedral angle in the pdb file. Sample files (sample-input.pdb and sample-ids.txt) are included in the repository. The corresponding output is given in sample-output.txt. 


## Usage

### Usage:      Rscript 4_3.r sample-ids.txt sample-input.pdb output.txt

### Examples:
```{r NCBI_python, engine="python", highlight=TRUE}
# Usage from terminal:
Rscript 4_3.r sample-ids.txt sample-input.pdb output.txt
```

## Input and Output formats
### Input Formats:
####R:
 * Argument 1: The Input ID file
 * Argument 2: The PDB file
 * Argument 3: The output file
  
### Format of Input ID file:
  * Text file containing integers in the format of four atom ids per row. The atom ids of each row define each dihedral angle.

| Atom ID1  | Atom ID2 | Atom ID3 | Atom ID4 |
|---|---|---|---|
| 1 | 2 | 3 | 12 |
| 1 | 2 | 5 | 6 |
| 2 | 5 | 6 | 7 |
| 5 | 6 | 7 | 8 |
| ... | ... | ... | ... |
| ... | ... | ... | ... |

### Order of Output file formatting:

| --------Atom  Names------- | ----------Atom IDs---------- | Angle | 
|---|---|---|

| A1 | A2 | A3 | A4 | A1 | A2 | A3 | A4 |	Angle (°) |
|---|---|---|---|---|---|---|---|---|
| N	| CA | C | N | 1 | 2 | 3 | 12 | 99.285 |
| N | CA | CB | CG | 1 | 2 | 5 | 6 | -69.859 |
| CA | CB | CG | CD |	2 | 5 | 6 | 7 | -176.012 |
| CB | CG | CD | NE | 5 | 6 | 7 | 8 | 176.505 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |
