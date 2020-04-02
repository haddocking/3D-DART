# 3D-DART

Standalone version of 3D-DART (originally developed by [Marc van Dijk](https://github.com/marcvdijk)) for Python 2.x series.

## Installation

* Install [numpy](https://numpy.org/) library for Python 2.x:

```
pip install numpy
```

* Clone this repository:

```bash
git clone https://github.com/haddocking/3D-DART.git
```

* Create [X3DNA](http://x3dna.org/) software folder structure:

```bash
cd 3D-DART/software
mkdir X3DNA-linux X3DNA-mac
```

* Download from [http://x3dna.org/](http://x3dna.org/) your specific architecture version and place it `X3DNA-linux` if GNU/Linux or `X3DNA-mac` if macOS.

* Test running the main script, you should get a similar output to this:

```bash
cd /path/to/3D-DART
./RunDART.py
```

```
--> Performing System Checks:
   * Python version is: 2.7.1
   * Could not import Numeric package trying NumPy
   * Importing NumPy package succesfull
--> Your current working directory is: /path/to/3D-DART
--------------------------------------------------------------------------------------------------------------
Welcome to 3D-DART version 1.2   Thu Apr  2 12:03:14 2020
--------------------------------------------------------------------------------------------------------------
--> Parsing command-line arguments:
Wrong command line

```

If you see the message above, then 3D-DART is ready.


## Example

```bash
cd /path/to/3D-DART
./RunDART.py -w NAensemblebuild -f example/struct_*.pdb
```


## Output description

A brief description of the directories:

| Folder  | Description  |
|---|---|
| jobnr1-FileSelector  | Contains the files you used as input. |
| jobnr2-BuildNucleicAcids  | Generates a canonical starting structure of the DNA with the desired sequence. |
| jobnr3-X3DNAAnalyze | Runs a 3DNA analysis over the structure generated in job 2 to obtain initial starting parameters for modeling. |
| jobnr4-ModelNucleicAcids | The actual modeling step. Models are represented as base-pair(step) parameter files. |
| jobnr5-BuildNucleicAcids | Converts the parameter files from job 4 into PDB structure files. |
| jobnr6-X3DNAanalyze | Runs a 3DNA analysis over the modeled structures. |
| jobnr7-NABendAnalyze | Runs a global bend analysis over the modeled structures. |
| jobnr8-PDBeditor | Makes some last-minute changes to the PDB files. This directory contains the final structures. |

In the main directory you will find the log file of your job named 'dart.out'
In case you encounter unexpected errors please look into the log file.

## Reference

If you use 3D-DART please cite:

M. van Dijk and A.M.J.J. Bonvin (2009) "3D-DART: a DNA structure modelling server", Nucl. Acids Res.
37 (Web Server Issue):W235-W239, [doi:10.1093/nar/gkp287](https://doi.org/10.1093/nar/gkp287)

