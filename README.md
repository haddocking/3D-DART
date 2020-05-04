# 3D-DART

Standalone version of 3D-DART (originally developed by [Marc van Dijk](https://github.com/marcvdijk)) for Python 2.x series.

## Installation

* We recommend to use a virtual environment to install 3D-DART:

```
virtualenv --python=python2.7 dart
cd dart
source bin/activate
```

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

* Download from [http://x3dna.org/](http://x3dna.org/) your specific architecture version and place it `X3DNA-linux` if GNU/Linux or `X3DNA-mac` if macOS. For example, `X3DNA-mac` folder should contain:

```
[3D-DART/software/X3DNA-mac]$ ls
BASEPARS FIBER    bin
```

and more in deep:

```
BASEPARS/:
ATOMIC           Atomic_G.pdb     Atomic_U.pdb     Block_BP.alc     Pxyz.dat         RNA_BASES        col_mname.dat    misc_3dna.par    ps_image.par     trans_pep.alc
Atomic_A.pdb     Atomic_P.pdb     Atomic_i.pdb     Block_R.alc      PxyzH.dat        baselist.dat     fig_image.par    my_header.r3d    raster3d.par     trans_pep.pdb
Atomic_C.pdb     Atomic_T.pdb     BLOCK            Block_Y.alc      README           col_chain.dat    help3dna.dat     ndb_raster3d.par rotz90

FIBER/:
Data   Str02  Str05  Str08  Str11  Str14  Str17  Str20  Str23  Str26  Str29  Str32  Str35  Str38  Str41  Str44  Str47  Str50  Str53
README Str03  Str06  Str09  Str12  Str15  Str18  Str21  Str24  Str27  Str30  Str33  Str36  Str39  Str42  Str45  Str48  Str51  Str54
Str01  Str04  Str07  Str10  Str13  Str16  Str19  Str22  Str25  Str28  Str31  Str34  Str37  Str40  Str43  Str46  Str49  Str52  Str55

bin/:
EnergyPDNA.exe anyhelix       cehs           dcmnfile       fiber          get_part       nmr_ensemble   pdb2img        regular_dna    std_base
alc2img        block_atom     comb_str       del_ms         find_pair      manalyze       nmr_strs       r3d_atom       rotate_mol     step_hel
analyze        blocview       cp_std         ex_str         frame_mol      mstack2img     o1p_o2p        rebuild        stack2img
```

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

