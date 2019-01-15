#!/usr/bin/env python2.7

USAGE = """
==========================================================================================

Author:				Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center
					for Biomolecular Research, Utrecht university, The Netherlands
Copyright (C):		2007 (DART project)
DART version:		1.2  (25-11-2008)
DART plugin: 		PDBeditor.py
Input:				PDB data file and any of the allowed options, any order and 
					combined.
Output:				A new PDB file or XML representation.
Plugin excecution:	Either command line driven (use -h/--help for the option) or as 
					part of a DART batch sequence.
Plugin function:	A suite of functions to modify PDB files. Features include:
					change nucleic acid nomenclature from a 1-letter to a 3-letter
					code and vice versa; support for the new wwwPDB nucleic-acid abreviation
					scheme; set chain-ID; renumber residues and/or atom
					numbering; check if PDB is valid for HADDOCK (TER statement,
					END statement, CNS nomenclature); place chain ID to location of
					seg ID; split ensemble files or concate PDB files to an ensemble;
					convert PDB to an XML representation.
Examples:			PDBeditor.py -f test.pdb -kn test_fixed.pdb
					PDBeditor.py -f test.pdb -r 1 -c B -adg
Dependencies:		Standard python2.3 or higher. DART package (XMLwriter,
					and Constants modules)

==========================================================================================
"""

"""Import modules"""
import os, sys, re, glob

"""Setting pythonpath variables if run from the command line"""
base, dirs = os.path.split(os.path.dirname(os.path.join(os.getcwd(), __file__)))

if base in sys.path:
    pass
else:
    sys.path.append(base)

from system.XMLwriter import Node
from system.Constants import *


def PluginXML():
    PluginXML = """
<metadata>
 <name>PDB formating options</name>
 <input type="Filetype">.pdb</input>
 <output type="Filetype">.pdb</output>
</metadata>
<parameters>
 <option type="useplugin" form="hidden" text="None">True</option>
 <option type="inputfrom" form="hidden" text="None">1</option>
 <option type="NA1to3" form="checkbox" text="Convert nucleic acid 1 letter to 3 letter notation">False</option>
 <option type="NA3to1" form="checkbox" text="Convert nucleic acid 3 letter to 1 letter notation">False</option> 
 <option type="setchainID" form="text" text="Set PDB chain id"></option>
 <option type="IUPACtoCNS" form="checkbox" text="Convert IUPAC to CNS notation">False</option>
 <option type="reres" form="text" text="Renumber residues starting from"></option>
 <option type="reatom" form="text" text="Renumber atoms starting from"></option>
 <option type="xsegchain" form="checkbox" text="Move chain id to segment id">False</option>
 <option type="noheader" form="checkbox" text="PDB without header lines">False</option>
 <option type="nohetatm" form="checkbox" text="PDB without HETATM records">False</option>
 <option type="nofooter" form="checkbox" text="PDB without footer lines">False</option>
 <option type="pdb2haddock" form="checkbox" text="Make PDB HADDOCK ready">True</option>
 <option type="joinpdb" form="checkbox" text="Join PDB files to one">False</option>
 <option type="splitpdb" form="text" text="Split PDB files based on TER or MODEL statement"></option>
 <option type="name" form="text" text="Give your structure a name"></option>
 <option type="pdb2xml" form="checkbox" text="Convert PDB to DART XML representation">False</option>
</parameters>"""

    return PluginXML


def PluginCore(paramdict, inputlist):
    print "--> Starting PDBeditor"

    """Split ensemble of PDB files in separate PDB files"""
    if not paramdict['splitpdb'] == None:
        pdb = PDBeditor()
        for files in inputlist:
            print("    * Spliting ensemble PDB in individual PDB files on %s statement" % paramdict['splitpdb'])
            pdb.SplitPDB(ensemble=files, mode=paramdict['splitpdb'])
        sys.exit(0)

    filecount = 1  # to keep track of processed files in the PDB joining process

    for files in inputlist:

        pdb = PDBeditor()
        pdb.ReadPDB(files)

        """Perform fixes to the pdb file to make it suitable for HADDOCK"""
        if paramdict['pdb2haddock']:
            paramdict['NA1to3'] = True
            paramdict['NA3to1'] = False
            paramdict['IUPACtoCNS'] = True
            paramdict['reatom'] = int(1)
            paramdict['noheader'] = True
            paramdict['nohetatm'] = True

        """Convert Nucleic-Acids residue one-letter-code to three-letter-code or vice versa"""
        if paramdict['NA1to3']:
            print "    * Convert Nucleic-Acids one-letter-code to three-letter-code"
            pdb.NAresid1to3()
        if paramdict['NA3to1']:
            print "    * Convert Nucleic-Acids three-letter-code to one-letter-code (wwwPDB notation)"
            pdb.NAresid3to1()

        """Convert IUPAC atom naming to CNS"""
        if paramdict['IUPACtoCNS']:
            print "    * Convert IUPAC atom notation to CNS atom notation"
            pdb.IUPACtoCNS()

        """Place HADDOCK chain ID in propper place"""
        if paramdict['xsegchain']:
            print "    * Set seg ID to position of chain ID"
            pdb.XsegChain()

        """Set the chain ID"""
        if paramdict['setchainID'] is not None:
            chainID = paramdict['setchainID'].split(',')
            try:
                old = chainID[0].upper()
                new = chainID[1].upper()
                print "    * Converting chain ID:", old, "to chain ID:", new
                pdb.SetchainID(old=old, new=new)
            except:
                new = chainID[0].upper()
                print "    * Converting all to chain ID:", new
                pdb.SetchainID(new=new)

        """Renumbering residues"""
        if paramdict['reres'] is not None:
            print "    * Renumber residues starting from:", paramdict['reres']
            pdb.Reres(paramdict['reres'])

        """Renumber atoms"""
        if paramdict['reatom'] is not None:
            print "    * Renumber atoms starting from:", paramdict['reatom']
            pdb.Reatom(paramdict['reatom'])
            pdb.CorrectConect(paramdict['reatom'])

        """Make XML representation of pdb"""
        if paramdict['pdb2xml']:
            xml = pdb.PDB2XML()

            root = os.path.basename(files)
            basename, extension = os.path.splitext(root)
            outfile = basename + ".xml"

            print "    * Generating DART XML representation of the PDB as:", outfile

            out = file(outfile, 'w')
            out.write(xml.xml())
            out.close

        """Write new PDB file"""
        if paramdict['pdb2xml'] == False and paramdict['joinpdb'] == False and paramdict['splitpdb'] == None:
            if paramdict['name'] == None:
                root = os.path.basename(files)
                basename, extension = os.path.splitext(root)
                outfile = basename + "_fixed" + extension
            else:
                files = glob.glob('*.pdb')
                basename, extension = os.path.splitext(paramdict['name'])
                if paramdict['name'] in files:
                    count = 1
                    while count < len(files):
                        newname = basename + '-' + str(count) + extension
                        if newname in files:
                            pass
                        else:
                            outfile = newname
                            break
                        count = count + 1
                else:
                    outfile = paramdict['name']

            print "    * Printing fixed pdb file as:", outfile
            pdb.WritePDB(file_out=outfile, join=False, modelnr=0, noheader=paramdict['noheader'],
                         nofooter=paramdict['nofooter'], nohetatm=paramdict['nohetatm'])

        elif paramdict['pdb2xml'] == False and paramdict['joinpdb'] == True and paramdict['splitpdb'] == None:
            if paramdict['name'] == None:
                outfile = 'joined.pdb'
            else:
                outfile = paramdict['name']
            print "    * Append", os.path.basename(files), "to concatenated file:", outfile
            pdb.WritePDB(file_out=outfile, join=True, modelnr=filecount, noheader=True, nofooter=True,
                         nohetatm=paramdict['nohetatm'])

        filecount = filecount + 1


# ================================================================================================================================#
# 										PLUGIN SPECIFIC DEFINITIONS BELOW THIS LINE												 #
# ================================================================================================================================#

class CommandlineOptionParser:
    """Parses command line arguments using optparse"""

    def __init__(self):

        self.option_dict = {}
        self.option_dict = self.CommandlineOptionParser()

    def CommandlineOptionParser(self):

        """Parsing command line arguments"""

        usage = "usage: %prog" + USAGE
        parser = OptionParser(usage)

        parser.add_option("-f", "--file", action="callback", callback=self.varargs, dest="inputfile", type="string",
                          help="Supply pdb inputfile(s). Standard UNIX selection syntax accepted")
        parser.add_option("-a", "--na1to3", action="store_true", dest="NA1to3", default=False,
                          help="Convert nucleic-acid residues from one-letter to three-letter code")
        parser.add_option("-b", "--na3to1", action="store_true", dest="NA3to1", default=False,
                          help="Convert nucleic-acid residues from three-letter to one-letter code")
        parser.add_option("-c", "--setchainid", action="store", dest="setchainID", type="string",
                          help="Convert chain ID. Input as A,B (old,new) or A (all to A)")
        parser.add_option("-i", "--iupactocns", action="store_true", dest="IUPACtoCNS", default=False,
                          help="Convert IUPAC NA atom naming to CNS atom naming")
        parser.add_option("-r", "--reres", action="store", dest="reres", type="int",
                          help="Renumber residues. Options: starting from (number)")
        parser.add_option("-p", "--reatom", action="store", dest="reatom", type="int",
                          help="Renumber atoms. Options: starting from (number)")
        parser.add_option("-s", "--xsegchain", action="store_true", dest="xsegchain", default=False,
                          help="Places chain ID in propper place")
        parser.add_option("-d", "--noheader", action="store_true", dest="noheader", default=False,
                          help="Write pdb file without header lines")
        parser.add_option("-e", "--nohetatm", action="store_true", dest="nohetatm", default=False,
                          help="Write pdb file without hetatm lines")
        parser.add_option("-g", "--nofooter", action="store_true", dest="nofooter", default=False,
                          help="Write pdb file without footer lines (CONECT)")
        parser.add_option("-k", "--pdb2haddock", action="store_true", dest="pdb2haddock", default=False,
                          help="Perform general pdb fixes for HADDOCK (-aipd)")
        parser.add_option("-l", "--joinpdb", action="store_true", dest="joinpdb", default=False,
                          help="Concatenate PDB files")
        parser.add_option("-m", "--splitpdb", action="store", dest="splitpdb", type="string",
                          help="Split ensemble PDB files on MODEL or TER statemend")
        parser.add_option("-n", "--name", action="store", dest="name", type="string", help="name for the new PDB file")
        parser.add_option("-x", "--pdb2xml", action="store_true", dest="pdb2xml", default=False,
                          help="Make DART XML representation of pdb")

        (options, args) = parser.parse_args()

        self.option_dict['input'] = options.inputfile
        self.option_dict['NA1to3'] = options.NA1to3
        self.option_dict['NA3to1'] = options.NA3to1
        self.option_dict['setchainID'] = options.setchainID
        self.option_dict['IUPACtoCNS'] = options.IUPACtoCNS
        self.option_dict['reres'] = options.reres
        self.option_dict['reatom'] = options.reatom
        self.option_dict['xsegchain'] = options.xsegchain
        self.option_dict['noheader'] = options.noheader
        self.option_dict['nohetatm'] = options.nohetatm
        self.option_dict['nofooter'] = options.nofooter
        self.option_dict['pdb2haddock'] = options.pdb2haddock
        self.option_dict['joinpdb'] = options.joinpdb
        self.option_dict['splitpdb'] = options.splitpdb
        self.option_dict['name'] = options.name
        self.option_dict['pdb2xml'] = options.pdb2xml

        if not self.option_dict['input'] == None:
            parser.remove_option('-f')
            arg = self.GetFirstArgument(parser, shorta='-f', longa='--file')
            self.option_dict['input'].append(arg)
            fullpath = self.GetFullPath(self.option_dict['input'])
            self.option_dict['input'] = fullpath

        if parser.has_option('-f'):
            pass
        else:
            parser.add_option("-f", "--file", action="store", dest="dummy2",
                              type="string")  # only needs to be here to complete the argument list, not used!

        return self.option_dict

    def GetFullPath(self, inputfiles):

        currdir = os.getcwd()
        filelist = []

        for files in inputfiles:
            path = os.path.join(currdir, files)
            filelist.append(path)

        return filelist

    def GetFirstArgument(self, parser, shorta, longa):

        """HACK, optparse has difficulties in variable argument lists. The varargs definition solves this but never reports the first
           argument of the list. This definition hacks this issue"""

        parser.add_option(shorta, longa, action="store", dest="temp", type="string",
                          help="Execute custom workflow assembled on the command line. You can execute a single plugin by typing '-p pluginname' or a sequence of plugins by typing '-p plugin1,plugin2...'")

        (options, args) = parser.parse_args()
        first_arg = options.temp
        parser.remove_option(shorta)

        return first_arg

    def varargs(self, option, opt_str, value, parser):

        """Deals with variable list of command line arguments"""

        value = []
        rargs = parser.rargs
        while rargs:
            arg = rargs[0]

            if ((arg[:2] == "--" and len(arg) > 2) or
                    (arg[:1] == "-" and len(arg) > 1 and arg[1] != "-")):
                break
            else:
                value.append(arg)
                del rargs[0]

        setattr(parser.values, option.dest, value)


class PDBeditor:
    def __init__(self, inputfile=None):

        self.title = []
        self.atcounter = 0
        self.allatom_line = []
        self.header = []
        self.footer = []
        self.end = []
        self.model = []
        self.label = []
        self.atnum = []
        self.elem = []
        self.atname = []
        self.atalt = []
        self.resname = []
        self.chain = []
        self.resnum = []
        self.resext = []
        self.coord = []
        self.occ = []
        self.b = []
        self.hdoc_chain = []
        self.sequence = {}
        self.firstatnr = 1

    def ReadPDB(self, inputfile, debug=0):

        # check if passed filename string or a file descriptor
        if type(inputfile) == type(sys.stdin):
            readfile = inputfile
        else:
            readfile = file(inputfile, 'r')

        lines = readfile.readlines()
        self.ReadPDBlines(lines, debug)

    def ReadPDBlines(self, lines, debug=0):

        """
        Reads a list of PDB-format file lines in to the Protein class object.
        Thus can be called by another routine that already has the lines in a list.
        Returns the number of atoms read in.
        """

        i = 0
        atom_hetatm = re.compile('(ATOM  |TER   |HETATM)')
        head = re.compile('^(HEADER|COMPND|SOURCE|JRNL|HELIX|REMARK|SEQRES|CRYST1|SCALE|ORIG)')
        title = re.compile('^TITLE')
        foot = re.compile('(CONECT|MASTER)')
        end = re.compile('(END)')
        model = re.compile('(MODEL)')
        element = re.compile('[A-Za-z ][A-Za-z]')
        for line in lines:
            if atom_hetatm.match(line):
                line = line[:-1]

                # Add TER statement if change in chainid
                if len(self.chain) and not self.chain[-1] == line[21]:
                    self.label.append('TER   ')
                    for b in (
                    self.allatom_line, self.atnum, self.atname, self.resnum, self.resname, self.chain, self.atalt,
                    self.resext, self.occ, self.b, self.hdoc_chain, self.elem):
                        b.append(blank)
                    self.coord.append((0.000, 0.000, 0.000))
                if not line.startswith("TER"):
                    self.allatom_line.append(line)
                    self.label.append(line[0:6])  # atom label
                    self.atnum.append(int(line[6:12]))  # atom number
                    self.atname.append(line[12:16])  # atom type
                    self.atalt.append(line[16:17])
                    self.resname.append(line[17:21])  # residu name
                    self.chain.append(line[21])  # chain
                    self.resnum.append(int(line[22:26]))  # residu number

                try:
                    self.resext.append(line[27])
                except:
                    self.resext.append(blank)

                try:
                    self.coord.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))  # X,Y,Z coordinates
                except:
                    if not line[0:3] == 'TER' or line[0:5] == 'MODEL':
                        print "    * ERROR: coordinate error in line:"
                        print "     ", line

                try:
                    self.occ.append(float(line[54:60]))
                except:
                    self.occ = ([1.00] * len(self.atnum))

                try:
                    self.b.append(float(line[60:66]))  # B factor
                except:
                    self.b = ([0.00] * len(self.atnum))

                try:
                    self.hdoc_chain.append(line[72])  # SEGID
                except:
                    self.hdoc_chain = ([blank] * len(self.atnum))

                if element.match(line[76:78]):  # Get first element in elementlist
                    self.elem.append(line[76:78])
                else:
                    self.elem.append(line[12:14])

                self.atcounter += 1
                i += 1

            elif head.match(line):
                self.header.append(line[:-1])
            elif foot.match(line):
                self.footer.append(line[:-1])
            elif end.match(line):
                self.end.append(line[:-1])
            elif model.match(line[:-1]):
                self.model.append(line)
            elif title.match(line):
                self.title.append(line[:-1])

        self.firstatnr = self.atnum[
            0]  # Need to know original number of first atom for possible CONECT statement correction when renumbering atoms

        if debug:
            return len(self.atnum), self.atcounter

    def WritePDB(self, file_out, join=False, modelnr=0, noheader=False, nofooter=False, nohetatm=False):

        """
        Saves the Protein class object to a PDB-format file
        if noheader = True, no header (REMARK etc.) or footer lines are written
        if nohetatm = True, no hetero atoms are written
        """

        if join == True:
            out = open(file_out, 'a')
        else:
            out = file(file_out, 'w')

        if noheader == False:
            for i in range(len(self.title)):
                out.write('%s\n' % self.title[i])
            for i in range(len(self.header)):
                out.write("%s\n" % self.header[i])

        if join == True:
            out.write('MODEL ' + str(modelnr) + '\n')

        for i in xrange(len(self.resnum)):
            if self.label[i] == 'ATOM  ':
                self.WritePDBline(out, i)
            elif self.label[i] == 'TER   ':
                out.write('TER   \n')
            elif self.label[i] == 'HETATM' and not nohetatm:
                out.write("%s\n" % self.allatom_line[i])

        if nofooter == False:
            for i in range(len(self.footer)):
                out.write("%s\n" % self.footer[i])

        if join == False:
            if len(self.end) == 3:
                for i in range(len(self.end)):
                    out.write("%s\n" % self.end[i])
            else:
                self.end = ['END']
                for i in range(len(self.end)):
                    out.write("%s\n" % self.end[i])
        else:
            self.end = ['ENDMDL']
            for i in range(len(self.end)):
                out.write("%s\n" % self.end[i])

        out.close()

    def WritePDBline(self, FD, i):

        """
        Writes a single line of data in the PDB-format
            called by writePDB
        """

        FD.write('%-6s%5i %-4s%1s%-4s%1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%10s%2s\n' %
                 (self.label[i], self.atnum[i], self.atname[i], self.atalt[i],
                  self.resname[i], self.chain[i], self.resnum[i], self.resext[i],
                  self.coord[i][0], self.coord[i][1], self.coord[i][2], self.occ[i], self.b[i], blank, self.elem[i]))

    def SplitPDB(self, ensemble=None, mode=None):

        """
        Split ensemble PDB files in seperate PDB files based on MODEL or TER tag
        """

        # check if passed filename string or a file descriptor
        if type(ensemble) == type(sys.stdin):
            readfile = ensemble
        else:
            readfile = file(ensemble, 'r')

        lines = readfile.readlines()

        mode = mode.upper()
        atom_hetatm = re.compile('(ATOM  |HETATM)')
        model = re.compile('(' + mode + ')')
        modelcount = 1
        models = {}
        models[modelcount] = []
        linecount = len(lines)
        linenr = 0
        while linenr < linecount:
            line = lines[linenr].strip()
            if model.match(line):
                if not len(models[modelcount]) == 0:
                    modelcount += 1
                    if models.has_key(modelcount) == False:
                        models[modelcount] = []
                linenr += 1
            if atom_hetatm.match(line):
                models[modelcount].append(linenr)
            linenr = linenr + 1

        if len(models) == 1:
            print "    * No splitting occured, splitting statement not found"
        else:
            for model in models.keys():
                outfile = os.path.splitext(ensemble)[0] + '_' + str(model) + '.pdb'
                out = file(outfile, 'w')
                print "    * Writing model %s as %s" % (model, outfile)
                for line in models[model]:
                    out.write(lines[line])
                out.write('END')
                out.close()

    def NAresid1to3(self):

        """
        Convert list of 1-letter nucleic-acid code sequence to 3-letter code and update resname
        """

        seq3 = []

        for resid1 in self.resname:
            try:
                resid3 = NAres3[
                    NAres1.index(resid1.upper())]  # If NAresid is one-letter code, convert to three-letter code
                seq3.append(resid3)
            except ValueError, err:
                if resid1.upper() in AAres3:  # If resid is amino-acid three letter code, just append
                    seq3.append(resid1.upper())  # Amino-acid one letter code in PDB not accepted(expected)
                elif resid1.upper() == 'HOH ':  # Waters are neglected, just append.
                    seq3.append(resid1.upper())
                elif resid1.upper() in NAres3:  # If NAresid allready in three letter code, just append
                    seq3.append(resid1.upper())
                else:
                    print "      - WARNING: no match for residue: %s" % (
                    resid1)  # If not of the above, raise exception.
                    seq3.append(resid1.upper())

        if len(seq3) == len(self.resname):
            self.resname = seq3
        else:
            pass

    def NAresid3to1(self):

        """
        Convert list of 3-letter nucleic-acid code sequence to 1-letter code and update resname. The 1-letter code is the new
        (2006) wwwPDB notation. This is DA,DT,DC,DG for DNA and RA,RU,RG,RC for RNA.
        """

        print "      - WARNING: The conversion of nucleic-acid three-letter code to two-letter code does not check for ribose or"
        print "                 deoxy-ribose. If Uracil is found the structure is regarded as RNA otherwise as DNA. Please check"
        print "                 your structure in case of mixed conformations."

        seq1 = []
        THREELETTER = ['--- ', 'CYT ', 'THY ', 'GUA ', 'ADE ', 'URI ']
        DNA1LETTER = ['  - ', ' DC ', ' DT ', ' DG ', ' DA ', ' RU ']
        RNA1LETTER = ['  - ', ' RC ', ' DT ', ' RG ', ' RA ', ' RU ']

        if 'URI ' in self.resname:
            RNA = True
            DNA = False
        else:
            DNA = True
            RNA = False

        for resid3 in self.resname:
            try:
                if RNA == True:
                    resid1 = RNA1LETTER[THREELETTER.index(resid3.upper())]
                    seq1.append(resid1)
                elif DNA == True:
                    resid1 = DNA1LETTER[THREELETTER.index(resid3.upper())]
                    seq1.append(resid1)
            except ValueError, err:
                print "      - WARNING: no match for residue:", resid3
                seq1.append(resid3.upper())

        if len(seq1) == len(self.resname):
            self.resname = seq1
        else:
            pass

    def SetchainID(self, old=None, new=None):

        """
        Convert the chain ID from the old ID to the user supplied new ID if None do nothing.
        Option examples: (A) all to A, (A,B) all A to B. Lower case is converted to upper case.
        """

        newchainseq = []

        if not old and new == None:
            for chainid in self.chain:
                if chainid == old:
                    newchainseq.append(new)
                else:
                    newchainseq.append(chainid)

        else:
            for chainid in self.chain:
                newchainseq.append(new)

        if len(newchainseq) == len(self.chain):  # in case of any errors that yield non-equal arrays
            self.chain = newchainseq
        else:
            pass

    def IUPACtoCNS(self):

        """
        Convert IUPAC atom type notation to CNS atom type notation. Get info from IUPAC and CNS lists from Constants.py.
        Currently only conversion of nucleic-acid atom types.
        """

        newatomseq = []

        for atom in self.atname:
            try:
                newatom = CNS[IUPAC.index(atom)]
                newatomseq.append(newatom)
            except ValueError:
                newatomseq.append(atom)

        if len(newatomseq) == len(self.atname):
            self.atname = newatomseq
        else:
            pass

    def PDB2XML(self):

        """
        Makes a XML representation of the PDB. Needs system.XMLwriter
        """

        main = Node("DART_pdbx")

        acount = 0
        lastchain = ' '
        lastresnum = ' '

        for i in xrange(len(self.atnum)):
            if i == 0:
                lastchain = self.chain[i]
                lastresnum = self.resnum[i]

                chain = Node("chain", ID=lastchain)
                resid = Node("resid", ID=self.resname[i], nr=str(lastresnum))
                atom = Node("atom", ID=self.atname[i], nr=str(self.atnum[i]), corx=str(self.coord[i][0]),
                            cory=str(self.coord[i][1]), corz=str(self.coord[i][2]), occ=str(self.occ[i]),
                            b=str(self.b[i]))

                resid += atom
                chain += resid
                main += chain
            else:
                if self.chain[i] == lastchain:
                    lastchain = self.chain[i]
                    if self.resnum[i] == lastresnum:
                        atom = Node("atom", ID=self.atname[i], nr=str(self.atnum[i]), corx=str(self.coord[i][0]),
                                    cory=str(self.coord[i][1]), corz=str(self.coord[i][2]), occ=str(self.occ[i]),
                                    b=str(self.b[i]))
                        lastresnum = self.resnum[i]
                        resid += atom
                    else:
                        lastresnum = self.resnum[i]
                        resid = Node("resid", ID=self.resname[i], nr=str(lastresnum))
                        atom = Node("atom", ID=self.atname[i], nr=str(self.atnum[i]), corx=str(self.coord[i][0]),
                                    cory=str(self.coord[i][1]), corz=str(self.coord[i][2]), occ=str(self.occ[i]),
                                    b=str(self.b[i]))
                        resid += atom
                        chain += resid
                else:
                    lastchain = self.chain[i]
                    lastresnum = self.resnum[i]
                    chain = Node("chain", ID=lastchain)
                    resid = Node("resid", ID=self.resname[i], nr=str(lastresnum))
                    atom = Node("atom", ID=self.atname[i], nr=str(self.atnum[i]), corx=str(self.coord[i][0]),
                                cory=str(self.coord[i][1]), corz=str(self.coord[i][2]), occ=str(self.occ[i]),
                                b=str(self.b[i]))
                    resid += atom
                    chain += resid
                    main += chain

        return main

    def Reres(self, start):

        """
        Renumber residues. Option example: (4) renumber starting from 4.
        """
        start = int(start)
        lastresnum = -9999
        lastchain = ' '
        lastresname = ' '
        lastext = ' '
        idres = start - 1
        icount = 0

        for i in xrange(len(self.resnum)):
            if i == 0:
                icount += 1
                lastchain = self.chain[i]
                lastresname = self.resname[i]
                lastresnum = self.resnum[i]
                lastext = self.resext[i]
                self.resnum[i] = start
                idres += 1
            else:
                if (self.chain[i] != lastchain or lastresnum != self.resnum[i] or
                            lastresname != self.resname[i] or lastext != self.resext[i]):
                    icount += 1
                    idres += 1
                    lastchain = self.chain[i]
                    lastresname = self.resname[i]
                    lastresnum = self.resnum[i]
                    self.resnum[i] = idres
                else:
                    self.resnum[i] = idres

    def Reatom(self, start):

        """
        Renumber atoms. Option example: (4) renumber complete list starting from 4.
        """

        start = int(start)

        newatomnum = range(start, (len(self.atnum) + start))

        if len(newatomnum) == len(self.atnum):
            self.atnum = newatomnum
        else:
            pass

    def CorrectConect(self, number):

        """
        Correct the CONECT statement when renumbering atoms.
        """

        diff = number - self.firstatnr
        correctconect = []

        for line in self.footer:
            p = line.split()
            if p[0] == 'CONECT':
                for atomnr in xrange(1, len(p)):
                    p[atomnr] = int(p[atomnr]) + diff
                correct = "CONECT"
                for n in p[1:]:
                    correct = correct + ("%5i" % n)
                correctconect.append(correct)
            else:
                correctconect.append(line)

        self.footer = correctconect

    def XsegChain(self):

        """
        Copy SEGID to CHAIN location.
        """

        hdoc_chain = []
        for chainid in self.hdoc_chain:
            chainid.strip()
            if chainid == ' ':
                pass
            else:
                hdoc_chain.append(chainid)

        if len(hdoc_chain) > 0:
            self.chain = self.hdoc_chain
        else:
            pass


if __name__ == '__main__':

    """Running from the command line"""
    from optparse import *

    """Parse command line arguments"""
    option_dict = CommandlineOptionParser().option_dict

    """Check for input"""
    if option_dict['input'] == None:
        print "    * Please supply pdb file using option -f or use option -h/--help for usage"
        sys.exit(0)
    else:
        inputlist = option_dict['input']

    """Envoce main functions"""
    PluginCore(option_dict, inputlist)
    sys.exit(0)
