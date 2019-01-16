#!/usr/bin/env python2.7

USAGE= """
==========================================================================================

Author:				Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center for 
					Biomolecular Research, Utrecht university, The Netherlands.
Copyright (C):		2006 (DART project)
DART version:		1.2 (25-11-2008)
DART plugin: 		X3DNAanalyze.py
Input:				PDB files, 3DNA .out files or files containing a list of 
					these files. Any of the allowed options, any order or combined.
Output:				.par,.out,.dat,.aux,.cf7,.helical,.curves for each analysed 
   					structure. A multianalysis.stat file and a selection.list file. 			
Plugin excecution:	Either command line driven (use -h/--help for the option) or as 
					part of a DART batch sequence.
Plugin function:	This plugin handels all 3DNA analysis routines. It requiers 3DNA 
					to be present in the main /DART/third-party/ directory. The plugin 
					takes PDB files as input.PDB files may contain only nucleic acids 
					or mixed systems. The plugin generates either the input for the 
					'analyze' routine of 3DNA or runs the complete analysis route 
					generating various output files. To compare structures a 
					multistructure analysis routine is implemented. This routine gathers 
					various parameters from the output file and base-pair/baise-pair 
					step parameter files of all individual structures and calculateds 
					statistical meaningfull data. This multistructure analysis routine 
					can also be run without running the 3DNA analysis routines by 
					either supplying the program with 3DNA .par files or .out files or 
					a list file containing the names of the input files. 
Examples:			X3DNAanalyze -f *.pdb
					X3DNAanalyze -f selection.list			
Dependencies:		Standard python2.3 or higher modules, 3DNA

==========================================================================================
"""

"""Import modules"""
import sys, os, glob, re, copy, string
from operator import itemgetter
from time import ctime
from numpy import *

"""Setting pythonpath variables if run from the command line"""
base, dirs = os.path.split(os.path.dirname(os.path.join(os.getcwd(), __file__)))

if base in sys.path:
	pass
else:
	sys.path.append(base)

"""Import DART specific modules"""
from system.Utils import FileRootRename, TransformDash
from system.IOlib import InputOutputControl
from system.Constants import *
from QueryPDB import GetSequence, NAsummery
from PDBeditor import PDBeditor

def PluginXML():
	PluginXML = """ 
<metadata>
 <name>3DNA nucleic acid analysis routines</name>
 <function>This plugin handels all 3DNA analysis routines. It requiers 3DNA to be 
  present in the main /DART/third-party/ directory. The plugin takes PDB files as 
  input.PDB files may contain only nucleic acids or mixed systems. The plugin 
  generates either the input for the 'analyze' routine of 3DNA or runs the complete 
  analysis route generating various output files. To compare structures a 
  multistructure analysis routine is implemented. This routine	gathers various 
  parameters from the output file and base-pai/baise-pair step parameter files of 
  all individual structures and calculateds statistical meaningfull data. This 
  multistructure analysis routine can also be run without running the 3DNA analysis 
  routines by either supplying the program with 3DNA .out files or a 
  list file containing the names of the input files.</function>
 <input type="filetype">.pdb,.out</input>
 <output type="filetype">.par,.out,.dat,.aux,.cf7,.helical,.list,.curves,multiout.stat</output>
</metadata>
<parameters>
 <option type="useplugin" form="hidden" text="None">True</option>
 <option type="inputfrom" form="hidden" text="None">1</option>
 <option type="onlyinput" form="checkbox" text="Run only find_pair command, nothing parsed to analyze">False</option>
 <option type="singlehelix" form="checkbox" text="treat the whole structure as a continuous single helix. Useful for get all backbone torsion angles">False</option>
 <option type="curvesinput" form="checkbox" text="get Curves input for a duplex">False</option>
 <option type="helregion" form="checkbox" text="generate a separate output file for each helical region">False</option>
 <option type="allbasepairs" form="checkbox" text="find all base-pairs and higher base associations">False</option>
 <option type="deformener" form="checkbox" text="Select structures with the least number of unpairjng/mispairing events for multistructure analysis">False</option>
 <option type="multistructure" form="checkbox" text="Perform multistructure analysis">True</option>
 <option type="master" form="file" text="Provide a master file for multistructure analysis"></option>
</parameters>"""
	
	return PluginXML
	
def PluginCore(paramdict, inputlist):
	
	"""Checking inputlist"""
	checked = InputOutputControl()
	checked.CheckInput(inputlist)
	
	if checked.checkedinput.has_key('.pdb'):
		x3dna = X3DNAanalyze(paramdict)
		x3dna.Run3DNA(checked.checkedinput['.pdb'])
		checked.InputUpdate(".pdb",".out")
		x3dna.RunEnerCalc(checked.checkedinput['.out'])
		
	"""Running Multistructure analysis"""
	if paramdict['multistructure'] == True:
		print "--> Performing multistructure analysis"
		if len(checked.checkedinput['.out']) == 1:
			print "    * WARNING: only 1 out file in input. Not performing multi-structure analysis"
		elif len(checked.checkedinput['.out']) > 1:
			print "    * Performing multistructure analysis on", len(checked.checkedinput['.out']), "parameter files" 
			multiout = MultiStructureAnalysis(checked.checkedinput['.out'])
			multiout.ReadOutfiles()
			multiout.ReadEnerfiles()
			if not paramdict['master'] == None:
				print "    * Using", paramdict['master'], "as sequence master file"
				multiout.GetMasterSeq(paramdict['master'])
			else:
				multiout.GetBaseseq()
				multiout.GetMasterSeq()
			print "    * Writing nucleic-acid pairing information to the file 'napairing.stat'"
			multiout.PairStats()
			if paramdict['deformener'] == True:
				print "    * Set sequence with the least number of unpairing/mispairing events as sequence master."	
				multiout.AutoMaster()
			print "    * Writing statistics to the file 'multiout.stat'"
			multiout.WriteStats()
			print "    * Writing the file 'selection.list' with the structures that match the master file"
			multiout.FileList()
		else:		
			print "    * WARNING: no par files in input, passing"
			
	"""Rounding up"""
	print "--> Finished X3DNAanalyze Core jobs"

#================================================================================================================================#
# 									PLUGIN SPECIFIC DEFINITIONS BELOW THIS LINE													 #
#================================================================================================================================#

class CommandlineOptionParser:
	
	"""Parses command line arguments using optparse"""
	
	def __init__(self):
		
		self.option_dict = {}
		self.option_dict = self.CommandlineOptionParser()
	
	def CommandlineOptionParser(self):
	
		"""Parsing command line arguments"""
	
		usage = "usage: %prog" + USAGE
		parser = OptionParser(usage)

		parser.add_option( "-f", "--file", action="callback", callback=self.varargs, dest="inputfile", type="string", help="Supply pdb inputfile(s)")
		parser.add_option( "-o", "--onlyinput", action="store_true", dest="onlyinput", default=False, help="Run only find_pair command, nothing parsed to analyze")
		parser.add_option( "-s", "--singlehelix", action="store_true", dest="singlehelix", default=False, help="treat the whole structure as a continuous single helix. Useful for get all backbone torsion angles (.outs)")
		parser.add_option( "-c", "--curvesinput", action="store_true", dest="curvesinput", default=False, help="get Curves input for a duplex")
		parser.add_option( "-r", "--helregion", action="store_true", dest="helregion", default=False, help="generate a separate output file for each helical region")
		parser.add_option( "-a", "--allbasepairs", action="store_true", dest="allbasepairs", default=False, help="find all base-pairs and higher base associations")
		parser.add_option( "-e", "--deformener", action="store_true", dest="deformener", default=False, help="Selecting structures based on base-pair and base-pair step deformation energy")
		parser.add_option( "-m", "--multistructure", action="store_true", dest="multistructure", default=False, help="Perform multistructure analysis")
		parser.add_option( "-p", "--master", action="store", dest="master", help="Use sequence of given .out file as master sequence for multi-structure analysis")
		
		(options, args) = parser.parse_args()
		
		self.option_dict['input'] = options.inputfile
		self.option_dict['onlyinput'] = options.onlyinput
		self.option_dict['singlehelix'] = options.singlehelix	
		self.option_dict['curvesinput'] = options.curvesinput
		self.option_dict['helregion'] = options.helregion
		self.option_dict['deformener'] = options.deformener
		self.option_dict['allbasepairs'] = options.allbasepairs
		self.option_dict['multistructure'] = options.multistructure
		self.option_dict['master'] = options.master
			
		if not self.option_dict['input'] == None:
			parser.remove_option('-f')
			arg = self.GetFirstArgument(parser, shorta='-f', longa='--file')
			self.option_dict['input'].append(arg)
			fullpath = self.GetFullPath(self.option_dict['input'])
			self.option_dict['input'] = fullpath

		if parser.has_option('-f'):
			pass
		else:
			parser.add_option( "-f", "--file", action="store", dest="dummy2", type="string") #only needs to be here to complete the argument list, not used!
	
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

		parser.add_option( shorta, longa, action="store", dest="temp", type="string", help="Execute custom workflow assembled on the command line. You can execute a single plugin by typing '-p pluginname' or a sequence of plugins by typing '-p plugin1,plugin2...'")

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

class X3DNAanalyze:
	
	"""
	Wrapper around the 3DNA analysis program. Program needs to be present in the ../DART/third-party/ directory.
	""" 
	
	def __init__(self, paramdict=None):
		
		self.paramdict = paramdict
	
	def _ConstructOptionString(self):	

		"""Constructing command line option string"""
	        
		option = []
	        option.append('-tz')
	        print "--> Running dna analysis with default option set:"
	        print "    *(-t) read HETATM records"
	        print "    *(-z) more detailed base-pairing information in the output"
	        if self.paramdict['singlehelix'] ==  True:
	        	option.append('s')
	        	print "    *(-s) treat the whole structure as a continuous single helix. Useful for get all backbone torsion angles"
	        if self.paramdict['helregion'] ==  True:
	        	option.append('d')
	        	print "    *(-d) generate a separate output file for each helical region"
	        if self.paramdict['allbasepairs'] == True:
	        	option.append('p')
	        	print "    *(-p) find all base-pairs and higher base associations"
	        if self.paramdict['curvesinput'] == True:
	        	option.append('c')
	        	print "    *(-c) get Curves input for a duplex"
	        
	        self.optionstring = ''.join(option)

	def _CleanUp(self):
	
		"""Cleaning up"""
		
		print "--> Cleaning up temporary files"
		trash = ['bestpairs.pdb','bp_order.dat','hel_regions.pdb','hstacking.pdb','poc_haxis.r3d','stacking.pdb','col_chains.scr','bp_step.alc',
	        	 'col_helices.scr','ref_frames.dat','tmp_file','multiplets.pdb','mulbp.inp','mref_frames.dat','allpairs.pdb','bp_step.pdb']	
		for n in trash:
			if os.path.isfile(n):
				os.remove(n)	
			else:
				pass

	def Run3DNA(self, inputlist):

		"""Running X3DNA analysis command"""
		
		self._ConstructOptionString()
		
		for files in inputlist:
			if self.paramdict['onlyinput'] == True:
				print "--> Only running the 3DNA find_pair command thus only generating input file for 3DNA analysis routine for the file:", files
				basename,extension = os.path.splitext(files)
				outfile = basename+".inp"
				files2 = os.path.split(files)[-1]
				if files != files2:
				   os.system("ln -s %s" % files)
				   files = files2
				cmd = "find_pair "+(self.optionstring)+" "+(files)+" "+(outfile)
				output = commands.getoutput(cmd)
				print output
			elif self.paramdict['curvesinput'] == True:
				print "--> Getting Curves input for the file", files
				basename,extension = os.path.splitext(files)
				outfile = basename+".curves"
				files2 = os.path.split(files)[-1]
				if files != files2:
				   os.system("ln -s %s" % files)
				   files = files2
				cmd = "find_pair "+(self.optionstring)+" "+(files)+" "+(outfile)
				os.system(cmd)
			else:
				print "--> Running both the 3DNA find_pair and analysis routine for the file:", files
				files2 = os.path.split(files)[-1]
				if files != files2:
				   os.system("ln -s %s" % files)
				   files = files2
				cmd = "find_pair "+(self.optionstring)+" "+(files)+" stdout | analyze"
				os.system(cmd)
			x3dna_files = ['bp_step.par','auxiliary.par','cf_7methods.par','ref_frames.dat','bp_helical.par','bp_step.pdb','bp_step.alc']			 	 
			extensions = ['.par','.aux','.cf7','.dat','.helical','_dna.pdb','.alc']
			for n in x3dna_files:
				if os.path.isfile(n):
					basename,extension = os.path.splitext(files)	 
					ext = extensions[x3dna_files.index(n)]
					FileRootRename(n,ext,basename)	   
				else:
					pass 
		self._CleanUp()
	
	def RunEnerCalc(self,inputlist):
	
		for files in inputlist:
			print "--> Calculating base-pair and base-pair step deformation energy for the file:", files
			outfile = os.path.splitext(files)[0]+".ener"
			cmd1 = "EnergyPDNA.exe -s "+files+" >"+outfile
			cmd2 = "EnergyPDNA.exe -b "+files+" >>"+outfile
			os.system(cmd1)
			os.system(cmd2)

class MultiStructureAnalysis:

	"""
	Performs a multistructure analysis by pulling all data from the individual runs together
	and calculating various statistical parameters. The class destinguishes between structures
	with a similar sequence and unique structures
	"""
	
	def __init__(self, outfiles=None):
		
		self.outfiles = outfiles
		self.origin = {}		
		self.pairs = {}
		self.bp = {}
		self.bpstep = {}
		self.groove = {}
		self.helical = {}
		self.lamdangl = {}
		self.globalp = {}
		self.angle11 = {}
		self.angle12 = {}
		self.angle21 = {}
		self.angle22 = {}
		self.sstvirtb = {}
		self.helixrad = {}
		self.position= {}
		self.bpstepener = {}
		self.bpener = {}
		
		self.basenr = 0
		self.basemoltype = {}
		self.basechainlib = {}
		self.basesequence = {}
		self.selected = []
		self.rejected = []
		
		self.average = {}
		self.zero = 0.0001
	
	def _ReadTable(self, linenr=None, outfile=None):	
		
		"""Reading all table lines of seleted table. Readstart = linenr stop reading at first 
		   instance of *, newline or ~. Store data in tablelines and return"""
		
		tablelines = []
		
		count = linenr+1
		while count < len(outfile):
			if outfile[count][0] == "*":
				break
			elif outfile[count][0] == "\n":
				break	
			elif	outfile[count][10] == "~":
				break
			else:	
				tablelines.append(outfile[count].split())
			count = count+1
	
		return tablelines
	
	def _SortTable(self,tablein=None,tableout=None,infile=None):
		
		"""Sorting tablelines from Readtable and extract colums. Temporary store in columns that 
		   append to library (self) with filename as key"""
		
		columns = []
		
		for n in range(len(tablein[0])):
			columns.append([])
			
		for lines in tablein:
			for n in range(len(lines)):
				columns[n].append(lines[n])
		
		tableout[infile] = columns				
	
	def _MatchMaster(self):
		
		"""Match all read-in outfiles to the (supplied) master file"""
		
		for structures in self.pairs:
			if self.pairs[structures][1] == self.basesequence[self.chainid][0]:
				self.selected.append(structures)
			else:
				self.rejected.append(structures)
	
	def _MatchPairs(self, linenr=None,outfile=None,infile=None):	
		
		tablelines = []
		
		count = linenr+1
		while count < len(outfile):
			if outfile[count][0] == "*":
				break
			elif outfile[count][0] == "\n":
				break	
			elif	outfile[count][10] == "~":
				break
			else:	
				tablelines.append(outfile[count].split())
			count = count+1
	
		intnr = []
		residnr1 = []
		residnr2 = []
		resid1 = []
		resid2 = []
		
		splitter1 = re.compile('[.]')
		splitter2 = re.compile('_')
		splitter3 = re.compile('\]|\[')
		for lines in tablelines:
			intnr.append(lines[0])
			tmp1 = splitter1.split(lines[2])
			tmp2 = []
			for n in tmp1:	
				if n == '':
					pass
				else:
					tmp2.append(n)	
			residnr1.append(int(splitter2.split(tmp2[1])[0]))
			try:													#split is different for one and three letter code
				residnr2.append(int(splitter2.split(tmp2[2])[0]))	#one letter code
			except:	
				residnr2.append(int(splitter2.split(tmp2[4])[0]))	#three letter code
			try:													#split is different for one and three letter cod
				if tmp2[2][0] in BASELIST1_T:						#one letter code
					resid1.append(BASELIST3_T[BASELIST1_T.index(tmp2[2][0])])							
				if tmp2[3][0] in BASELIST1_T:
					resid2.append(BASELIST3_T[BASELIST1_T.index(tmp2[3][0])])
			except:
				resid1.append(splitter3.split(tmp2[1])[1])		#three letter code
				resid2.append(splitter3.split(tmp2[1])[3])
		
		self.pairs[infile] = []
		self.pairs[infile].append(intnr)
		self.pairs[infile].append(resid1)
		self.pairs[infile].append(resid2)
		self.pairs[infile].append(residnr1)
		self.pairs[infile].append(residnr2)
		
	def _AverageOnType(self,ntype,intable=None,selrow=0,selrange=None):
		
		"""Average values based on type"""
		
		self.average.clear()
		selection = {}
		
		for n in ntype:
			selection[n] = []
			self.average[n]= []
		
		for structures in self.selected:
			for p in selection:
				for row in intable[structures][selrow]:
					if row == p:
						tmp = []
						for n in selrange:
							tmp.append(intable[structures][n][intable[structures][selrow].index(row)])
						selection[p].append(tmp)
		
		for n in ntype:
			if len(selection[n]) == 0:
				del selection[n]
				del self.average[n]
			else:
				pass
		
		for k in range(0,6):
			for p in selection:
				tmp = []
				for n in selection[p]:
					tmp.append(TransformDash(n[k]))
				if len(tmp) == 1:
					self.average[p].append(tmp[0])
					self.average[p].append(float(0))
					self.average[p].append(len(tmp))
				else:	
					self.average[p].append(mean(tmp))
					self.average[p].append(std(tmp))
					self.average[p].append(len(tmp))
		
	def _AverageOnSequence(self,intable=None,selrange=None,basenr=None,seq=None):
		
		"""Average values based on sequence"""
		
		self.average.clear()
		selection = {}
		self.average['n'] = []
		
		for keys in self.basechainlib.keys():
			for n in range(1,basenr+1):
				self.average[str(n)] = []
				selection[str(n)] = []
		
		for structures in self.selected:
			for p in selection:
				newp = self._MatchBase(p,structures)
				if not newp == False:
					tmp = []
					for n in selrange:
						try:	#Try if data for given structure is present or not
							tmp.append(intable[structures][n][intable[structures][0].index(newp)])
						except:
							pass
					selection[p].append(tmp)
		
		for k in range(0,len(selrange)):
			for p in selection:
				tmp = []
				for n in selection[p]:
					if len(n) < len(selrange):
						pass
					else:
						tmp.append(TransformDash(n[k]))
				if len(tmp) == 0:
					pass
				elif len(tmp) == 1:
					self.average[p].append(tmp[0])
					self.average[p].append(float(0))
				else:	
					self.average[p].append(mean(tmp))
					self.average[p].append(std(tmp))
		
		for p in selection:
			if len(self.average[p]): self.average['n'].append(len(selection[p]))
			else: del self.average[p]      
		
		if seq == 'bp':
			self.average['sequence'] = self._ConvertSeq(self.basesequence[self.chainid][0],'bp')
		elif seq == 'bpstep':
			self.average['sequence'] = self._ConvertSeq(self.basesequence[self.chainid][0],'bpstep')
		else:	
			self.average['sequence'] = self._ConvertSeq(self.basesequence[self.chainid][0],'b')
	
	def _ConvertSeq(self,inlist,btype):
		
		baselist3 = ['ADE','THY','GUA','CYT','URI']
		baselist1 = ['A','T','G','C','U']
		baselist1c = ['T','A','C','G','A']
		pairlist1 = ['A-T','T-A','G-C','C-G','U-A']	
		steplist1 = []
		
		newlist = []
	
		for na in range(len(inlist)):
			if btype == 'b':
				if inlist[na] in baselist3:
					newlist.append(baselist1[baselist3.index(inlist[na])])
			elif btype == 'bp':
				if inlist[na] in baselist3:
					newlist.append(pairlist1[baselist3.index(inlist[na])])
			elif btype == 'bpstep':
				try:
					if inlist[na] in baselist3 and inlist[na+1] in baselist3:
						string = baselist1[baselist3.index(inlist[na])]+baselist1[baselist3.index(inlist[na+1])]+'/'+baselist1c[baselist3.index(inlist[na+1])]+baselist1c[baselist3.index(inlist[na])]
						newlist.append(string)
				except:
					pass		
					
		if len(newlist) > 0:
			return newlist
		else:
			return inlist		
	
	def _MatchBase(self,p,structures):
		
		match = False
		
		try:
			residbase1 = int(self.basechainlib[self.chainid][0][int(p)-1]) 
			residbase2 = int(self.basechainlib[self.chainid][1][int(p)-1])
			
			solbase1 = self.pairs[structures][1][self.pairs[structures][3].index(residbase1)]
			solbase2 = self.pairs[structures][2][self.pairs[structures][4].index(residbase2)]
			
			newp = self.pairs[structures][0][self.pairs[structures][3].index(residbase1)]
			if residbase1 in self.pairs[structures][3]:
				if self.basesequence[self.chainid][0][int(p)-1] == solbase1 and self.basesequence[self.chainid][1][int(p)-1] == solbase2:
					match = newp
			else:
				match = False
		except:
			pass
		
		return match
	
	def _AverageSstvirtb(self,intable=None,selrange=None,basenr=None):
		
		"""Average same strand P-P and C1'-C1' virtual bond distances"""
		
		self.average.clear()
		selection = {}
		self.average['n'] = []
		
		for keys in self.basechainlib.keys():
			for n in range(1,basenr+1):
				self.average[str(n)] = []
				selection[str(n)] = []
				
		for structures in self.selected:
			for p in selection:
				newp = self._MatchBase(p,structures)
				if not newp == False:
					tmp = []
					for n in selrange:
						try:	#Try if data for given structure is present or not
							tmp.append(intable[structures][n][intable[structures][0].index(newp)])
						except:
							pass
					selection[p].append(tmp)
		
		for k in [0,1,4,5]:
			for p in selection:
				tmp = []
				for n in selection[p]:
					if len(n) < 4:
						pass
					else:
						tmp.append(TransformDash(n[k]))
				if len(tmp) == 0:
					pass
				elif len(tmp) == 1:
					self.average[p].append(tmp[0])
					self.average[p].append(float(0))
				else:	
					self.average[p].append(mean(tmp))
					self.average[p].append(std(tmp))
			
		for p in selection:
			if len(self.average[p]): self.average['n'].append(len(selection[p]))
			else: del self.average[p]
		
		self.average['seq1'] = self.basesequence[self.chainid][0]
		self.average['seq2'] = self.basesequence[self.chainid][1]
	
	def _pucker(self,average,std):
		
		"""Classify the sugar puckering based on average phase angle and standard deviation"""
		
		puckerrange = []
		up = int(average)+int(std)
		dw = int(average)-int(std)
		
		if dw in range(0,36):
			puckerrange.append("C3'-endo -")
		elif dw in range(36,72):
			puckerrange.append("C4'-exo -")
		elif dw in range(72,108):
			puckerrange.append("O4'-endo -")
		elif dw in range(108,144):
			puckerrange.append("C1'-exo -")
		elif dw in range(144,180):
			puckerrange.append("C2'-endo -")	
		elif dw in range(180,216) or range(-180,-144):
			puckerrange.append("C3'-exo -")
		elif dw in range(216,252) or range(-144,-108):
			puckerrange.append("C4'-endo -")
		elif dw in range(252,288) or range(-108,-72):
			puckerrange.append("O4'-exo -")
		elif dw in range(288,324) or range(-72,-36):
			puckerrange.append("C1'-endo -")
		elif dw in range(324,360) or range(-36,0):
			puckerrange.append("C2'-exo -")								
	
		if up in range(0,36):
			puckerrange.append("C3'-endo")
		elif up in range(36,72):
			puckerrange.append("C4'-exo")
		elif up in range(72,108):
			puckerrange.append("O4'-endo")
		elif up in range(108,144):
			puckerrange.append("C1'-exo")
		elif up in range(144,180):
			puckerrange.append("C2'-endo")	
		elif up in range(180,216) or range(-180,-144):
			puckerrange.append("C3'-exo")
		elif up in range(216,252) or range(-144,-108):
			puckerrange.append("C4'-endo")
		elif up in range(252,288) or range(-108,-72):
			puckerrange.append("O4'-exo")
		elif up in range(288,324) or range(-72,-36):
			puckerrange.append("C1'-endo")
		elif up in range(324,360) or range(-36,0):
			puckerrange.append("C2'-exo")			
	
		return string.join(puckerrange)
	
	def _UnpairingReport(self):
		
		self.average.clear()
		self.pairsort = {}
	
		for structures in self.pairs:
			match = []
			unpair = 0
			for resid in range(len(self.basechainlib[self.chainid][0])):
				try:
					residbase1 = int(self.basechainlib[self.chainid][0][int(resid)]) 
					residbase2 = int(self.basechainlib[self.chainid][1][int(resid)])
			
					solbase1 = self.pairs[structures][1][self.pairs[structures][3].index(residbase1)]
					solbase2 = self.pairs[structures][2][self.pairs[structures][4].index(residbase2)]
			
					if residbase1 in self.pairs[structures][3]:
						if self.basesequence[self.chainid][0][int(resid)] == solbase1 and self.basesequence[self.chainid][1][int(resid)] == solbase2:
							match.append('  | ')
					else:
						match.append('  X ')
						unpair += 1
				except:
					match.append('  X ')
					unpair += 1
			self.average[structures] = match
			self.pairsort[structures] = unpair
		
	def AutoMaster(self):
	
		"""Set master file based on least number of un-pairing/mispairing events"""
		
		valuesort = sorted(self.pairsort.items(), key=itemgetter(1))
		
		maxpair = valuesort[1]
		selected = []
		for keys in valuesort:
			if keys[1] == maxpair[1]:
				selected.append(keys[0])
		
		self.selected = selected
		
		print ("    * Selecting structures based on minimum unpairing/mispairing events")
		print ("    * %s structures match to a minimum of %s unpairing/mispairing events" % (len(selected),maxpair[1])) 
		
	def ReadOutfiles(self, debug=0):
		
		"""Extract data from tables in all selected 3DNA .out files. Extract with self._ReadTable and
		   self._SortTable"""
		
		rmsd = re.compile("RMSD of the bases")
		origin = re.compile("bp        Ox        Oy        Oz        Nx        Ny        Nz")
		bp = re.compile("bp        Shear    Stretch   Stagger    Buckle  Propeller  Opening")
		bpstep = re.compile("step       Shift     Slide      Rise      Tilt      Roll     Twist")
		helical = re.compile("step       X-disp    Y-disp   h-Rise     Incl.       Tip   h-Twist")
		lamdangl = re.compile("bp     lambda")
		globalp = re.compile("bp       disp.    angle     twist      rise")
		groove = re.compile("P-P     Refined     P-P     Refined")
		angle1 = re.compile("base    alpha    beta   gamma   delta  epsilon   zeta    chi")
		angle2 = re.compile("base       v0      v1      v2      v3      v4      tm       P    Puckering")
		sstvirtb = re.compile("base      P--P")
		helixrad = re.compile("step         P        O4'       C1'        P        O4'        C1'")
		position= re.compile("bp        Px        Py        Pz        Hx        Hy        Hz")
		strand2 = re.compile("Strand II")
		
		for files in self.outfiles:
			readfile = file(files,'r')
			lines = readfile.readlines()	
			linecount = 1
			countstart = []
			for line in lines:
				line = line.strip()
				result1 = bp.match(line)
				result2 = bpstep.match(line)
				result3 = groove.match(line)
				result4 = angle1.match(line)
				result5 = angle2.match(line)
				result6 = strand2.match(line)
				result7 = origin.match(line)
				result8 = helical.match(line)
				result9 = lamdangl.match(line)
				result10 = globalp.match(line)
				result11 = sstvirtb.match(line)
				result12 = helixrad.match(line)
				result13 = position.match(line)
				result14 = rmsd.match(line)
				if result1:
					tablelines = self._ReadTable(linenr=len(countstart),outfile=lines)
					self._SortTable(tablein=tablelines,tableout=self.bp,infile=files)
				elif result14:
					self._MatchPairs(linenr=len(countstart)+2,outfile=lines,infile=files)	
				elif result2:
					tablelines = self._ReadTable(linenr=len(countstart),outfile=lines)
					self._SortTable(tablein=tablelines,tableout=self.bpstep,infile=files)
				elif result3:
					tablelines = self._ReadTable(linenr=len(countstart),outfile=lines)
					self._SortTable(tablein=tablelines,tableout=self.groove,infile=files)
				elif result7:
					tablelines = self._ReadTable(linenr=len(countstart),outfile=lines)
					self._SortTable(tablein=tablelines,tableout=self.origin,infile=files)
				elif result8:
					tablelines = self._ReadTable(linenr=len(countstart),outfile=lines)
					self._SortTable(tablein=tablelines,tableout=self.helical,infile=files)	
				elif result9:
					tablelines = self._ReadTable(linenr=len(countstart),outfile=lines)
					self._SortTable(tablein=tablelines,tableout=self.lamdangl,infile=files)	
				elif result10:
					tablelines = self._ReadTable(linenr=len(countstart),outfile=lines)
					self._SortTable(tablein=tablelines,tableout=self.globalp,infile=files)	
				elif result11:
					tablelines = self._ReadTable(linenr=len(countstart),outfile=lines)
					self._SortTable(tablein=tablelines,tableout=self.sstvirtb,infile=files)	
				elif result12:
					tablelines = self._ReadTable(linenr=len(countstart),outfile=lines)
					self._SortTable(tablein=tablelines,tableout=self.helixrad,infile=files)	
				elif result13:
					tablelines = self._ReadTable(linenr=len(countstart),outfile=lines)
					self._SortTable(tablein=tablelines,tableout=self.position,infile=files)								
				elif not result6:
					try:
						fline = lines[linecount].strip()
						if angle1.match(fline):
							tablelines = self._ReadTable(linenr=len(countstart)+1,outfile=lines)
							self._SortTable(tablein=tablelines,tableout=self.angle11,infile=files)
					except:
						pass
					
					try:
						fline = lines[linecount].strip()
						if angle2.match(fline):
							tablelines = self._ReadTable(linenr=len(countstart)+1,outfile=lines)
							self._SortTable(tablein=tablelines,tableout=self.angle21,infile=files)
					except:
						pass	
				elif result5:
					tablelines = self._ReadTable(linenr=len(countstart),outfile=lines)
					self._SortTable(tablein=tablelines,tableout=self.angle21,infile=files)
				elif result6:
					try:
						fline = lines[linecount].strip()
						if angle1.match(fline):
							tablelines = self._ReadTable(linenr=len(countstart)+1,outfile=lines)
							self._SortTable(tablein=tablelines,tableout=self.angle12,infile=files)
					except:
						pass
					
					try:
						fline = lines[linecount].strip()
						if angle2.match(fline):
							tablelines = self._ReadTable(linenr=len(countstart)+1,outfile=lines)
							self._SortTable(tablein=tablelines,tableout=self.angle22,infile=files)
					except:
						pass	
				linecount += 1
        			countstart.append(linecount)			
	
	def ReadEnerfiles(self):
	
		"""Read EnerPDNA output (base-pair and base-pair step deformation energy)"""
		
		for files in self.outfiles:
			step = re.compile("Total Step Energy")
			pair = re.compile("Total Pair Energy")
			steplines = []
			pairlines = []
			enefile = os.path.splitext(files)[0]+'.ener'
			readfile = file(enefile,'r')
			lines = readfile.readlines()
			
			tlines = len(lines)
			linecount = 0
			Pair = False
			while linecount < tlines:
				if pair.match(lines[linecount].strip()):
					break
				if step.match(lines[linecount].strip()):
					Pair = True
					linecount += 1
				if Pair == False:
					try:
						steplines.append(lines[linecount].strip())
						linecount += 1
					except:
						break
				elif Pair == True:
					try:
						pairlines.append(lines[linecount].strip())
						linecount += 1
					except:
						break
					
			bpstepnr = []
			bpstep = []
			bpstepener = []
			for step in steplines:
				n = step.split()
				bpstepnr.append(n[2])
				bpstep.append(n[3])
				bpstepener.append(float((n[6].split('+'))[1]))
			
			bpnr = []
			bp = []
			bpener = []
			for pair in pairlines:
				p = pair.split()
				bpnr.append(p[2])
				bp.append(p[3])
				bpener.append(float((p[6].split('+'))[1]))
			
			self.bpstepener[files] = []
			self.bpstepener[files].append(bpstepnr)
			self.bpstepener[files].append(bpstep)
			self.bpstepener[files].append(bpstepener)
			
			self.bpener[files] = []
			self.bpener[files].append(bpnr)
			self.bpener[files].append(bp)
			self.bpener[files].append(bpener)	
		
		self.bpenersum = {}
		for files in self.bpener:
			self.bpenersum[files] = (sum(self.bpener[files][2]))
		self.bpstepenersum = {}
		for files in self.bpstepener:
			self.bpstepenersum[files] = (sum(self.bpstepener[files][2]))			
		
	def GetBaseseq(self):
			
		"""Get the base sequence and pairing as it should be in a fully paired structure (input)
		   Calls GetSequence and NAsummery class from QueryPDB plugin. The output must be of a
		   fixed format: the chains in chainlib and sequence in pairs must can only have 2 lists
		   more than two chains are not supported. The two chains must be written as 5'->3' for 
		   the template strand (first list) and 3'->5' for the complementary strand (second list)"""
		
		master = os.path.splitext(os.path.basename(self.outfiles[0]))[0]+'.pdb'
		
		pdb = PDBeditor()
		pdb.ReadPDB(master)	
		xml = pdb.PDB2XML().xml()
		
		sequence = GetSequence()
		sequence.GetSequence(pdbxml=xml)
		
		naeval = NAsummery(pdbxml=xml,sequence=sequence.seqlib)
		naeval.Evaluate()	
		
		self.basemoltype = naeval.moltype
		self.basechainlib = naeval.chainlib
		self.basesequence = naeval.pairs
	
		for keys in self.basechainlib.keys():
			self.basenr += len(self.basechainlib[keys][0])
			self.chainid = keys
		
	def GetMasterSeq(self, master=None):
		
		"""Use excact match with supplied master (.out) file as selection criterium if a 
		   master file is supplied. A excact match is stored in self.selected non-matches 
		   in self.rejected. If no master file is supplied all files are selected and 
		   stored in self.selected"""
	  
		if master == None:
			self.selected = self.outfiles
			self.master = None
		else:	
			currdir = os.getcwd()
			master = os.path.join(currdir, master)
			self.chainid = 'B'	#Dummy chainid, type does not matter
			
			self.basechainlib.clear()
			self.basesequence.clear()
			
			self.basechainlib[self.chainid] = []
			self.basechainlib[self.chainid].append(self.pairs[master][3])
			self.basechainlib[self.chainid].append(self.pairs[master][4])
			
			self.basesequence[self.chainid] = []
			self.basesequence[self.chainid].append(self.pairs[master][1])
			self.basesequence[self.chainid].append(self.pairs[master][2])
			
			self.basenr = len(self.basechainlib[self.chainid][0])
			
			self._MatchMaster()
			self.master = master
				
	def WriteStats(self):
	
		"""Write all statistics to the multiout.stat file"""
			
		outfile = file('multiout.stat','w')
		
		outfile.write('*************************************************************************************************************************\n')
		outfile.write('Multi-structure analysis for %i structures\n' % (len(self.selected)+len(self.rejected)))
		outfile.write('Time and date: %s\n' % ctime())
		outfile.write("\n1. The list of the parameters given below correspond to the 5' to 3' direction\n")
		outfile.write("  of strand I and 3' to 5' direction of strand II.\n")
		outfile.write("\n2. All angular parameters, except for the phase angle of sugar pseudo-\n")
		outfile.write("   rotation, are measured in degrees in the range of [-180, +180], and all\n")
		outfile.write("   displacements are measured in Angstrom units.\n")
		outfile.write('*************************************************************************************************************************\n')
		outfile.write('Structure sequence information:\n\n')
		
		if not len(self.basemoltype) == 0:
			chains = self.basemoltype.keys()
			for chain in xrange(len(chains)):
				try:
					outfile.write('Found chain: %s of type %s\n' % (chains[chain],self.basemoltype[chains[chain]]))
					if self.basemoltype[chains[chain]] == 'PROT':
						pass
					else:
						outfile.write("Sequence of chain %s: 5'" % chains[chain])	
						for resid in self.basesequence[chains[chain]][0]:
							outfile.write(' %s' % resid)
						outfile.write(" 3'\n")	
				except:
					if len(self.basesequence[chains[chain-1]]) == 2:
						for resid in self.basesequence[chains[chain-1]][1]:
							outfile.write(' %s' % resid)
						outfile.write(" 3'\n")		
							
		else:
			outfile.write("Nucleic acid sequence: 5'")	
			for resid in self.basesequence[self.chainid][1]:
				outfile.write(' %s' % resid)
			outfile.write(" 3'\n")				
		
		if not self.master == None:
			outfile.write('\nSelected and rejected files based on a match in base-pairing with the master file: %s\n' % os.path.basename(self.master))
			outfile.write('\n%i structures were selected:\n' % (len(self.selected)))
		
			for n in self.selected:
				outfile.write('%s\n' % (os.path.basename(n)))
		
			outfile.write('\n%i structures were rejected:\n' % (len(self.rejected)))
		
			for n in self.rejected:
				outfile.write('%s\n' % (os.path.basename(n)))
		
		outfile.write('\n*************************************************************************************************************************\n')
		outfile.write('Statistical data for each type of local base-pair parameter in %i equal structures (average/standard deviation)\n' % (len(self.selected)))
		outfile.write('\nbp      freq.      shear          stretch         stagger         buckle          proptw         opening\n')
		
		self._AverageOnType(ntype=['A-T','T-A','G-C','C-G'],intable=self.bp,selrow=1,selrange=range(2,8))
		
		for ntype in self.average:
			outfile.write('%1s%8i%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n' % 
			((ntype),self.average[ntype][2],self.average[ntype][0],self.average[ntype][1],self.average[ntype][3],self.average[ntype][4],self.average[ntype][6],
			self.average[ntype][7],self.average[ntype][9],self.average[ntype][10],self.average[ntype][12],self.average[ntype][13],self.average[ntype][15],
			self.average[ntype][16]))
		
		outfile.write('\n*************************************************************************************************************************\n')
		outfile.write('Statistical data for each type of local base-pair step parameter in %i equal structures (average/standard deviation)\n' % (len(self.selected)))
		outfile.write('\nbpstep  freq.      shift          slide            rise            tilt            roll           twist\n')
		
		self._AverageOnType(ntype=PAIRSTEPS,intable=self.bpstep,selrow=1,selrange=range(2,8))
		
		for ntype in self.average:
			outfile.write('%1s%6i%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n' % 
			((ntype),self.average[ntype][2],self.average[ntype][0],self.average[ntype][1],self.average[ntype][3],self.average[ntype][4],self.average[ntype][6],
			self.average[ntype][7],self.average[ntype][9],self.average[ntype][10],self.average[ntype][12],self.average[ntype][13],self.average[ntype][15],
			self.average[ntype][16]))
		
		outfile.write('\n*************************************************************************************************************************\n')
		outfile.write('Statistical data for origin (Ox, Oy, Oz) and mean normal vector (Nx, Ny, Nz) of each base-pair in\n') 
		outfile.write('the coordinate system of the given structure in %i equal structures (average/standard deviation)\n' % (len(self.selected)))
		outfile.write('\nnr    bp     freq.            Ox              Oy              Oz              Nx              Ny              Nz\n')
		
		self._AverageOnSequence(intable=self.origin,selrange=range(2,8),basenr=self.basenr,seq='bp')
		
		baselist = []
		for n in self.average:
			try: baselist.append(int(n))
			except: pass
		
		count= min(baselist)
		while count < len(self.average)-1:
			seqnr = str(count)
			outfile.write('%2d%6s%6s%10.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n' % 
			((int(seqnr)),self.average['sequence'][count-1],self.average['n'][count-1],self.average[seqnr][0],self.average[seqnr][1],self.average[seqnr][2],
			self.average[seqnr][3],self.average[seqnr][4],self.average[seqnr][5],self.average[seqnr][6],self.average[seqnr][7],self.average[seqnr][8],
			self.average[seqnr][9],self.average[seqnr][10],self.average[seqnr][11]))
			count = count+1
		
		outfile.write('\n*************************************************************************************************************************\n')
		outfile.write('Statistical data for each local base-pair parameter in %i equal structures (average/standard deviation)\n' % (len(self.selected)))
		outfile.write('\nnr   bp      freq.        shear           stretch       stagger          buckle           proptw         opening\n')
		
		self._AverageOnSequence(intable=self.bp,selrange=range(2,8),basenr=self.basenr,seq='bp')
		
		baselist = []
		for n in self.average:
			try: baselist.append(int(n))
			except: pass
		
		count= min(baselist)
		while count < len(self.average)-1:
			seqnr = str(count)
			outfile.write('%2d%6s%6s%10.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n' % 
			((int(seqnr)),self.average['sequence'][count-1],self.average['n'][count-1],self.average[seqnr][0],self.average[seqnr][1],self.average[seqnr][2],self.average[seqnr][3],
			self.average[seqnr][4],self.average[seqnr][5],self.average[seqnr][6],self.average[seqnr][7],self.average[seqnr][8],self.average[seqnr][9],
			self.average[seqnr][10],self.average[seqnr][11]))
			count = count+1
		
		outfile.write('\n*************************************************************************************************************************\n')
		outfile.write('Statistical data for each loacl base-pair step parameter in %i equal structures (average/standard deviation)\n' % (len(self.selected)))
		outfile.write('\nnr   bpstep    freq.      shift           slide           rise            tilt             roll           twist\n')
		
		self._AverageOnSequence(intable=self.bpstep,selrange=range(2,8),basenr=self.basenr-1,seq='bpstep')
		
		baselist = []
		for n in self.average:
			try: baselist.append(int(n))
			except: pass
		
		count= min(baselist)
		while count < len(self.average)-1:
			seqnr = str(count)
			outfile.write('%2d%8s%8s%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n' % 
			((int(seqnr)),self.average['sequence'][count-1],self.average['n'][count-1],self.average[seqnr][0],self.average[seqnr][1],self.average[seqnr][2],self.average[seqnr][3],
			self.average[seqnr][4],self.average[seqnr][5],self.average[seqnr][6],self.average[seqnr][7],self.average[seqnr][8],self.average[seqnr][9],
			self.average[seqnr][10],self.average[seqnr][11]))
			count = count+1
		
		outfile.write('\n*************************************************************************************************************************\n')
		outfile.write('Statistical data for each local base-pair helical parameter in %i equal structures (average/standard deviation)\n' % (len(self.selected)))
		outfile.write('\nnr   bpstep     freq.      X-disp           Y-disp         h-rise            Incl.             Tip         h-twist\n')
		
		self._AverageOnSequence(intable=self.helical,selrange=range(2,8),basenr=self.basenr-1,seq='bpstep')
		
		baselist = []
		for n in self.average:
			try: baselist.append(int(n))
			except: pass
		
		count= min(baselist)
		while count < len(self.average)-1:
			seqnr = str(count)
			outfile.write('%2d%8s%8s%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n' % 
			((int(seqnr)),self.average['sequence'][count-1],self.average['n'][count-1],self.average[seqnr][0],self.average[seqnr][1],self.average[seqnr][2],self.average[seqnr][3],
			self.average[seqnr][4],self.average[seqnr][5],self.average[seqnr][6],self.average[seqnr][7],self.average[seqnr][8],self.average[seqnr][9],
			self.average[seqnr][10],self.average[seqnr][11]))
			count = count+1
		
		outfile.write('\n*************************************************************************************************************************\n')
		outfile.write('Statistical data for each lambda viryal angle in %i equal structures (average/standard deviation)\n' % (len(self.selected)))
		outfile.write("lambda: virtual angle between C1'-YN1 or C1'-RN9 glycosidic bonds and the\n")
		outfile.write("        base-pair C1'-C1' line\n")
		outfile.write("\nC1'-C1': distance between C1' atoms for each base-pair\n")
		outfile.write("RN9-YN1: distance between RN9-YN1 atoms for each base-pair\n")
		outfile.write("RC8-YC6: distance between RC8-YC6 atoms for each base-pair\n")
		outfile.write("\nnr   bp     freq.      lambda(I)      lambda(II)         C1'-C1'        RN9-YN1         RC8-YC6\n")
		
		self._AverageOnSequence(intable=self.lamdangl,selrange=range(2,7),basenr=self.basenr,seq='bp')
		
		baselist = []
		for n in self.average:
			try: baselist.append(int(n))
			except: pass
		
		count= min(baselist)
		while count < len(self.average)-1:
			seqnr = str(count)
			outfile.write('%2d%6s%6s%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n' % 
			((int(seqnr)),self.average['sequence'][count-1],self.average['n'][count-1],self.average[seqnr][0],self.average[seqnr][1],self.average[seqnr][2],self.average[seqnr][3],
			self.average[seqnr][4],self.average[seqnr][5],self.average[seqnr][6],self.average[seqnr][7],self.average[seqnr][8],self.average[seqnr][9]))
			count = count+1
		
		outfile.write('\n*************************************************************************************************************************\n')
		outfile.write('Statistical data for the minor and major groove width in %i equal structures (average/standard deviation)\n' % (len(self.selected)))
		outfile.write('\nMinor and major groove widths: direct P-P distances and refined P-P distances\n')
		outfile.write('  which take into account the directions of the sugar-phosphate backbones\n')
		outfile.write('\n  (Subtract 5.8 Angstrom from the values to take account of the vdw radii\n')
		outfile.write('  of the phosphate groups, and for comparison with FreeHelix and Curves.)\n')
		outfile.write('\nRef: M. A. El Hassan and C. R. Calladine (1998). ``Two Distinct Modes of\n')
		outfile.write('  Protein-induced Bending in DNA.'' J. Mol. Biol., v282, pp331-343.\n')
		outfile.write('\n                           Minor Groove                     Major Groove\n')
		outfile.write('nr   basepair   freq.    (P-P)         (Refined)         (P-P)         (Refined)\n')
		
		self._AverageOnSequence(intable=self.groove,selrange=range(2,6),basenr=self.basenr-1,seq='bpstep')
		
		baselist = []
		for n in self.average:
			try: baselist.append(int(n))
			except: pass
		
		count= min(baselist)
		while count < len(self.average)-1:
			seqnr = str(count)
			outfile.write('%2d%8s%8s%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n' % 
			((int(seqnr)),self.average['sequence'][count-1],self.average['n'][count-1],self.average[seqnr][0],self.average[seqnr][1],self.average[seqnr][2],self.average[seqnr][3],
			self.average[seqnr][4],self.average[seqnr][5],self.average[seqnr][6],self.average[seqnr][7]))
			count = count+1
		
		outfile.write('\n*************************************************************************************************************************\n')
		outfile.write("Statistical data for global parameters based on C1'-C1' vectors in %i equal structures (average/standard deviation)\n" % (len(self.selected)))
		outfile.write("\ndisp.: displacement of the middle C1'-C1' point from the helix\n")
		outfile.write("angle: inclination between C1'-C1' vector and helix (subtracted from 90)\n")
		outfile.write("twist: helical twist angle between consecutive C1'-C1' vectors\n")
		outfile.write("rise : helical rise by projection of the vector connecting consecutive\n")
		outfile.write("       C1'-C1' middle points onto the helical axis\n")
		outfile.write('\nnr   bp      freq.       disp.          angle            twist           rise\n')
		
		self._AverageOnSequence(intable=self.globalp,selrange=range(2,6),basenr=self.basenr,seq='bp')
		
		baselist = []
		for n in self.average:
			try: baselist.append(int(n))
			except: pass
		
		if len(baselist): 
			count= min(baselist)
			while count < len(self.average)-1:
				seqnr = str(count)
				outfile.write('%2d%6s%6s%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n' % 
				((int(seqnr)),self.average['sequence'][count-1],self.average['n'][count-1],self.average[seqnr][0],self.average[seqnr][1],self.average[seqnr][2],self.average[seqnr][3],
				self.average[seqnr][4],self.average[seqnr][5],self.average[seqnr][6],self.average[seqnr][7]))
				count = count+1
		else:
			outfile.write('\n	No data available\n')
		
		outfile.write('\n*************************************************************************************************************************\n')
		outfile.write('Statistical data for main chain and chi torsion angles in %i equal structures (average/standard deviation)\n' % (len(self.selected)))
		outfile.write("\nNote: alpha:   O3'(i-1)-P-O5'-C5'\n")
		outfile.write("      beta:    P-O5'-C5'-C4'\n")
		outfile.write("      gamma:   O5'-C5'-C4'-C3'\n")
		outfile.write("      delta:   C5'-C4'-C3'-O3'\n")
		outfile.write("      epsilon: C4'-C3'-O3'-P(i+1)\n")
		outfile.write("      zeta:    C3'-O3'-P(i+1)-O5'(i+1)\n")
		outfile.write("\n      chi for pyrimidines(Y): O4'-C1'-N1-C2\n")
		outfile.write("      chi for purines(R): O4'-C1'-N9-C4\n")
		outfile.write('\nStrand I\n')
		outfile.write('nr   base    freq.    alpha           beta            gamma           delta         epsilon           zeta            chi\n')
		
		self._AverageOnSequence(intable=self.angle11,selrange=range(2,9),basenr=self.basenr,seq='b')
		
		baselist = []
		for n in self.average:
			try: baselist.append(int(n))
			except: pass
		
		count= min(baselist)
		while count < len(self.average)-1:
			seqnr = str(count)
			outfile.write('%2d%6s%6s%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n' % 
			((int(seqnr)),self.average['sequence'][count-1],self.average['n'][count-1],self.average[seqnr][0],self.average[seqnr][1],self.average[seqnr][2],self.average[seqnr][3],
			self.average[seqnr][4],self.average[seqnr][5],self.average[seqnr][6],self.average[seqnr][7],self.average[seqnr][8],self.average[seqnr][9],
			self.average[seqnr][10],self.average[seqnr][11],self.average[seqnr][12],self.average[seqnr][13]))
			count = count+1
		
		outfile.write('\nStrand II\n')
		outfile.write('nr   base    freq.     alpha           beta            gamma           delta         epsilon           zeta            chi\n')
		
		self._AverageOnSequence(intable=self.angle12,selrange=range(2,9),basenr=self.basenr,seq='b')
		
		baselist = []
		for n in self.average:
			try: baselist.append(int(n))
			except: pass
		
		count= min(baselist)
		while count < len(self.average)-1:
			seqnr = str(count)
			outfile.write('%2d%6s%6s%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n' % 
			((int(seqnr)),self.average['sequence'][count-1],self.average['n'][count-1],self.average[seqnr][0],self.average[seqnr][1],self.average[seqnr][2],self.average[seqnr][3],
			self.average[seqnr][4],self.average[seqnr][5],self.average[seqnr][6],self.average[seqnr][7],self.average[seqnr][8],self.average[seqnr][9],
			self.average[seqnr][10],self.average[seqnr][11],self.average[seqnr][12],self.average[seqnr][13]))
			count = count+1
		
		outfile.write('\n*************************************************************************************************************************\n')
		outfile.write('Statistical data for sugar conformational parameters in %i equal structures (average/standard deviation)\n' % (len(self.selected)))
		outfile.write("\nNote: v0: C4'-O4'-C1'-C2'\n")
		outfile.write("      v1: O4'-C1'-C2'-C3'\n")
		outfile.write("      v2: C1'-C2'-C3'-C4'\n")
		outfile.write("      v3: C2'-C3'-C4'-O4'\n")
		outfile.write("      v4: C3'-C4'-O4'-C1'\n")
		outfile.write("\n      tm: amplitude of pseudorotation of the sugar ring\n")
		outfile.write("       P:  phase angle of pseudorotation of the sugar ring\n")
		outfile.write('\nStrand I\n')
		outfile.write('nr   base    freq.     v0            v1            v2            v3            v4            tm            P        Puckering (max-min)\n')
		
		self._AverageOnSequence(intable=self.angle21,selrange=range(2,9),basenr=self.basenr,seq='b')
		
		baselist = []
		for n in self.average:
			try: baselist.append(int(n))
			except: pass
		
		count= min(baselist)
		while count < len(self.average)-1:
			seqnr = str(count)
			outfile.write('%2d%6s%6s%8.2f%6.2f%8.2f%6.2f%8.2f%6.2f%8.2f%6.2f%8.2f%6.2f%8.2f%6.2f%8.2f%6.2f%22s\n' % 
			((int(seqnr)),self.average['sequence'][count-1],self.average['n'][count-1],self.average[seqnr][0],self.average[seqnr][1],self.average[seqnr][2],self.average[seqnr][3],
			self.average[seqnr][4],self.average[seqnr][5],self.average[seqnr][6],self.average[seqnr][7],self.average[seqnr][8],self.average[seqnr][9],
			self.average[seqnr][10],self.average[seqnr][11],self.average[seqnr][12],self.average[seqnr][13],self._pucker(self.average[seqnr][12],self.average[seqnr][13])))
			count = count+1
		
		outfile.write('\nStrand II\n')
		outfile.write('nr   base    freq.     v0            v1            v2            v3            v4            tm            P        Puckering (max-min)\n')
		
		self._AverageOnSequence(intable=self.angle22,selrange=range(2,9),basenr=self.basenr,seq='b')
		
		baselist = []
		for n in self.average:
			try: baselist.append(int(n))
			except: pass
		
		count= min(baselist)
		while count < len(self.average)-1:
			seqnr = str(count)
			outfile.write('%2d%6s%6s%8.2f%6.2f%8.2f%6.2f%8.2f%6.2f%8.2f%6.2f%8.2f%6.2f%8.2f%6.2f%8.2f%6.2f%22s\n' % 
			((int(seqnr)),self.average['sequence'][count-1],self.average['n'][count-1],self.average[seqnr][0],self.average[seqnr][1],self.average[seqnr][2],self.average[seqnr][3],
			self.average[seqnr][4],self.average[seqnr][5],self.average[seqnr][6],self.average[seqnr][7],self.average[seqnr][8],self.average[seqnr][9],
			self.average[seqnr][10],self.average[seqnr][11],self.average[seqnr][12],self.average[seqnr][13],self._pucker(self.average[seqnr][12],self.average[seqnr][13])))
			count = count+1
		
		outfile.write('\n*************************************************************************************************************************\n')
		outfile.write("Statistical data for same strand P-P and C1'-C1' virtual bond distances in %i equal structures (average/standard deviation)\n" % (len(self.selected)))
		outfile.write("\n                          Strand I                                  Strand II\n")
		outfile.write("nr   base     freq.      P--P           C1'-C1'      nr   base       P--P           C1'-C1'\n")
		
		self._AverageSstvirtb(intable=self.sstvirtb,selrange=range(2,8),basenr=self.basenr-1)
		
		baselist = []
		for n in self.average:
			try: baselist.append(int(n))
			except: pass
		
		count= min(baselist)
		while count < len(self.average)-2:
			seqnr = str(count)
			outfile.write('%2d%6s%6s%8.2f%8.2f%8.2f%8.2f%6d%6s%8.2f%8.2f%8.2f%8.2f\n' % 
			((int(seqnr)),self.average['seq1'][count-1],self.average['n'][count-1],self.average[seqnr][0],self.average[seqnr][1],self.average[seqnr][2],self.average[seqnr][3],
			(int(seqnr)),self.average['seq2'][count-1],self.average[seqnr][4],self.average[seqnr][5],self.average[seqnr][6],self.average[seqnr][7]))
			count = count+1
		
		outfile.write('\n*************************************************************************************************************************\n')
		outfile.write("Statistical data for helix radius (radial displacement of P, O4', and C1' atoms in local helix frame of each dimer) in\n") 
		outfile.write('%i equal structures (average/standard deviation)\n' % (len(self.selected)))
		outfile.write('\n                                  Strand I                                     Strand II\n')              
		outfile.write("nr   bpstep    freq.          P               O4'             C1'             P              O4'              C1'\n")
		
		self._AverageOnSequence(intable=self.helixrad,selrange=range(2,8),basenr=self.basenr-1,seq='bpstep')
		
		baselist = []
		for n in self.average:
			try: baselist.append(int(n))
			except: pass
		
		count= min(baselist)
		while count < len(self.average)-1:
			seqnr = str(count)
			outfile.write('%2d%8s%6s%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n' % 
			((int(seqnr)),self.average['sequence'][count-1],self.average['n'][count-1],self.average[seqnr][0],self.average[seqnr][1],self.average[seqnr][2],self.average[seqnr][3],
			self.average[seqnr][4],self.average[seqnr][5],self.average[seqnr][6],self.average[seqnr][7],self.average[seqnr][8],self.average[seqnr][9],
			self.average[seqnr][10],self.average[seqnr][11]))
			count = count+1
		
		outfile.write('\n*************************************************************************************************************************\n')
		outfile.write('Statistical data for position (Px, Py, Pz) and local helical axis vector (Hx, Hy, Hz)  for each dinucleotide step\n')
		outfile.write('in %i equal structures (average/standard deviation)\n' % (len(self.selected)))
		outfile.write('\nnr   bp      freq           Px              Py              Pz              Hx              Hy              Hz\n')
		
		self._AverageOnSequence(intable=self.position,selrange=range(2,8),basenr=self.basenr-1,seq='bpstep')
		
		baselist = []
		for n in self.average:
			try: baselist.append(int(n))
			except: pass
		
		count= min(baselist)
		while count < len(self.average)-1:
			seqnr = str(count)
			outfile.write('%2d%8s%6s%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n' % 
			((int(seqnr)),self.average['sequence'][count-1],self.average['n'][count-1],self.average[seqnr][0],self.average[seqnr][1],self.average[seqnr][2],self.average[seqnr][3],
			self.average[seqnr][4],self.average[seqnr][5],self.average[seqnr][6],self.average[seqnr][7],self.average[seqnr][8],self.average[seqnr][9],
			self.average[seqnr][10],self.average[seqnr][11]))
			count = count+1	
		
		outfile.write('\n*************************************************************************************************************************\n')
		
		outfile.close()	
	
	def PairStats(self):
		
		"""Writing nucleic acid pairing information to the file 'napairing.stat'"""
		
		outfile = file('napairing.stat','w')
		
		outfile.write('\n*************************************************************************************************************************\n')
		outfile.write('Nucleic acid (un)-pairing for %i structures\n' % (len(self.selected)+len(self.rejected)))
		outfile.write('Time and date: %s\n' % ctime())
		outfile.write('Sequences are compared to the full length paired sequence or the supplied master file. Unpaired\n') 
		outfile.write('or mispaired bases are shown with an X. Consult the associated .out file(s) for more information about the nature of the\n')
		outfile.write('unpairing or mispairing. Structures are sorted from low to high number of unpairing/mispatching events\n') 
		outfile.write('*************************************************************************************************************************\n') 
		 
		self._UnpairingReport()
		valuesort = sorted(self.pairsort.items(), key=itemgetter(1))
		
		for key in valuesort:
			outfile.write("\nStructure %s, total bp energy %0.2f, total bp-step energy %0.2f, unpairing/mispairing %i:\n" % 
			(os.path.basename(key[0]),self.bpenersum[key[0]],self.bpstepenersum[key[0]],key[1]))
			outfile.write(" 5'-")	
			for n in self.basesequence[self.chainid][0]:
				outfile.write(' %s' % n)
			outfile.write(" -3'\n")
			outfile.write("    ")	
			for n in self.average[key[0]]:
				outfile.write('%s' % n)
			outfile.write("\n 5'-")	
			for n in self.basesequence[self.chainid][1]:
				outfile.write(' %s' % n)		
			outfile.write(" -3'\n")	
		
		outfile.write('\n*************************************************************************************************************************\n')
		
		outfile.close()
		
	def FileList(self):
		
		"""Write selected files to selection.list"""
		
		outfile = file('selection.list','w')
		
		for n in self.selected:
			outfile.write("%s\n" % n)
		
		outfile.close()
	
if __name__ == '__main__':

	"""Running this plugin from the command line"""
	
	"""Import command line specific modules"""
	from optparse import *
	
	"""Setting up parameter dictionary"""
	paramdict = CommandlineOptionParser().option_dict
	
	"""Check for input"""
	if paramdict['input'] == None:
		print "    * Please supply pdb file using option -f or use option -h/--help for usage"
		sys.exit(0)
	else:
		inputlist = paramdict['input']
	
	"""Running plugin main functions"""
	PluginCore(paramdict, inputlist)
	sys.exit(0)

	
