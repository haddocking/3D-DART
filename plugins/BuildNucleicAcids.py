#!/usr/bin/env python2.7

USAGE = """
==========================================================================================

Author:				Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center
					for Biomolecular Research, Utrecht university, The Netherlands
Copyright (C):		2006 (DART project)
DART version:		1.2  (25-11-2008)
DART plugin: 		BuildNucleicAcids.py
Input:				A 3DNA PAR data file or a sequence and any of the allowed options,
					any order and combination.
Output:				A PDB or alchemy file of the generated structure.
Plugin excecution:	Either command line driven (use -h/--help for the option) or as 
					part of a DART batch sequence.	
Plugin function:	Build a nucleic acid structure of a user defined type and sequence 
					(uses the "fiber" module of 3DNA) or rebuild a nucleic acid 
					structure from a base-pair(step) parameter file in various modes 
					(uses the "rebuild" module of 3DNA).
Examples:			BuildNucleicAcids -f test.par -at ADNA
					BuildNucleicAcids -s AAGTCGGTC -dn test.alc
Dependencies:		Standard python2.3 or higher modules, 3DNA

==========================================================================================
"""

"""Import modules"""
import os, sys, glob, shutil, commands

"""Setting pythonpath variables if run from the command line"""
base, dirs = os.path.split(os.path.dirname(os.path.join(os.getcwd(), __file__)))

if base in sys.path: pass
else: sys.path.append(base)

def PluginXML():
	
	PluginXML = """ 
<metadata>
 <name>Build Nucleic Acid structure from base-pair/base-pair step parameter file</name>
 <input type="Filetype">.par</input>
 <output type="Filetype">.pdb,.alc</output>
</metadata>
<parameters>
 <option type="useplugin" form="hidden" text="None">True</option>
 <option type="inputfrom" form="hidden" text="None">1</option>
 <option type="name" form="text" text="Give your structure a name">dna</option>
 <option type="sequence" form="text" text="Provide the nucleic acid sequence (5'->3')"></option>  
 <option type="repeat" form="text" text="Number of repeats">1</option>  
 <option type="type" form="list" default="BDNA,ADNA,NDB96" text="Nucleic acid type">BDNA</option>
 <option type="listfiber" form="checkbox" text="Show a list of the 55 fiber models">False</option>
 <option type="atomic" form="checkbox" text="Full atomic model in PDB format">True</option>
 <option type="basep" form="checkbox" text="Only base and P atoms in PDB format">False</option>
 <option type="block1" form="checkbox" text="One block per base-pair/base in ALCHEMY format">False</option>
 <option type="block2" form="checkbox" text="Two blocks per base-pair in ALCHEMY format">False</option>
 <option type="negx" form="checkbox" text="reverse the direction of x- and z-axes (for Z-DNA)">False</option>
</parameters>"""
	
	return PluginXML

def PluginCore(paramdict, inputlist):
	
	print "--> Starting BuildNucleicAcids"
	
	if paramdict['name']: paramdict['name'] = ''.join(paramdict['name'].split())	#Strip any whitespaces from the name	
	
	if paramdict['listfiber']: ListFiber()	#Only list the 55 available fiber models	
	
	if paramdict['sequence']:		

		#Building nucleic acid structure from scratch, This is basically a wrapper around the "fiber" command
		#of 3DNA.
		validseq = ''; nonvalid = ''
		for base in paramdict['sequence']:
			if base.upper() in ['A','C','T','G']: validseq += base.upper()
			else: nonvalid += base	 
		
		if len(validseq): 
			paramdict['sequence'] = validseq
			if len(nonvalid):
				print "    * The sequence contained the following non-valid bases: %s. They were removed" % nonvalid
			print "    * Building structure from scratch with sequence: %s" % paramdict['sequence']
			print "    * Repeating this sequence %i times" % paramdict['repeat'] 	
			print "    * Building nucleic acid as type: %s" % paramdict['type']
			print "    * Structure will be given the name: %s.pdb" % paramdict['name']
		
			FiberModule(paramdict['sequence'], paramdict['repeat'], paramdict['type'], paramdict['name']) 
		else:
			print "    * ERROR: the complete sequence '%s' is not valid. Stopping" % paramdict['sequence']
			sys.exit()
		
	elif inputlist:

		#Building nucleic acid structure from .par file, This is basically a wrapper around the "rebuild" command
		#of 3DNA.
		
		for inputfile in inputlist:
			
			if not paramdict['name']: basename, extension = os.path.splitext(os.path.basename(inputfile))
			else: basename, extension = os.path.splitext(paramdict['name'])
			
			if paramdict['type'] in ['BDNA','ADNA','NDB96']: natype = paramdict['type']
			elif paramdict['type'] == 'custom': natype = 'custom'
			else: natype = 'BDNA'

			option = None
			
			#Building structures in a atomic representation
			if paramdict['atomic']:
				print "    * Building full atomic model"	
				option = '-atomic'
				outputfile = basename+".pdb"	
			elif paramdict['basep']:
				print "    * Building model with only base and P atoms"
				option = '-base_p'
				outputfile = basename+".pdb"	
			else:
				pass
			
			if option:
				print "    * Rebuilding a structure from supplied par file:", os.path.basename(inputfile)
				if paramdict['negx']:
					print "    * reverse the direction of x- and z-axes (for Z-DNA)"
					option = option + ' -negx'
				if natype == 'custom':
					print "    * 'custom' option selected for nucleic acid type. Looking for supplied pdb files"
					use_Custom(curdir=os.getcwd())
				else:
					get_Atomic(natype)
				
				RebuildNA(option, inputfile, outputfile)
			
			option = None
			
			#Building structures representations in a ALCHEMY format
			if paramdict['block1']:
				print "    * Building model with one block per base-pair/base in ALCHEMY format"
				option = '-block1'
				outputfile = basename+".alc"
			elif paramdict['block2']:
				print "    * Building model with two block per base-pair/base in ALCHEMY format"
				option = '-block2'
				outputfile = basename+".alc"
			else:
				pass
			
			if option:
				print "    * Rebuilding a structure from supplied par file:", os.path.basename(inputfile)
				if paramdict['negx']:
					print "    * reverse the direction of x- and z-axes (for Z-DNA)"
					option = option + ' -negx'
				if natype == 'custom':
					print "    * 'custom' option selected for nucleic acid type. Looking for supplied pdb files"
					use_Custom(curdir=os.getcwd())
				else:
					get_Atomic(natype)	
				
				RebuildNA(option, inputfile, outputfile)
				
	else:
		print "    * ERROR: You have neither given a sequence to build or a .par file to rebuild. Stopping"
		sys.exit(0)		

	CleanUp()

#================================================================================================================================#
# 										PLUGIN SPECIFIC DEFINITIONS BELOW THIS LINE												 #
#================================================================================================================================#

class CommandlineOptionParser:
	
	"""Parses command line arguments using optparse"""
	
	def __init__(self):
		
		self.option_dict = {}
		self.option_dict = self.CommandlineOptionParser()
	
	def CommandlineOptionParser(self):
	
		"""Parsing command line arguments"""
	
		usage = "usage: %prog" + USAGE
		parser = OptionParser(usage=usage)
		
		typeusage = """Define the type of nucleic acid you whant your structure to be 
                       build or rebuild with. The options are: 
			       
			       standard = BDNA, ADNA, ZDNA or CDNA for build 
			       	          BDNA, ADNA or NDB96 for rebuild 
			       fiber    = Any of the 55 fiber models can be used for building
			                  a structure. Use the -l option to get a list of models
			       custom   = Supply your own pdb files for te individual G,C,T,A and U
			       	          bases in the CUSTOM directory of X3DNA at /DART/system/
					  		  third-party/X3DNA/. Supply the files in the form "A.pdb" ...
			       """	
		
		parser.add_option( "-f", "--file", action="callback", callback=self.varargs, dest="inputfile", type="string", help="Supply par inputfile(s)")
		parser.add_option( "-n", "--name", action="store", dest="name", default='dna', help="Define the name for your structure")
		parser.add_option( "-s", "--sequence", action="store", dest="sequence", help="Enter the desired sequence")
		parser.add_option( "-r", "--repeat", action="store", dest="repeat", default=1, type="int", help="State how many times you would like to repeat the sequence (default=1)")
		parser.add_option( "-t", "--type", action="store", dest="type", default='BDNA', help=typeusage)
		parser.add_option( "-l", "--listfiber", action="store_true", dest="listfiber", default=False, help="Display a short description of al 55 fiber models")
		parser.add_option( "-a", "--atomic", action="store_true", dest="atomic", default=True, help="full atomic model in PDB format")
		parser.add_option( "-b", "--basep", action="store_true", dest="basep", default=False, help="with only base and P atoms in PDB format")
		parser.add_option( "-c", "--block1", action="store_true", dest="block1", default=False, help="one block per base-pair/base in ALCHEMY format (default)")
		parser.add_option( "-d", "--block2", action="store_true", dest="block2", default=False, help="two blocks per base-pair in ALCHEMY format")
		parser.add_option( "-e", "--negx", action="store_true", dest="negx", default=False, help="reverse the direction of x- and z-axes (for Z-DNA)")
		
		(options, args) = parser.parse_args()
		
		self.option_dict['input'] = options.inputfile
		self.option_dict['name'] = options.name
		self.option_dict['sequence'] = options.sequence
		self.option_dict['repeat'] = options.repeat
		self.option_dict['type'] = options.type
		self.option_dict['listfiber'] = options.listfiber
		self.option_dict['atomic'] = options.atomic
		self.option_dict['basep'] = options.basep
		self.option_dict['block1'] = options.block1
		self.option_dict['block2'] = options.block2
		self.option_dict['negx'] = options.negx
			
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

def ListFiber():
	
	"""Print list to screen of all 55 available fiber structures""" 
	
	cmd = "fiber -l" 
	os.system(cmd)	

def FiberModule(sequence, repeat, natype, name):
	
	"""Make nucleic acid structures from using on of the 55 available fiber structures"""
	
	try:
		x = int(natype)
		if x in range(1,56): option = "-"+str(x)
		else:
			print "    * WARNING:", natype, "is no valid nucleic acid option, set to default BDNA"
			option = -4
	except:
		if natype == 'ADNA': option = -1
		elif natype == 'BDNA': option = -4
		elif natype == 'CDNA': option = -47
		else:	
			print "    * WARNING:", natype, "is no valid nucleic acid option, set to default BDNA"
			option = -4
		
	cmd = "fiber " + str(option)+" " + (name)+".pdb"" << _End_ \n 2 \n %s \n %d \n _End_\n" % (sequence, repeat)
	os.system(cmd)

def RebuildNA(option, inputfile, outputfile):

	"""Make nucleic acid structure using the "rebuild" module of 3DNA"""	
	
	cmd = "rebuild " + str(option)+" " + (inputfile)+" " + (outputfile) 
	output = commands.getoutput(cmd)
	
def use_Custom(curdir):
	
	"""Use custom supplied atomic models for the bases "A,T,C,G,U". Must be present in /DART/third-party/X3DNA/CUSTOM/"""
	
	customdir = base+'/third-party/X3DNA/CUSTOM/'
	pdbs = ['A.pdb','T.pdb','C.pdb','G.pdb','U.pdb']
	
	if os.path.isdir(customdir): os.chdir(customdir)
	else: print "    * CUSTOM directory is not present. Custom option cannot be used"
	
	tmppdbs = []
	for n in pdbs:
		if os.path.isfile(n):
			print "    * Found pdb file for nucleic acid base", n, "including"
			tmppdbs.append(n)
		else:
			print "    * WARNING: pdb file", n, "not present in CUSTOM directory"
	
	for n in tmppdbs:
		print "    * Use 3DNA std_base function to correct", n, "to match standard reference frame"
		cmd = "std_base " + n + " " + "Atomic_"+n
		os.system(cmd)
	
	correct = glob.glob('Atomic_*.pdb')
	for n in correct:
		print "    * Copy corrected pdb", n, "to working directory", curdir
		shutil.copy(n,curdir)
	
	os.chdir(curdir)
	
def get_Atomic(natype):
	
	"""Copy Atomic_*.pdb files to current directory for rebuilding allatom models"""

	output = commands.getoutput("cp_std %s > /dev/null" % natype)
	
def CleanUp():
	
	"""Cleanup temporary files"""
	
	trash = ['Atomic_A.pdb','Atomic_a.pdb','Atomic_G.pdb','Atomic_g.pdb','Atomic_C.pdb','Atomic_c.pdb','Atomic_T.pdb',
		     'Atomic_t.pdb','Atomic_i.pdb','Atomic_I.pdb','Atomic_U.pdb','Atomic_u.pdb','ref_frames.dat']		
	
	for files in trash:
		if os.path.isfile(files): os.remove(files)

if __name__ == '__main__':

	"""Running from the command line"""
	from optparse import *
	
	"""Parse command line arguments"""
	option_dict = CommandlineOptionParser().option_dict
	inputlist = option_dict['input']
	
	"""Envoce main functions"""
	PluginCore(option_dict,inputlist)
	sys.exit(0)
	
