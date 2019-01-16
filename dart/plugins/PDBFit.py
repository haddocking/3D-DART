#!/usr/bin/env python2.7

USAGE = """
====================================================================================================

Author:				Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center for 
					Biomolecular Research, Utrecht university, The Netherlands.
Copyright (C):		2006 (DART project)
DART version:		1.2 (25-11-2008)
DART plugin: 		PDBFit.py
Plugin excecution:	Either command line driven (use -h/--help for the option) or as part of a 
					DART batch sequence.
Plugin function:	This plugin is a basic wrapper around the functionality of the least-square
					fitting program 'Profit'. 
Dependencies:		PROFIT in ../DART/third-party/profit/

====================================================================================================
"""

"""Import modules"""
import os, sys, commands
from time import ctime

"""Setting pythonpath variables if run from the command line"""
base, dirs = os.path.split(os.path.dirname(os.path.join(os.getcwd(), __file__)))

if base in sys.path:
	pass
else:
	sys.path.append(base)
	
def PluginXML():
	PluginXML = """ 
<metadata>
 <name>Profit least-square fitting wrapper</name>	
 <function>This plugin is a basic wrapper around the functionality of the least-square
  fitting program 'Profit'.</function>		       
 <input type="filetype">.pdb</input>
 <output type="filetype">.list</output>
</metadata>
<parameters>
 <option type="useplugin" form="hidden" text="None">True</option>
 <option type="inputfrom">1</option>
 <option type="reference"></option>
 <option type="atoms"></option>
 <option type="zone"></option>
 <option type="writefit">False</option>
 <option type="default">True</option>
</parameters>"""
	
	return PluginXML

def PluginGui():
	pass
		
def PluginCore(paramdict, metadict, inputlist):
	
	"""Some definitions"""
	
	DNABASE = "^P,O1P,O2P,O5',O4',O3',C5',C4',C3',C2',C1'"
	DNABB = "P,O1P,O2P,O5',O4',O3',C5',C4',C3',C2',C1'"
	PROTBB = "N,CA,C,O"
	PROTSIDE = "^N,CA,C,O"
	ALLHEAVY = "P,N*,C*,O*"
	
	if paramdict['default'] == 'True' or paramdict['default'] == True:
		print "--> Performig default set of protein-DNA fittings"
		print "    * Calculating rmsd full"
		rmsdfull = ProFitting(input1=paramdict['reference'], input2=inputlist, atoms=ALLHEAVY, writefit=paramdict['writefit'], metadict=metadict).rmsd
		print "    * Calculating rmsd dna full"
		rmsddnaall = ProFitting(input1=paramdict['reference'], input2=inputlist, atoms=ALLHEAVY, zone='B*', writefit=paramdict['writefit'], metadict=metadict).rmsd
		print "    * Calculating rmsd protein all"
		rmsdprotall = ProFitting(input1=paramdict['reference'], input2=inputlist, atoms=ALLHEAVY, zone='A*', writefit=paramdict['writefit'], metadict=metadict).rmsd
		print "    * Calculating rmsd dna backbone"
		rmsddnabb = ProFitting(input1=paramdict['reference'], input2=inputlist, atoms=DNABB, zone='B*', writefit=paramdict['writefit'], metadict=metadict).rmsd
		print "    * Calculating rmsd protein backbone"
		rmsdprotbb = ProFitting(input1=paramdict['reference'], input2=inputlist, atoms=PROTBB, zone='A*', writefit=paramdict['writefit'], metadict=metadict).rmsd
		print "    * Calculating rmsd dna base-pairs"
		rmsddnabase = ProFitting(input1=paramdict['reference'], input2=inputlist, atoms=DNABASE, zone='B*', writefit=paramdict['writefit'], metadict=metadict).rmsd
		print "    * Calculating rmsd protein side-chains"   
		rmsdprotside = ProFitting(input1=paramdict['reference'], input2=inputlist, atoms=PROTSIDE, zone='A*', writefit=paramdict['writefit'], metadict=metadict).rmsd
		
		print "    * Write output to rmsd.stat"
		WriteOutput(inputlist,rmsdfull,rmsddnaall,rmsdprotall,rmsddnabb,rmsdprotbb,rmsddnabase,rmsdprotside)
	
#================================================================================================================================#
# 					PLUGIN SPECIFIC DEFINITIONS BELOW THIS LINE						 #
#================================================================================================================================#

def WriteOutput(inputlist,rmsdfull,rmsddnaall,rmsdprotall,rmsddnabb,rmsdprotbb,rmsddnabase,rmsdprotside):
	
	"""Write RMSD values to rmsd.stat"""
	
	outfile = file('rmsd.stat','w')
		
	outfile.write('*************************************************************************************************************************\n')	               
	outfile.write('Root means square fitting data for %i structures\n' % (len(inputlist))) 		               
	outfile.write('Time and date: %s\n' % ctime())	
	outfile.write('Legende:h.a = heavy atom, bb = backbone atoms\n')										               
	outfile.write('When there is no rmsd value given but --- instead, something has gone wrong during the fitting\n')	               
	outfile.write('*************************************************************************************************************************\n')	               									               
	outfile.write('structure        full(h.a.)   dna(h.a.)   prot(h.a.)    dna(b.b.)     prot(b.b.)     dna(base)      prot(side)\n')
	
	i = 0
	while i < len(inputlist):
		outfile.write('%2s%8s%8s%8s%8s%8s%8s%8s\n' % 
		(os.path.split(inputlist[i])[-1],rmsdfull[inputlist[i]],rmsddnaall[inputlist[i]],rmsdprotall[inputlist[i]],
		 rmsddnabb[inputlist[i]],rmsdprotbb[inputlist[i]],rmsddnabase[inputlist[i]],rmsdprotside[inputlist[i]]))
		i = i+1
		
	outfile.close()

class ProFitting:
	
	"""Full features wrapper around ProFit"""
	
	def __init__(self, input1=None, input2=None, atoms=None, zone=None, writefit=False, metadict=None):
		
		self.profit = metadict['dependencies']+"/profit"
		self.input1 = input1
		self.input2 = input2
		self.atoms = atoms
		self.zone = zone
		self.writefit = writefit
		
		self.rmsd = {}
		self.RunProfit()
		
	def _ConstructOptionString(self, files):
		
		if not self.writefit == False:
			write = '\nwrite '+os.path.splitext(os.path.split(files)[-1])[0]+'_fit.pdb'
		else:
			write = ''
		
		if not self.zone == None:
			zone = '\nzone '+self.zone
		else:
			zone = ''
			
		atoms = '\natoms '+self.atoms
		quit = '\nquit'
		fit = '\nfit '
		
		profit_cmd = self.profit + " " + self.input1 + " " + files
		options_cmd = 'echo "'+ atoms + zone + fit + write + quit + '" | ' + profit_cmd + '| grep RMS'
		
		return options_cmd
		
	def _FormatOutput(self, rmsd):
		
		try:
			l = (rmsd.strip()).split(' ')
			return (float(l[1]))
		except:
			return ('----')
		
	def RunProfit(self):
		
		for files in self.input2:
			rmsd = commands.getoutput(self._ConstructOptionString(files))	
			self.rmsd[files] = self._FormatOutput(rmsd)
		
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
		parser.add_option( "-r", "--reference", action="store", dest="reference", type='string', help="Define reference structure for fitting")
		parser.add_option( "-a", "--atoms", action="store", dest="atoms", help="Define atoms to fit on (ProFit syntax)")
		parser.add_option( "-z", "--zone", action="store", dest="zone", help="Define zone to fit on (ProFit syntax")
		parser.add_option( "-w", "--writefit", action="store_true", dest="writefit", default=False, help="Output a pdb file of the fitted structure")
		parser.add_option( "-d", "--default", action="store_true", dest="default", default=False, help="Perform default HADDOCK fittings")
		
		(options, args) = parser.parse_args()
		
		self.option_dict['input'] = options.inputfile
		self.option_dict['reference'] = options.reference
		self.option_dict['atoms'] = options.atoms
		self.option_dict['zone'] = options.zone
		self.option_dict['writefit'] = options.writefit
		self.option_dict['default'] = options.default
			
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
		
if __name__ == '__main__':

	"""Running this plugin from the command line"""
	
	"""Import command line specific modules"""
	from optparse import *
	
	"""Setting up parameter dictionary"""
	metadict = {}
	metadict['dependencies'] = base+"/third-party/profit"
	paramdict = {'showonexec':'False','inputfrom':'self','autogenerateGui':'False'}
	option_dict = CommandlineOptionParser().option_dict
	
	for key in option_dict:
		paramdict[key] = option_dict[key]
	del option_dict
	
	"""Check for input"""
	if paramdict['input'] == None:
		print "    * Please supply pdb file using option -f or use option -h/--help for usage"
		sys.exit(0)
	else:
		inputlist = paramdict['input']
	
	"""Running plugin main functions"""
	PluginCore(paramdict, metadict, inputlist)
	
	"""Say goodbye"""
	print "--> Thanks for using PDBFit, bye"
	sys.exit(0)
