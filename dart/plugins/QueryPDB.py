#!/usr/bin/env python2.7

USAGE = """
==========================================================================================

Author:				Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center
					for Biomolecular Research, Utrecht university, The Netherlands
Copyright (C):		2006 (DART project)
DART version:		1.2 (01-01-2007)
DART plugin: 		QueryPDB.py
Plugin function:	Makes a hirargical representation of the PDB (same as QueryPDB) 
					and and provides various function to perform calculations on the
					data in the PDB
Dependencies:		None

==========================================================================================
"""

"""Import modules"""
import os, sys, re

"""Setting pythonpath variables if run from the command line"""
base, dirs = os.path.split(os.path.dirname(os.path.join(os.getcwd(), __file__)))

if base in sys.path:
	pass
else:
	sys.path.append(base)

from PDBeditor import PDBeditor
from system.Xpath import Xpath
from system.NAfunctionLib import CalculateDistance
from system.Constants import *
from numpy import *

def PluginXML():
	PluginXML = """ 
<metadata>
 <name>QueryPDB</name>
 <function>Makes a hirargical representation of the PDB (same as QueryPDB) 
  and and provides various function to perform calculations on the data in 
  the PDB</function>
 <input type="Filetype">.pdb</input>
 <output type="Filetype">.txt</output>
</metadata>
<parameters>
 <option type="useplugin" form="hidden" text="None">True</option>
 <option type="inputfrom" form="hidden" text="None">1</option>
 <option type="contact" form="checkbox" text="Calcualte contacts">False</option>
 <option type="cutoff" form="text" text="Contacts calculation upper distance cutoff">5</option>
</parameters>"""
	
	return PluginXML

def PluginGui():
	pass

def PluginCore(paramdict, inputlist):

	for files in inputlist:
		
		"""Converting PDB file to XML in memory"""
		if (os.path.splitext(files))[1] == '.xml':
			xml = files
		else:
			pdb = PDBeditor()
			pdb.ReadPDB(files)	
			xml = pdb.PDB2XML().xml()
		
		if paramdict['sequence'] == True:
		
			sequence = GetSequence()
			sequence.GetSequence(pdbxml=xml)
			sequence.FormatOutput()
		
		if paramdict['NAsummery'] == True:
		
			print("--> Starting nucleic-acid structure evaluation process on structure %s" % os.path.basename(files))
			print "    * Getting sequence information"
	
			sequence = GetSequence()
			sequence.GetSequence(pdbxml=xml)
			
			naeval = NAsummery(pdbxml=xml,sequence=sequence.seqlib)
			naeval.Evaluate()
				
#================================================================================================================================#
# 					PLUGIN SPECIFIC DEFINITIONS BELOW THIS LINE						 #
#================================================================================================================================#

"""Some defaults"""
DNATHREE = ['GUA','CYT','THY','ADE'] 
DNAONE = ['G','C','T','A']
RNATHREE = ['URI']
RNAONE = ['U']

class CommandlineOptionParser:
	
	"""Parses command line arguments using optparse"""
	
	def __init__(self):
		
		self.option_dict = {}
		self.option_dict = self.CommandlineOptionParser()
	
	def CommandlineOptionParser(self):
	
		"""Parsing command line arguments"""
	
		usage = "usage: %prog [options] arg"
		parser = OptionParser(usage)

		parser.add_option( "-f", "--file", action="callback", callback=self.varargs, dest="inputfile", type="string", help="Supply pdb inputfile(s)")
		parser.add_option( "-n", "--NAsummery", action="store_true", dest="NAsummery", default=False, help="Print summery of nucleic acid structure data")
		parser.add_option( "-s", "--sequence", action="store_true", dest="sequence", default=False, help="Return the sequence of all chains in the PDB")
		
		(options, args) = parser.parse_args()
		
		self.option_dict['input'] = options.inputfile
		self.option_dict['NAsummery'] = options.NAsummery
		self.option_dict['sequence'] = options.sequence
			
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

class GetSequence:
			
	"""Return the sequence of the supplied PDB"""
	
	def __init__(self):
	
		self.seqlib = {}
	
	def GetSequence(self,pdbxml=None):
	 			
		"""Append sequence and resid nr of all chains in pdb to seqlib"""
		query = Xpath(pdbxml)
		
		query.Evaluate(query={1:{'element':'chain','attr':None}})
		for chain in query.nodeselection[1]:
			query.getAttr(node=chain,selection='ID',export='string')
		chains = query.result
		query.ClearResult()
		
		for chain in chains:
			query.Evaluate(query={1:{'element':'chain','attr':{'ID':chain}},2:{'element':'resid','attr':None}})
			for resid in query.nodeselection[2]:
				query.getAttr(node=resid,selection=['ID','nr'])
			residues = query.result 
			resid = []
			resnr = []
			for residue in residues:
				resid.append(residue[0])
				resnr.append(residue[1])
			self.seqlib[chain] = []	
			self.seqlib[chain].append(resid)
			self.seqlib[chain].append(resnr)
			query.ClearResult()	

	def FormatOutput(self):
	
		for chain in self.seqlib.keys():
			print("Chain: %s" % chain)
			print("start at resid nr %i end at resid nr %i" % (min(self.seqlib[chain][0]),max(self.seqlib[chain][0])))
			print("sequence:")
			for resid in self.seqlib[chain][1]:
				print("%s " % resid), 

class NAsummery:

	"""Evaluate the structure of a nucleic acid on: type, chains and pairing"""

	def __init__(self,pdbxml=None,sequence=None):
		
		self.pdbxml = pdbxml
		self.sequence = sequence
		
		self.moltype = {}
		self.chainlib = {}
		self.pairs = {}
		self.cutoff = []
		
	def Evaluate(self):
	 	
		"""First indentify type of chain in sequence"""
	 	chains = self.sequence.keys()
		if len(chains) == 0:
			print "    * ERROR: No chains found in structure, stopping"
			sys.exit(0)
		else:
		 	print("    * Found %i chains in structure" % len(chains))
			self._EvalMolType(chains)
		
		"""Find segments (if any) in chains"""
		for chain in self.moltype.keys():
			if self.moltype[chain] == 'DNA' or self.moltype[chain] == 'RNA':
				self._BackboneTrace(chain)
		
		"""Extract paired chains/segments"""
		for chain in self.chainlib.keys():
			self._FindPairs(chain)
	
		"""TEMPORARY HACK TO MAKE X3DNAANALYZE WORK"""
		if len(self.chainlib) == 1:
			for keys in self.chainlib.keys():
				self.chainlib[keys][1].reverse()
				self.pairs[keys][1].reverse()
		elif len(self.chainlib) == 2:
			chain = {}
			chain[self.chainlib.keys()[0]] = []
			for keys in self.chainlib.keys():
				chain[self.chainlib.keys()[0]].append(self.chainlib[keys][0])
			chain[self.chainlib.keys()[0]][1].reverse()
			pair = {}
			pair[self.chainlib.keys()[0]] = []
			for keys in self.pairs.keys():
				pair[self.chainlib.keys()[0]].append(self.pairs[keys][0])
			pair[self.chainlib.keys()[0]][1].reverse()
			self.chainlib = chain
			self.pairs = pair
			
	def _EvalMolType(self,chains):
		
		print "    * Evaluate chain molecule type (RNA, DNA, protein)"
		for chain in chains:
			for resid in self.sequence[chain][1]:
				if resid in RNATHREE:
					self.moltype[chain] = 'RNA'
					break				#First occurence of URI, claim RNA. Little messy
				elif resid in RNAONE:
					self.moltype[chain] = 'RNA'
					break				#First occurence of U, claim RNA. Little messy
				elif resid in DNATHREE:					
					self.moltype[chain] = 'DNA'
				elif resid in DNAONE:					
					self.moltype[chain] = 'DNA'	
				elif resid in PROTTHREE:
				   	self.moltype[chain] = 'PROT'
				elif resid in PROTONE:
				   	self.moltype[chain] = 'PROT'	
				else:
					pass		
		
		for chain in self.moltype.keys():
			print( "      Chain %s is indentified as moltype %s" % (chain,self.moltype[chain]))
	
	def _BackboneTrace(self,chainid):
				
		"""Calculate same-strand C5' to C5' distance"""		
		
		print("    * Indentify segments for chain %s" % chainid)
		print "    * Calculating same-strand C5' to C5' distance to extract segments from structure. Segment indentified"
		print "      when C5'-C5' distance is larger than dynamic average + standard deviation + 1 = cutoff"
		print "      Distance  Residue  Residue+1  Cutoff"
		
		query = Xpath(self.pdbxml)
		query.Evaluate(query={1:{'element':'chain','attr':{'ID':chainid}},2:{'element':'resid','attr':None},3:{'element':'atom','attr':{'ID':"C5'"}}})
		
		for resid in query.nodeselection[2]:
		 	query.getAttr(node=resid,selection='nr',export='string')
		residues = query.result
		query.ClearResult()
		
		for atom in query.nodeselection[3]:
			query.getAttr(node=atom,selection=['corx','cory','corz'])
		atoms = query.result
		query.ClearResult()
		
		chain = []
		self.chainlib[chainid] = []
		for residue in range(len(residues)):
			try:
				distance = CalculateDistance(atoms[residue],atoms[residue+1])
				cutoff = self._CalcCutoff(distance)
				print("      %1.4f %6i %6i %15.4f" % (distance,(residues[residue]),(residues[residue+1]),cutoff)) 
				if distance < cutoff:
					chain.append(residues[residue])
				else:
					chain.append(residues[residue])
					self.chainlib[chainid].append(chain)
					self.cutoff = []
					chain = []
			except:
				chain.append(residues[residue])
				self.chainlib[chainid].append(chain)
			
		if len(self.chainlib[chainid]) == 0:
			print "    * ERROR: no segments indentified, stopping"
			sys.exit(0)
		else:
			print("    * Identified %i segment(s):" % len(self.chainlib[chainid]))	
			for chain in range(len(self.chainlib[chainid])):
				print("      Segment %i range: %i to %i" % (chain+1,min(self.chainlib[chainid][chain]),max(self.chainlib[chainid][chain])))		
		
	def _CalcCutoff(self,distance):
		
		self.cutoff.append(distance)
	
		if len(self.cutoff) == 1:
			cutoff = self.cutoff[0]+1
		else:	
			average = mean(self.cutoff)
			stdev = std(self.cutoff)
			cutoff = average+stdev+1
		
		return cutoff			
	
	def _FindPairs(self, chain):
		
		self.pairs[chain] = []
		
		if len(self.chainlib[chain]) > 1:
			minl = 0
			maxl = 0
			for segment in self.chainlib[chain]:
				maxl = maxl + len(segment)
				self.pairs[chain].append(self.sequence[chain][1][minl:maxl])
				minl = maxl
				maxl += 1	
		else:
			self.pairs[chain].append(self.sequence[chain][1])		

if __name__ == '__main__':

	"""For testing purposes"""
	from optparse import *
	
	"""Parse command line arguments"""
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
	
	"""Envoce main functions"""
	PluginCore(paramdict, inputlist)

	sys.exit(0)
