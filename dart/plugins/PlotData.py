#!/usr/bin/env python2.7

#Plotting base-pair and base-pair step data fro  multistat files

USAGE = """
==========================================================================================

Author:				Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center
					for Biomolecular Research, Utrecht university, The Netherlands
Copyright (C):		2006 (DART project)
DART version:		1.2  (25-11-2008)
DART plugin: 		PlotData.py
Input:		   
Output:				Gnuplot graphs in postscript file format or to screen.
Plugin excecution:	Either command line driven (use -h/--help for the option) or as 
					part of a DART batch sequence.	
Plugin function: 	This utility scripts allows for the easy generation of graphs from Nucleic Acid
			        analysis data. It uses Gnuplot for the graph generation and the Gnuplot module
			        for python to call Gnuplot from python.
Examples:		    PlotData -vf multibend.stat
					PlotData -f multiout.stat -r 3CRO.par -n twist
Dependencies:		Standard python installation, Gnuplot, Python Gnuplot bindings

==========================================================================================
"""

"""Import modules"""
import Gnuplot, re, os, sys

"""Setting pythonpath variables if run from the command line"""
base, dirs = os.path.split(os.path.dirname(os.path.join(os.getcwd(), __file__)))

if base in sys.path:
	pass
else:
	sys.path.append(base)

def PluginCore(paramdict, inputlist):
	
	plot = PlotData()
	
	for files in inputlist:
		base = os.path.basename(files)
		if base == 'multiout.stat':
			plot.ReadMultiout(files)
			if not paramdict['reference'] == None and os.path.splitext(paramdict['reference'])[1] == '.par':
				plot.ReadPar(paramdict['reference'])
			plot.PlotParamData(paramdict['verbose'],paramdict['name'])
		elif base == 'multibend.stat':
			plot.ReadMultibend(files)
			if not paramdict['reference'] == None and os.path.splitext(paramdict['reference'])[1] == '.bend':
				plot.ReadBend(paramdict['reference'])
			plot.PlotBendData(paramdict['verbose'],paramdict['name'])
		elif os.path.splitext(files)[1] == '.par':
			plot.ReadPar(files)
			plot.PlotParamData(paramdict['verbose'],paramdict['name'])
		elif os.path.splitext(files)[1] == '.bend':
			plot.ReadBend(files)		
			plot.PlotBendData(paramdict['verbose'],paramdict['name'])
		else:
			print("The File: %s is not a allowd input file" % base)	

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
		
		parser.add_option( "-f", "--file", action="callback", callback=self.varargs, dest="inputfile", type="string", help="Supply inputfile(s)")
		parser.add_option( "-n", "--name", action="store", dest="name", help="Only generate a graph for a specific parameter")
		parser.add_option( "-r", "--reference", action="store", dest="reference", help="Include data from a reference structure to the graph")
		parser.add_option( "-v", "--verbose", action="store_true", dest="verbose", default=True, help="Print the graph to screen instead of a file")
	
		(options, args) = parser.parse_args()
		
		self.option_dict['input'] = options.inputfile
		self.option_dict['name'] = options.name
		self.option_dict['reference'] = options.reference
		self.option_dict['verbose'] = options.verbose
			
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

class PlotData:

	"""Read data from multiout.stat, multibend.stat or individual base-pair(step) parameter files
	   and present the data in a wel formated graph using Gnuplot."""

	def __init__(self):

		self.sequence = None
		self.bpparm = {}
		self.bpparmsd = {}
		self.bpparmr = {}
		self.bend = {}
		self.bendsd = {}
		self.bendref = {}

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
			elif outfile[count][10] == "~":
				break
			else:	
				tablelines.append(outfile[count].split())
			count = count+1

		return tablelines

	def TransformDash(self,n):

		"""Transform dashes often found in 3DNA tables to floats"""
		
		if n == '---' or n == '----':
			return float(0)
		else:
			return float(n)

	def ReadPar(self,infile):

		"""Read 3DNA generated base-pair and base-pair step parameter file (*.par) and append
		   to bpparr database"""

		readfile = file(infile, 'r')
		lines = readfile.readlines()

		lib = {}
		for n in range(13):
			lib[n] = []

		for line in lines[3:]:
			params = []
			for n in range(0,13):	
		 		params.append(line.split()[n])	
			lib[0].append(str(params[0]))
			for value in range(1,13):
				lib[value].append(float(params[value]))

		if self.sequence == None:
			self.sequence = lib[0]
			self.bpparmr['shear'] = lib[1]
			self.bpparmr['stretch'] = lib[2]
			self.bpparmr['stagger'] = lib[3]
			self.bpparmr['buckle'] = lib[4]
			self.bpparmr['proptw'] = lib[5]
			self.bpparmr['opening'] = lib[6]
			self.bpparmr['shift'] = lib[7]
			self.bpparmr['slide'] = lib[8]
			self.bpparmr['rise'] = lib[9]
			self.bpparmr['tilt'] = lib[10]
			self.bpparmr['roll'] = lib[11]
			self.bpparmr['twist'] = lib[12]
		elif len(self.sequence) == len(lib[0]):
			self.bpparmr['shear'] = lib[1]
			self.bpparmr['stretch'] = lib[2]
			self.bpparmr['stagger'] = lib[3]
			self.bpparmr['buckle'] = lib[4]
			self.bpparmr['proptw'] = lib[5]
			self.bpparmr['opening'] = lib[6]
			self.bpparmr['shift'] = lib[7]
			self.bpparmr['slide'] = lib[8]
			self.bpparmr['rise'] = lib[9]
			self.bpparmr['tilt'] = lib[10]
			self.bpparmr['roll'] = lib[11]
			self.bpparmr['twist'] = lib[12]
		else:
			print("Number of base-pairs in reference does not match the number of base-pairs in the model(s): %i versus %i" % (len(lib[0]),len(self.sequence)))
			print("\nNo reference used")	
		
		lib.clear()

	def ReadBend(self,infile):

		"""Read nucleic acid bend analysis file (*.bend) and append to bendref database"""

		readfile = file(infile, 'r')
		lines = readfile.readlines()
		
		lib = {}
		for n in range(9):
			lib[n] = []

		for line in lines[18:]:
			params = line.split()
			try:
				lib[0].append(str(params[1]))
				for value in range(2,9):
					lib[value].append(float(params[value]))
			except:
				break		

		if self.sequence == None:
			self.sequence = lib[0]
			self.bendref['fractilt'] = lib[2]
			self.bendref['fracroll'] = lib[3]
			self.bendref['bpangle'] = lib[4]
			self.bendref['orient'] = lib[5]
			self.bendref['gtilt'] = lib[6]
			self.bendref['groll'] = lib[7]
			self.bendref['acctwist'] = lib[8]
		elif len(self.sequence) == len(lib[0]):
			self.bendref['fractilt'] = lib[2]
			self.bendref['fracroll'] = lib[3]
			self.bendref['bpangle'] = lib[4]
			self.bendref['orient'] = lib[5]
			self.bendref['gtilt'] = lib[6]
			self.bendref['groll'] = lib[7]
			self.bendref['acctwist'] = lib[8]
		else:
			print("Number of base-pairs in reference does not match the number of base-pairs in the model(s): %i versus %i" % (len(lib[0]),len(self.sequence)))
			print("\nNo reference used")
				
		lib.clear()
		
	def ReadMultiout(self, infile):
	
		"""Read statistics for base-pair and base-pair step parameters from multiout file and 
		   append to database"""
	   
		lib = {}
		for n in range(26):
			lib[n] = []
	
		bp = re.compile("nr   bp      freq.        shear           stretch")
		bpstep = re.compile("nr   bpstep    freq.      shift           slide")
		readfile = file(infile, 'r')
		lines = readfile.readlines()
		linecount = 1
		countstart = []
		for line in lines:
			line = line.strip()
			result1 = bp.match(line)
			result2 = bpstep.match(line)
			if result1:
				bplines = self._ReadTable(linenr=len(countstart),outfile=lines)	
			elif result2:
				bpsteplines = self._ReadTable(linenr=len(countstart),outfile=lines)
			linecount += 1
	    		countstart.append(linecount)	

		bpsteplines.insert(0,([0.00001]*15))

		for lines in bplines:
			lines = lines[1:15]
			for value in range(14):
				try:
					lib[value].append(self.TransformDash(lines[value]))
				except:
					lib[value].append(lines[value])
		for lines in bpsteplines:
			lines = lines[3:15]
			for value in range(12):
				lib[value+14].append(self.TransformDash(lines[value]))	
	
		del bpsteplines,bplines
	
		self.sequence = lib[0]
		self.bpparm['shear'] = lib[2]
		self.bpparmsd['shearsd'] = lib[3]
		self.bpparm['stretch'] = lib[4]
		self.bpparmsd['stretchsd'] = lib[5]
		self.bpparm['stagger'] = lib[6]
		self.bpparmsd['staggersd'] = lib[7]
		self.bpparm['buckle'] = lib[8]
		self.bpparmsd['bucklesd'] = lib[9]
		self.bpparm['proptw'] = lib[10]
		self.bpparmsd['proptwsd'] = lib[11]		
		self.bpparm['opening'] = lib[12]
		self.bpparmsd['openingsd'] = lib[13]
		self.bpparm['shift'] = lib[14]
		self.bpparmsd['shiftsd'] = lib[15]
		self.bpparm['slide'] = lib[16]
		self.bpparmsd['slidesd'] = lib[17]
		self.bpparm['rise'] = lib[18]
		self.bpparmsd['risesd'] = lib[19]
		self.bpparm['tilt'] = lib[20]
		self.bpparmsd['tiltsd'] = lib[21]
		self.bpparm['roll'] = lib[22]
		self.bpparmsd['rollsd'] = lib[23]		
		self.bpparm['twist'] = lib[24]
		self.bpparmsd['twistsd'] = lib[25]
	
		lib.clear()

	def ReadMultibend(self,infile):

		"""Read multibend file and append parameters to database"""

		ref = re.compile("Reference base-pair:")
		start = re.compile("index  bp-step")
		readfile = file(infile, 'r')
		lines = readfile.readlines()
		linecount = 1
		countstart = []
		for line in lines:
			line = line.strip()
			result1 = ref.match(line)
			result2 = start.match(line)
			if result1:
				refbp = float(line.split()[2])	
			elif result2:
				bendlines = self._ReadTable(linenr=len(countstart),outfile=lines)
			linecount += 1
        		countstart.append(linecount)	

		lib = {}
		for n in range(17):
			lib[n] = []
		for lines in bendlines:
			for value in range(17):
				try:
					lib[value].append(self.TransformDash(lines[value]))
				except:	
					lib[value].append(lines[value])

		self.sequence = lib[1]
		self.refbp = refbp
		self.bend['fractilt' ] = lib[3]
		self.bendsd['fractiltsd'] = lib[4]
		self.bend['fracroll'] = lib[5]
		self.bendsd['fracrollsd'] = lib[6]
		self.bend['bpangle'] = lib[7]
		self.bendsd['bpanglesd'] = lib[8]
		self.bend['orient'] = lib[9]
		self.bendsd['orientsd'] = lib[10]
		self.bend['gtilt'] = lib[11]
		self.bendsd['gtiltsd'] = lib[12]
		self.bend['groll'] = lib[13]
		self.bendsd['grollsd'] = lib[14]		
		self.bend['acctwist'] = lib[15]
		self.bendsd['acctwistsd'] = lib[16]

		del bendlines
		lib.clear()

	def PlotParamData(self,verbose,param=None):	

		if not param == None:
			paramlist = [param]
		else:
			paramlist = self.bpparm.keys()	
	
		for param in paramlist:
			index = range(1,len(self.sequence)+1)
			xtics = ""
			for n in range(len(self.sequence)):
				if n+1 == len(self.sequence):
					xtics = xtics+'"'+self.sequence[n]+'" '+str(n+1)
				else:	
					xtics = xtics+'"'+self.sequence[n]+'" '+str(n+1)+','
		
			command = "g.plot(\\"
			if param in self.bpparm:
				command = command+"\nGnuplot.Data(index, self.bpparm[param], with='boxes')"	
			if (param+"sd") in self.bpparmsd:
				command = command+","
				command = command+"\nGnuplot.Data(index, self.bpparm[param], self.bpparmsd[param+'sd'], with='errorbars lt 3')"
			if param in self.bpparmr and param in self.bpparm:
				command = command+","
				command = command+"\nGnuplot.Data(index, self.bpparmr[param], with='points 3')"
			elif param in self.bpparmr:
				command = command+"\nGnuplot.Data(index, self.bpparmr[param], with='boxes')"
			else:
				pass
			command = command+"\n)"	
			
			if command == "g.plot(\\\n)":
				print("The requested data to plot was not found: %s" % param)
				sys.exit(0)
		
			g = Gnuplot.Gnuplot()
			g.title('Statistical data for base-pair(step) parameter: %s' % param)
			g.xlabel('Base-pair')
			g.ylabel(param)
			g('set grid')
			g('set xtics (%s)' % xtics)
			g('set boxwidth 0.5')
			g('set style fill solid 1.00 border -1')
			
			exec(command)
			
			if verbose == False:
				g.hardcopy(param+'.ps', enhanced=1, color=1)
	
	def PlotBendData(self,verbose,param=None):
	
		if not param == None:
			paramlist = [param]
		else:
			paramlist = self.bend.keys()	
	
		for param in paramlist:
			index = range(1,len(self.sequence)+1)
			xtics = ""
			for n in range(len(self.sequence)):
				if n+1 == len(self.sequence):
					xtics = xtics+'"'+self.sequence[n]+'" '+str(n+1)
				else:	
					xtics = xtics+'"'+self.sequence[n]+'" '+str(n+1)+','
			
			command = "g.plot(\\"
			if param in self.bend:
				command = command+"\nGnuplot.Data(index, self.bend[param], with='boxes')"	
			if (param+"sd") in self.bendsd:
				command = command+","
				command = command+"\nGnuplot.Data(index, self.bend[param], self.bendsd[param+'sd'], with='errorbars lt 3')"
			if param in self.bendref and param in self.bend:
				command = command+","
				command = command+"\nGnuplot.Data(index, self.bendref[param], with='points 3')"
			elif param in self.bendref:
				command = command+"\nGnuplot.Data(index, self.bendref[param], with='boxes')"
			else:
				pass
			command = command+"\n)"	
		
			if command == "g.plot(\\\n)":
				print("The requested data to plot was not found: %s" % param)
				sys.exit(0)
		
			g = Gnuplot.Gnuplot()
			g.title('Statistical data for nucleic acid bend parameter: %s' % param)
			g.xlabel('Base-pair')
			g.ylabel(param)
			g('set grid')
			g('set xtics (%s)' % xtics)
			g('set boxwidth 0.5')
			g('set style fill solid 1.00 border -1')
			
			exec(command)
			
			if verbose == False:
				g.hardcopy(param+'.ps', enhanced=1, color=1)
	
if __name__ == '__main__':
	
	"""Running from the command line"""
	from optparse import *
	option_dict = option_dict = CommandlineOptionParser().option_dict
	inputlist = option_dict['input']
    
	"""Envoce main functions"""
	PluginCore(option_dict, inputlist)
	sys.exit(0)
	
	
