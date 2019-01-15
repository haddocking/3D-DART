#!/usr/bin/env python2.7

USAGE="""
==========================================================================================

Author:			Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center for 
			Biomolecular Research, Utrecht university, The Netherlands.
Copyright (C):		2006 (DART project)
DART version:		1.0 (01-01-2007)
DART plugin: 		ModelNucleicAcids.py
Input:			PAR data file and any of the allowed options, any order and 
			combined
Output:			One or more new PAR file(s) and a summery file of the generated 
			PAR files.
Plugin excecution:	Either command line driven (use -h/--help for the option) or as 
			part of a DART batch sequence.
Plugin function:	This plugin allows you to model Nucleic Acids by introduce local 
			and global conformational changes. Modelling is accomplished by 
			manipulation of the 12 base-pair and base-pair step parameters. To 
			use this plugin you need to have a base-pair/base-pair step 
			parameter file. The global conformational changes you can introduce 
			cover nucleic acid bending and twisting and control over the 
			direction of the bend in space.
Plugin dependencies:	None

for further information, please contact:
			- DART website (http://www.nmr.chem.uu.nl/DART)
			- email: abonvin@chem.uu.nl

If you are using this software for academic purposes please quoting the following 
reference:

===========================================================================================
"""

"""Import modules"""
import sys, os, re
from math import sqrt
from time import ctime
from numpy import *
from copy import deepcopy

"""Setting pythonpath variables if run from the command line"""
base, dirs = os.path.split(os.path.dirname(os.path.join(os.getcwd(), __file__)))

if base in sys.path:
	pass
else:
	sys.path.append(base)

"""Import DART specific modules"""
from system.NAfunctionLib import *
from system.IOlib import *
from system.Utils import TransformDash,MakeBackup
from system.Constants import *
from BuildNucleicAcids import FiberModule

def PluginXML():
	PluginXML = """ 
<metadata>
 <name>Nucleic Acids Modeling Routines</name>
 <function>This plugin allows you to model Nucleic Acids by introduce local and global 
  conformational changes. Modelling is accomplished by manipulation of the 12 base-pair 
  and base-pair step parameters. To use this plugin you need to have a base-pair/base-pair 
  step parameter file. The global conformational changes you can introduce cover nucleic 
  acid bending and twisting and control over the direction of the bend in space.</function>		       
 <input type="filetype">.par,multibend.stat,multiout.stat</input>
 <output type="filetype">.par</output>
</metadata>
<parameters>
 <option type="useplugin" form="hidden" text="None">True</option>
 <option type="inputfrom" form="hidden" text="None">1</option>
 <option type="name" form="text" text="Give your structure a name"></option>
 <option type="input" form="file" text="Upload base-pair/base-pair step parameter file"></option>
 <option type="verbose" form="checkbox" text="Verbose output">False</option>
 <option type="automodel" form="checkbox" text="Automatic modeling">True</option>
 <option type="gltolerance" form="text" text="Deviation from ensamble average for global conformation, scale 0-1>">1.0</option>
 <option type="glvariance" form="text" text="Deviation from optimal base-pair(step) parameter for global conformation, scale 0-1>">0.8</option>
 <option type="glsmoothing" form="text" text="Smoothing for twist, roll and twist. scale 0-1">0.8</option>
 <option type="lcvariance" form="text" text="Deviation from optimal base-pair(step) parameter for local conformation, scale 0-1>">0.8</option>
 <option type="number" form="text" text="Number of models to generate">10</option>
 <option type="refbp" form="text" text="Base-pair used as origin in global reference frame"></option>
 <option type="maxangle" form="text" text="Define maximum bend angle"></option>
 <option type="minangle" form="text" text="Define either the minimum bend angle or give a custom sequence of bend angles"></option>
 <option type="anglestep" form="text" text="Define bend-angle step size">1</option>
 <option type="startbp" form="text" text="Define start base-pair for bending"></option>
 <option type="endbp" form="text" text="Define end base-pair for bending"></option>
 <option type="minorient" form="text" text="Define either the minimum bend angle in space or give a custom sequence"></option>
 <option type="maxorient" form="text" text="Define a maximum angle of bend in space"></option>
 <option type="orientstep" form="text" text="Define step size of bend angle in space"></option>
 <option type="bpstep" form="text" text="Supply custom values for base-pair step parameters as a komma seperated list"></option>
 <option type="bp" form="text" text="Supply custom values for base-pair parameters as a komma seperated list"></option>
 <option type="mingroove" form="text" text="Define minimum major groove width"></option>
 <option type="maxgroove" form="text" text="Define minimum minor groove width"></option>
 <option type="groovestep" form="text" text="Define major groove width step size"></option>
 <option type="groovestart" form="text" text="Define start base-pair for major groove width change"></option>
 <option type="grooveend" form="text" text="Define end base-pair for major groove width change"></option>
 <option type="helicalphase" form="text" text="Set the DNA helical phasing (bp/turn)"></option>
</parameters>"""

	return PluginXML

def PluginCore(paramdict, inputlist):
	
	print "--> Initializing modelling process"
	
	base = []
	stats = []
 
	for files in inputlist:
		basename = os.path.basename(files)
		if basename == 'multibend.stat' or basename == 'multiout.stat':
			stats.append(files)
		elif os.path.splitext(files)[1] == '.par' or os.path.splitext(files)[1] == '.out':
		 	base.append(files)
	
	if len(stats) > 0:
		paramdict['stats'] = stats
	else:
		paramdict['stats'] = None

	if len(base) > 0:
	  # Read napairing.stat file to determine basefile without unpaired bases to use are reference file.
	  for files in inputlist:
	    if os.path.basename(files) == 'napairing.stat':
	      basefile = os.path.join(os.path.split(files)[0], readNApairing(files))
	      if '%s.par' % basefile in base:
	        paramdict['base'] = '%s.par' % basefile
	      elif '%s.out' % basefile in base:
	        paramdict['base'] = '%s.out' % basefile
	      else:
	        pass    
	      break
	  if not paramdict.has_key('base'):
	    paramdict['base'] = base[0]
	else:
		paramdict['base'] = None
		
	model = ModelNucleicAcids(paramdict)
					
	if paramdict['base'] == None and not paramdict['stats'] == None:
		print "    * WARNING: No .par file supplied, an attempt will be made to construct one from the sequence in one of the supplied statistics files"
	elif paramdict['base'] == None and paramdict['stats'] == None:
		print "    * ERROR: No input supplied."
		sys.exit(0)
	elif os.path.splitext(os.path.basename(paramdict['base']))[1] == '.par':
		print("    * Start modelprocess on base-pair(step) parameter file %s" % os.path.basename(paramdict['base']))
		model.ReadPar(paramdict['base'])
	elif os.path.splitext(os.path.basename(paramdict['base']))[1] == '.out':
		print("    * Start modelprocess on 3DNA analysis file %s" % os.path.basename(paramdict['base']))
		model.ReadOut(paramdict['base'])
	else:
		print("    * ERROR: the file %s is not an excepted file format, stopping" % os.path.basename(paramdict['base']))
		sys.exit(0)
	
	if not paramdict['stats'] == None:
		for files in paramdict['stats']:
			if os.path.basename(files) == 'multiout.stat':
				print "    * Importing statistics from multiout.stat file"
				model.ReadMultiout(files)
			if os.path.basename(files) == 'multibend.stat':
				print "    * Importing statistics from multibend.stat file"
				model.ReadMultibend(files)

	if paramdict['automodel'] == True and not paramdict['stats'] == None:
		print "    * Initializing automodel process"
		model.Automodel()
	else:
		print "    * Initializing manual modeling process"
		model.CheckParamdict()
		model.Manualmodel()	
				
#================================================================================================================================#
# 					PLUGIN SPECIFIC DEFINITIONS BELOW THIS LINE						 											 #
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

		parser.add_option( "-f", "--file", action="callback", callback=self.varargs, dest="input", type="string", help="Supply base parameter file and/or multiout.stat or multibend.stat")
		parser.add_option( "-n", "--name", action="store", dest="name", type="string", default='parfile', help="Supply a name for the output parameter file")
		parser.add_option( "-a", "--automodel", action="store_true", dest="automodel", default=False, help="Run automodel process")
		parser.add_option( "-g", "--gltolerance", action="store", dest="gltolerance", type="float", default=1.0, help="Deviation from ensamble average for global conformation, scale 0-1>")
		parser.add_option( "-b", "--glvariance", action="store", dest="glvariance", type="float", default=0.8, help="Deviation from optimal base-pair and base-pair step parameter for global conformation, scale 0-1>")		
		parser.add_option( "-i", "--glsmoothing", action="store", dest="glsmoothing", type="float", default=0.8, help="Smoothing for twist, roll and tilt, scale 0-1")
		parser.add_option( "-j", "--lcvariance", action="store", dest="lcvariance", type="float", default=0.8, help="Deviation from optimal base-pair and base-pair step parameter for local conformation, scale 0-1>")		
		parser.add_option( "-c", "--number", action="store", dest="number", type="int", default=10, help="Define the number of structures for the automatic modeling process to generate")		
		parser.add_option( "-r", "--refbp", action="store", dest="refbp", type="float", help="Define the reference base-pair over whitch to bend. default is the middel most base-pair")
		parser.add_option( "-o", "--minangle", action="store", dest="minangle", help="Define either the minimum bend angle or give a custom sequence of bend angles")
		parser.add_option( "-p", "--maxangle", action="store", dest="maxangle", type="float", help="Define maximum bend angle")
		parser.add_option( "-q", "--anglestep", action="store", dest="anglestep", type="float", help="Define bend-angle step size")
		parser.add_option( "-d", "--startbp", action="store", dest="startbp", type="int", help="Define start base-pair for bend")
		parser.add_option( "-e", "--endbp", action="store", dest="endbp", type="int", help="Define end base-pair for bend")
		parser.add_option( "-k", "--minorient", action="store", dest="minorient", help="Define either the minimum bend angle in space or give a custom sequence of bend angles")
		parser.add_option( "-l", "--maxorient", action="store", dest="maxorient", type="float", help="Define maximum angle of bend in space")
		parser.add_option( "-m", "--orientstep", action="store", dest="orientstep", type="float", help="Define step size of bend angle in space")
		parser.add_option( "-v", "--verbose", action="store_true", dest="verbose", default=False, help="Output to standard out")
		parser.add_option( "-s", "--bpstep", action="store", dest="bpstep", default="None", help="Supply custom values for base-pair step parameters as a komma seperated list. shift,slide,rise,tilt,roll,twist")
		parser.add_option( "-t", "--bp", action="store", dest="bp", default="None", help="Supply custom values for base-pair parameters as a komma seperated list.")
		parser.add_option( "-u", "--mingroove", action="store", dest="mingroove", type="float", help="Define minimum major groove width")
		parser.add_option( "-w", "--maxgroove", action="store", dest="maxgroove", type="float", help="Define maximum major groove width")
		parser.add_option( "-x", "--groovestep", action="store", dest="groovestep", type="float", help="Define major groove width step size")
		parser.add_option( "-y", "--groovestart", action="store", dest="groovestart", type="int", help="Define start base-pair for major groove width change")
		parser.add_option( "-z", "--grooveend", action="store", dest="grooveend", type="int", help="Define end base-pair for major groove width change")
		parser.add_option( "--helicalphase", action="store", dest="helicalphase", type="float", help="Set the DNA helical phasing (bp/turn)")
				
		(options, args) = parser.parse_args()
		
		self.option_dict['input'] = options.input
		self.option_dict['verbose'] = options.verbose
		self.option_dict['name'] = options.name
		self.option_dict['automodel'] = options.automodel
		self.option_dict['gltolerance'] = options.gltolerance
		self.option_dict['glvariance'] = options.glvariance
		self.option_dict['glsmoothing'] = options.glsmoothing
		self.option_dict['lcvariance'] = options.lcvariance
		self.option_dict['number'] = options.number
		self.option_dict['refbp'] = options.refbp
		self.option_dict['minangle'] = options.minangle	
		self.option_dict['maxangle'] = options.maxangle
		self.option_dict['anglestep'] = options.anglestep
		self.option_dict['startbp'] = options.startbp
		self.option_dict['endbp'] = options.endbp
		self.option_dict['minorient'] = options.minorient
		self.option_dict['maxorient'] = options.maxorient
		self.option_dict['orientstep'] = options.orientstep
		self.option_dict['bpstep'] = options.bpstep
		self.option_dict['bp'] = options.bp
		self.option_dict['mingroove'] = options.mingroove
		self.option_dict['maxgroove'] = options.maxgroove
		self.option_dict['groovestep'] = options.groovestep
		self.option_dict['groovestart'] = options.groovestart
		self.option_dict['grooveend'] = options.grooveend
		self.option_dict['helicalphase'] = options.helicalphase
			
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
		
		if not inputfiles[0] == None:
			for files in inputfiles:
				path = os.path.join(currdir, files)
				filelist.append(path)
			return filelist
		else:
			return inputfiles
	
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

def readNApairing(infile):

    pairdict = {}
    base = None

    infile = open(infile, 'r')
    for line in infile.readlines():
        line = line.strip()
        if line.startswith('Structure'):
            line = line.split()

            try:
                inputfile = os.path.splitext(line[1])[0]
                paircount = int(line[-1].strip(':'))
                pairdict[inputfile] = paircount
            except:
                pass

    for inputfile in pairdict:
        if pairdict[inputfile] == 0:
            base = inputfile
            break

    if not base:
        print "    * ERROR: napairing.stat file indicates that there is no structure without unpaired bases. Unable to set reference file for modelling"
        sys.exit(0)

    return base          

class ModelNucleicAcids:
	
	"""Major class for the modelling of nucleic acid structures"""
	
	def __init__(self,paramdict=None):
	
		self.paramdict = paramdict
		
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
	
	def _FormatParams(self, param, ptype=None):
	
		if param == None:
			return None
		else: 
			try:
				tmp = []
				for n in param.split(','):
					tmp.append(float(n))
				return tmp
			except:
				return [param]
	
	def _MakeDecisions(self, minp=None, maxp=None, pstep=None, def1=None, def2=None):
	
		"""Checks validitie and sets defaults for angle and orient parameters"""
		
		orgminp = self._FormatParams(self.paramdict[minp], minp)
		orgmaxp = self._FormatParams(self.paramdict[maxp], maxp)
		orgpstep = self._FormatParams(self.paramdict[pstep], pstep)
		
		tmpminp = []
		tmpmaxp = []
		tmppstep = []
		
		if orgminp == None and orgmaxp == None:
			print "    * WARNING: No value set for", minp, "and", maxp, "Set to default of:", def1
			tmpminp.append(def1)
			tmpmaxp.append(def1)
			tmppstep.append(None)
		elif orgminp == None and not orgmaxp == None:
			tmpminp = orgmaxp
			tmpmaxp = orgmaxp
			tmppstep.append(None)
			print "    * WARNING: Only value for", maxp, "supplied, set", minp, "equal to", maxp, ":", tmpminp[0] 
		elif not orgminp == None and orgmaxp == None:
			if len(orgminp) == 1:
				tmpminp = orgminp
				tmpmaxp = orgminp
				tmppstep.append(None)
				print "    * Value for",minp,"set to",tmpminp[0]
			elif len(orgminp) == self.paramdict['nrbp']:
				tmpminp = orgminp
				tmpmaxp = orgminp
				tmppstep.append(None)
				print "    * Custom sequence for",minp,"supplied. Set to:",tmpminp
			else:
				tmpminp.append(def1)
				tmpmaxp.append(def1)
				tmppstep.append(None)
				print "    * ERROR: Custom sequence for",minp,"not equal to sequence length. Set to default of", tmpminp[0]
		else:
			if len(orgminp) == 1:
				tmpminp = orgminp
				tmpmaxp = orgmaxp
				if orgpstep == None:
					tmppstep.append(def2)
				else:
					tmppstep = orgpstep
				print "    * Value for",minp," set to:", tmpminp[0], "Value for",maxp,"set to:", tmpmaxp[0]
				print "    * Value for",pstep,"set to:", tmppstep[0]
			elif len(orgminp) == self.nrbp:
				tmpminp = orgminp
				tmpmaxp = orgmaxp
				if orgpstep == None:
					tmppstep.append(def2)
				else:
					tmppstep = orgpstep
				print "    * Custom sequence for",minp," set to:", tmpminp, "Value for",maxp,"set to:", tmpmaxp[0]
				print "    * Value for",pstep,"set to:", tmppstep[0]
				print "    * The values in the custom sequence will be increased", tmppstep[0],"times with a value of:", tmpmaxp[0]	
			else:
				tmpminp.append(def1)
				tmpmaxp.append(def1)
				tmppstep.append(None)
				print "    * ERROR: Custom sequence for",minp,"not equal to sequence length. Set to default of", tmpminp[0]
		
		return tmpminp,tmpmaxp,tmppstep	
	
	def _MakeBaseparDatabase(self,lib):
	
		"""Create the actual database with baseparameter values"""
		
		self.baseparameter = DatabaseDeamon()
		self.baseparameter.Load('sequence',lib[0])
		self.baseparameter.Load('shear',lib[1])
		self.baseparameter.Load('stretch',lib[2])
		self.baseparameter.Load('stagger',lib[3])
		self.baseparameter.Load('buckle',lib[4])
		self.baseparameter.Load('proptw',lib[5])
		self.baseparameter.Load('opening',lib[6])
		self.baseparameter.Load('shift',lib[7])
		self.baseparameter.Load('slide',lib[8])
		self.baseparameter.Load('rise',lib[9])
		self.baseparameter.Load('tilt',lib[10])
		self.baseparameter.Load('roll',lib[11])
		self.baseparameter.Load('twist',lib[12])
	
	def _MakeSwapDatabase(self,sequence):
	
		self.swapdatabase = DatabaseDeamon()
		self.swapdatabase.Load('sequence',sequence)
		
		zerofill = [0.00001]*len(sequence)
		
		basepars = BASEPAIRS+BASEPAIR_STEPS
		
		for parameters in basepars:
			self.swapdatabase.Load(parameters,zerofill)
	
	def _GlobalVariation(self,base,value,sd,variance,scale,smooth=None):
		
		"""Apply variance to the parameter in a weighted manner. More conserved values (measured as amound of variation based on sd) will be introduced
		   with a larger frequency then less conserved parameters. The smoothing factor determines the final value"""
		
		value = array(value)
		sd = array(sd)
		base = array(base)
			
		tsd = array([std(value[1:len(value)])]*len(value))
		weight = (2-pow((sqrt(sd/tsd)),scale))*value
		
		if not smooth == None:
			tmp = []
			maxp = max(weight[1:len(weight)])*smooth 
			minp = min(weight[1:len(weight)])*((1.0-smooth)+1)
			for n in range(len(weight)):
				if weight[n] > maxp or weight[n] < minp:
					tmp.append(weight[n])
				else:
					tmp.append(base[n])
			return (base+(((tmp)-base)*variance)).tolist()
				
		return (base+(((weight)-base)*variance)).tolist()
		
	def _ZeroFill(self,partlist):
	
		"""Fill array shorter than length of sequence with zero's untill the match the length of the sequence"""
	
		begin = [0]*(self.paramdict['startbp']-1)
		end = [0]*(len(self.baseparameter['sequence'])-self.paramdict['endbp'])
	
		return begin+partlist+end
	
	def _MakeArray(self,minl=None,maxl=None,step=None,nondiv=False):
		
		if nondiv == False:
			if len(minl) == 1:
				tmparray = []
				if minl[0] == maxl[0]:
					listl = [minl[0]/(self.paramdict['nrbp']+self.paramdict['nrbp']/10)]*self.paramdict['nrbp']
					tmparray.append(self._ZeroFill(listl))
				else:
					while minl[0] <= maxl[0]:
						listl = [minl[0]/(self.paramdict['nrbp']+self.paramdict['nrbp']/10)]*self.paramdict['nrbp']
						tmparray.append(self._ZeroFill(listl))
						if not step[0] == None:
							minl[0] = minl[0]+step[0]

			else:
				tmparray = []
				if step[0] == None:
					tmparray.append(self._ZeroFill(minl))
				else:
					count = 0
					fixedmaxl = deepcopy(maxl[0])
					while count < step[0]:
						listl = minl+maxl[0]
						tmparray.append(self._ZeroFill(listl))
						count = count+fixedmaxl
						maxl[0] = maxl[0]+fixedmaxl
		else:
			if len(minl) == 1:
				tmparray = []
				if minl[0] == maxl[0]:
					tmparray.append(minl*len(self.baseparameter['sequence']))
				else:
					while minl[0] <= maxl[0]:
						tmparray.append(minl*len(self.baseparameter['sequence']))
						if not step[0] == None:
							minl[0] = minl[0]+step[0]
			else:
				tmparray = []
				tmparray.append(minl)
		
		return tmparray
	
	def _AdjustParams(self, newparam=None, oldparam=None):
	
		"""Replace zero's outside of the calculation zone (startbp-endbp) with the original values in input .par file"""
		
		olddw = oldparam[0:(self.paramdict['startbp']-1)]
		oldup = oldparam[self.paramdict['endbp']:len(oldparam)+1]
		
		return olddw+newparam[self.paramdict['startbp']-1:self.paramdict['endbp']]+oldup	
	
	def _LimitManualModels(self,models):
		
		if float(self.paramdict['number']) > float(MAXMODELS):				# Amount of models that can be generated is not allowed to exeed an maximum
			self.paramdict['number'] = MAXMODELS  							# Maximum is defined in Constants.py (MAXMODELS)	
			print("    * WARNING: number of requested models exceeds allowed maximum of %s. Set to maximum" % MAXMODELS)
		else:
			print("    * Number of models to generate set to: %i" % self.paramdict['number'])
		
		if float(len(models))/float(self.paramdict['number']) < 1.0:		# Amount of models that can be generated cannot be more than number of combinations
			return models													# between angle and orientation
		else:
			tmp = []
			modelrange = range(0,len(models),int(round(float(len(models))/float(self.paramdict['number']))))
			print("    * Take a representative selection of %i models from a total of %i models" % (len(modelrange),len(models)))
			for n in modelrange:
				tmp.append(models[n])
			return tmp

	def _LimitModels(self,anglem,orientm):

		if float(self.paramdict['number']) > float(MAXMODELS):				# Amount of models that can be generated is not allowed to exeed an maximum
			self.paramdict['number'] = MAXMODELS  							# Maximum is defined in Constants.py (MAXMODELS)	
			print("    * WARNING: number of requested models exceeds allowed maximum of %s. Set to maximum" % MAXMODELS)
		else:
			print("    * Number of models to generate set to: %i" % self.paramdict['number'])
			
		models = []	
		if float(len(anglem)*len(orientm))/float(self.paramdict['number']) < 1.0:		# Amount of models that can be generated cannot be more than number of combinations
			for angle in anglem:														# between angle and orientation
				for orient in orientm:
					models.append((angle,orient))
			return models															
		else:
			if int(len(orientm)/round(sqrt(self.paramdict['number']))) == 0:
				opart = 1
			else:
				opart = int(len(orientm)/round(sqrt(self.paramdict['number'])))
			if int(len(anglem)/round(sqrt(self.paramdict['number']))) == 0:
				apart = 1
			else:
				apart = int(len(anglem)/round(sqrt(self.paramdict['number'])))	
		
			anglerange = range(0,len(anglem),apart)
			orientrange = range(0,len(orientm),opart)
			print("    * Take a representative selection of %i models from a total of %i models" % (self.paramdict['number'],len(anglem)*len(orientm)))
			for a in anglerange:
				for o in orientrange:
					models.append((anglem[a],orientm[o]))
			return models

	def _AdjustGroove(self,slide,roll,slidefactor):
			
		"""Adjust the slide to accomondate the user defined major groove width. Set over complete sequence or a user
		   defined zone."""	
		
		rollmax = max(roll)
		newslide = []

		for n in roll:
			try:															# Try to correct for roll. If roll is 0.00 then append slidefactor
				newslide.append((0.0-((n*1.0)/rollmax)+slidefactor))
			except:
				newslide.append(slidefactor)
				
		olddw = slide[0:self.paramdict['groovestart']-1]
		oldup = slide[self.paramdict['grooveend']:len(slide)+1]

		return olddw+newslide[self.paramdict['groovestart']-1:self.paramdict['grooveend']]+oldup
		
	def ReadPar(self,files):
		
		"""Read 3DNA generated base-pair and base-pair step parameter file (*.par) and append
		   to baseparameter database"""
		
		readfile = file(files, 'r')
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
		
		self._MakeBaseparDatabase(lib)
		lib.clear()
	
	def ReadOut(self,files):
		
		"""Read 3DNA generated analysis file (*.out) and append to baseparameter database, Extract 
		   with self._ReadTable"""
		
		lib = {}
		for n in range(13):
			lib[n] = []
		
		bp = re.compile("bp        Shear    Stretch   Stagger    Buckle  Propeller  Opening")
		bpstep = re.compile("step       Shift     Slide      Rise      Tilt      Roll     Twist")
		readfile = file(files, 'r')
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
	
		bpsteplines.insert(0,([0.00001]*8))
	
		for lines in bplines:
			lines = lines[1:8]
			for value in range(7):
				lib[value].append(lines[value])
		for lines in bpsteplines:
			lines = lines[2:8]
			for value in range(6):
				lib[value+7].append(TransformDash(lines[value]))	
		
		del bpsteplines,bplines
		self._MakeBaseparDatabase(lib)
		lib.clear()
		
	def ReadMultibend(self,files):
		
		"""Read multibend file and append parameters to database"""
		
		ref = re.compile("Reference base-pair:")
		start = re.compile("index  bp-step")
		readfile = file(files, 'r')
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
		
		#bendlines.insert(0,([0.00001]*17))
		
		lib = {}
		for n in range(17):
			lib[n] = []
		for lines in bendlines:
			for value in range(17):
				try:
					lib[value].append(TransformDash(lines[value]))
				except:	
					lib[value].append(lines[value])
		self.bendref = DatabaseDeamon()
		self.bendref.Load('refbp', refbp)
		self.bendref.Load('globbend', mean(lib[7]))
		self.bendref.Load('globbendsd', std(lib[7]))
		self.bendref.Load('sequence', lib[1])
		self.bendref.Load('fractilt', lib[3])
		self.bendref.Load('fractiltsd', lib[4])
		self.bendref.Load('fracroll', lib[5])
		self.bendref.Load('fracrollsd', lib[6])
		self.bendref.Load('bpangle', lib[7])
		self.bendref.Load('bpanglesd', lib[8])
		self.bendref.Load('orient', lib[9])
		self.bendref.Load('orientsd', lib[10])
		self.bendref.Load('gtilt', lib[11])
		self.bendref.Load('gtiltsd', lib[12])
		self.bendref.Load('groll', lib[13])
		self.bendref.Load('grollsd', lib[14])
		self.bendref.Load('acctwist', lib[15])
		self.bendref.Load('acctwistsd', lib[16])
		self.bendref.Load('globorient', mean(lib[9]))
		self.bendref.Load('globorientsd', std(lib[9]))
		
		del bendlines
		lib.clear()
		
	def ReadMultiout(self, files):
		
		"""Read statistics for base-pair and base-pair step parameters from multiout file and 
		   append to database"""
		
		lib = {}
		for n in range(26):
			lib[n] = []
		
		bp = re.compile("nr   bp      freq.        shear           stretch")
		bpstep = re.compile("nr   bpstep    freq.      shift           slide")
		readfile = file(files, 'r')
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
					lib[value].append(TransformDash(lines[value]))
				except:
					lib[value].append(lines[value])
		for lines in bpsteplines:
			lines = lines[3:15]
			for value in range(12):
				lib[value+14].append(TransformDash(lines[value]))	
		
		del bpsteplines,bplines
		
		self.baseref = DatabaseDeamon()
		self.baseref.Load('sequence',lib[0])
		self.baseref.Load('shear',lib[2])
		self.baseref.Load('shearsd',lib[3])
		self.baseref.Load('stretch',lib[4])
		self.baseref.Load('stretchsd',lib[5])
		self.baseref.Load('stagger',lib[6])
		self.baseref.Load('staggersd',lib[7])
		self.baseref.Load('buckle',lib[8])
		self.baseref.Load('bucklesd',lib[9])
		self.baseref.Load('proptw',lib[10])
		self.baseref.Load('proptwsd',lib[11])		
		self.baseref.Load('opening',lib[12])
		self.baseref.Load('openingsd',lib[13])
		self.baseref.Load('shift',lib[14])
		self.baseref.Load('shiftsd',lib[15])
		self.baseref.Load('slide',lib[16])
		self.baseref.Load('slidesd',lib[17])
		self.baseref.Load('rise',lib[18])
		self.baseref.Load('risesd',lib[19])
		self.baseref.Load('tilt',lib[20])
		self.baseref.Load('tiltsd',lib[21])
		self.baseref.Load('roll',lib[22])
		self.baseref.Load('rollsd',lib[23])		
		self.baseref.Load('twist',lib[24])
		self.baseref.Load('twistsd',lib[25])
		
		lib.clear()	   	
	
	def CheckParamdict(self):
		
		"""Checking all parameters on there validitie and set defaults"""
		print "    * Checking all parameters on there validitie and set defaults"

		#Check startbp for bending
		if self.paramdict['startbp'] == None:
			self.paramdict['startbp'] = int(1)
			print "    * No start base-pair for bending set, use first base-pair:", self.paramdict['startbp']
		else:
			if self.paramdict['startbp'] > len(self.baseparameter['sequence']):
				self.paramdict['startbp'] = int(1)
				print "    * WARNING: Start base-pair for bending is larger that sequence length. Start base-pair is set to first base-pair:", self.paramdict['startbp']
			else:
				self.paramdict['startbp'] = int(self.paramdict['startbp'])
				print "    * Start base-pair for bending set to:", self.paramdict['startbp']

		#Check startbp for major groove width
		if self.paramdict['groovestart'] == None:
			self.paramdict['groovestart'] = int(1)
			print "    * No start base-pair for major groove width set, use first base-pair:", self.paramdict['groovestart']
		else:
			if self.paramdict['groovestart'] > len(self.baseparameter['sequence']):
				self.paramdict['groovestart'] = int(1)
				print "    * WARNING: Start base-pair for major groove width is larger that sequence length. Start base-pair is set to first base-pair:", self.paramdict['groovestart']
			else:
				self.paramdict['groovestart'] = int(self.paramdict['groovestart'])
				print "    * Start base-pair for major groove width set to:", self.paramdict['groovestart']

		#Check endbp for bending
		if self.paramdict['endbp'] == None:
			self.paramdict['endbp'] = int(len(self.baseparameter['sequence']))
			print "    * No end base-pair for bending set, use last base-pair:", self.paramdict['endbp']
		else:
			if self.paramdict['endbp'] > len(self.baseparameter['sequence']):
				self.endbp = int(len(self.baseparameter['sequence']))
				print "    * WARNING: End base-pair for bending is larger that sequence length. End base-pair is set to last base-pair:", self.endbp
			elif self.paramdict['endbp'] < self.paramdict['startbp']:
				self.paramdict['endbp'] = int(len(self.baseparameter['sequence']))
				print "    * WARNING: End base-pair for bending cannot be smaller than the start base-pair. End base-pair set to last base-pair:", self.endbp
			else:
				self.paramdict['endbp'] = int(self.paramdict['endbp'])
				print "    * End base-pair for bending set to:", self.paramdict['endbp']

		#Check endbp for major groove width
		if self.paramdict['grooveend'] == None:
			self.paramdict['grooveend'] = int(len(self.baseparameter['sequence']))
			print "    * No end base-pair for major groove width set, use last base-pair:", self.paramdict['grooveend']
		else:
			if self.paramdict['grooveend'] > len(self.baseparameter['sequence']):
				self.paramdict['grooveend'] = int(len(self.baseparameter['sequence']))
				print "    * WARNING: End base-pair for groove width is larger that sequence length. End base-pair is set to last base-pair:", self.paramdict['grooveend']
			elif self.paramdict['grooveend'] < self.paramdict['groovestart']:
				self.paramdict['grooveend'] = int(len(self.baseparameter['sequence']))
				print "    * WARNING: End base-pair for major groove width cannot be smaller than the start base-pair. End base-pair set to last base-pair:", self.paramdict['grooveend']
			else:
				self.paramdict['grooveend'] = int(self.paramdict['grooveend'])
				print "    * End base-pair for major groove width set to:", self.paramdict['grooveend']

		#Set the number of base-pairs used in the modelling
		self.paramdict['nrbp'] = (self.paramdict['endbp'] - self.paramdict['startbp'])+1 
		
		#Check reference base-pair
		if self.paramdict['automodel'] == False:
			if self.paramdict['refbp'] == None:
				self.paramdict['refbp'] = float((self.paramdict['nrbp']/2)+self.paramdict['startbp']-1)
				print "    * No value for the reference base-pair used for bending supplied. Set to middle of sequence:", self.paramdict['refbp']
			elif self.paramdict['refbp'] < self.paramdict['startbp'] or self.paramdict['refbp'] > self.paramdict['endbp']:
				self.paramdict['refbp'] = float((self.paramdict['nrbp']/2)+self.paramdict['startbp']-1)
				print "    * WARNING: Reference base-pair has to reside between the start and end base-pair. Set to middle of sequence:", self.paramdict['refbp']
			else:
				print "    * Reference base-pair for bending set to:", self.paramdict['refbp']	
		if self.paramdict['automodel'] == True:
		 	if self.paramdict['refbp'] == None and self.bendref['refbp'] == None:
				self.paramdict['refbp'] = float((self.paramdict['nrbp']/2)+self.paramdict['startbp']-1)
				print "    * No value for the reference base-pair used for bending supplied. Set to middle of sequence:", self.paramdict['refbp']
			elif self.paramdict['refbp'] == None and not self.bendref['refbp'] == None:
				self.paramdict['refbp'] = self.bendref['refbp'][0]
				print "    * Reference base-pair extracted from multibend.stat file set to:", self.paramdict['refbp']
			elif self.paramdict['refbp'] < self.paramdict['startbp'] or self.paramdict['refbp'] > self.paramdict['endbp']:
				self.paramdict['refbp'] = float((self.paramdict['nrbp']/2)+self.paramdict['startbp']-1)
				print "    * WARNING: Reference base-pair has to reside between the start and end base-pair. Set to middle of sequence:", self.paramdict['refbp']
			else:
				print "    * Reference base-pair for bending set to:", self.paramdict['refbp']	

		#Construct slidefactor list from major groove with range
		if not self.paramdict['mingroove'] == None:
			if float(self.paramdict['mingroove']) < 14.0 or float(self.paramdict['mingroove']) > 19.0:
				print "    * Allowed range for major groove width (14.0-19.0 A) violated width value:", self.paramdict['mingroove']
				print "    * Minimum major groove width set to P-P distance of 14 A not refined"	
				self.paramdict['mingroove'] = 14.0
		else:
			print "    * No value for minimum major groove size defined. Not changing anything"
		if not self.paramdict['maxgroove'] == None:
			if float(self.paramdict['maxgroove']) < 14.0 or float(self.paramdict['maxgroove']) > 19.0:
				print "    * Allowed range for major groove width (14.0-19.0 A) violated width value:", self.paramdict['maxgroove']
				print "    * Minimum major groove width set to P-P distance of 19 A not refined"
		else:
			print "    * No value for maximum major groove size defined. Not changing anything"	
			
		if not self.paramdict['mingroove'] == None and not self.paramdict['maxgroove'] == None:
			
			minslide = (1-((float(self.paramdict['mingroove'])-14.0)*2.0)/5.0)
			maxslide = (-1+((19-float(self.paramdict['maxgroove']))*2.0)/5.0)
		
			if not self.paramdict['groovestep'] == None:
				print "    * Major groove step size set to", self.paramdict['groovestep']
				stepsize = (maxslide-minslide)/self.paramdict['groovestep']
			else:
				self.paramdict['groovestep'] = 0.1
				stepsize = 0.1
				print "    * Both minimum and maximum major groove size set but no stepsize. Set to 1 "	
			
			self.paramdict['slidematrix'] = []
			while maxslide <= minslide:
				self.paramdict['slidematrix'].append(maxslide)
				maxslide += stepsize
			
		elif self.paramdict['maxgroove'] == None and not self.paramdict['mingroove'] == None:
			self.paramdict['slidematrix'] = [(1-((float(self.paramdict['mingroove'])-14.0)*2.0)/5.0)]
		elif not self.paramdict['maxgroove'] == None and self.paramdict['mingroove'] == None:
			self.paramdict['slidematrix'] = [(-1+((19-float(self.paramdict['maxgroove']))*2.0)/5.0)]					

		#Construct and validate base-range for angle and orientation
		if self.paramdict['automodel'] == False:
			self.paramdict['minangle'],self.paramdict['maxangle'],self.paramdict['anglestep'] = self._MakeDecisions(minp='minangle',maxp='maxangle',pstep='anglestep',def1=None,def2=1.0)
			self.paramdict['minorient'],self.paramdict['maxorient'],self.paramdict['orientstep'] = self._MakeDecisions(minp='minorient',maxp='maxorient',pstep='orientstep',def1=45.0,def2=1.0)
		if self.paramdict['automodel'] == True:
			if self.paramdict['minangle'] == None and self.paramdict['maxangle'] == None:
				if self.paramdict['gltolerance'] == 0.0:
					self.paramdict['minangle'] = self.bendref['globbend']
					self.paramdict['maxangle'] = self.bendref['globbend']
					self.paramdict['anglestep'] = [1.0]
				else:
					self.paramdict['minangle'] = [((self.bendref['globbend'][0]-self.bendref['globbendsd'][0])*self.paramdict['gltolerance'])]
					self.paramdict['maxangle'] = [((self.bendref['globbend'][0]+self.bendref['globbendsd'][0])*self.paramdict['gltolerance'])]
					self.paramdict['anglestep'] = [(self.paramdict['maxangle'][0]-self.paramdict['minangle'][0])/self.paramdict['number']]
				print "    * No bend-angle preferences supplied. Set automaticly based on range extracted from multibend.stat file. Settings:"
				print("      Minimum bend-angle: %1.2f, maximum bend-angle: %1.2f, angle step: %1.2f" % (self.paramdict['minangle'][0],
				      self.paramdict['maxangle'][0],self.paramdict['anglestep'][0]))
			else:
				print "    * User supplied angle preferences will override automatic defined values from multibend.stat file"
				self.paramdict['minangle'],self.paramdict['maxangle'],self.paramdict['anglestep'] = self._MakeDecisions(minp='minangle',maxp='maxangle',pstep='anglestep',def1=float(0),def2=1.0)
			if self.paramdict['minorient'] == None and self.paramdict['maxorient'] == None:
				if self.paramdict['gltolerance'] == 0.0:
					self.paramdict['minorient'] = self.bendref['globorient']
					self.paramdict['maxorient'] = self.bendref['globorient']
					self.paramdict['orientstep'] = [1.0]
				else:
					self.paramdict['minorient'] = [self.bendref['globorient'][0]-self.bendref['globorientsd'][0]]
					self.paramdict['maxorient'] = [self.bendref['globorient'][0]+self.bendref['globorientsd'][0]]
					self.paramdict['orientstep'] = [(self.paramdict['maxorient'][0]-self.paramdict['minorient'][0])/self.paramdict['number']]
				print "    * No bend-angle orientation preferences supplied. Set automaticly based on range extracted from multibend.stat file. Settings:"
				print("      Minimum bend-angle orientation: %1.2f, maximum bend-angle orientation: %1.2f, orientation step: %1.2f" % (self.paramdict['minorient'][0],
				      self.paramdict['maxorient'][0],self.paramdict['orientstep'][0]))
			else:
				print "    * User supplied angle-orientation preferences will override automatic defined values from multibend.stat file"
				self.paramdict['minorient'],self.paramdict['maxorient'],self.paramdict['orientstep'] = self._MakeDecisions(minp='minorient',maxp='maxorient',pstep='orientstep',def1=float(45),def2=1.0)

		# Setting global and local tolerance and variance parameters
		if self.paramdict['automodel'] == True:
		 	print("    * Global tolerance set to: %1.1f" % self.paramdict['gltolerance'])
			print("    * Global variance set to: %1.1f" % self.paramdict['glvariance'])
			print("    * Global smoothing set to: %1.1f" % self.paramdict['glsmoothing'])
			print("    * Local variance set to: %1.1f" % self.paramdict['lcvariance'])
	
		# Set helical phasing (number of base-pairs per turn)
		if not self.paramdict['helicalphase']:
			print "    * No custom value for helical phasing set (bp/turn)"
		else:
			if self.paramdict['helicalphase'] == 0.0:
				print "    * ERROR: Helical phase cannot be set to 0.0 base-pairs per turn. Skipping"
				self.paramdict['helicalphase'] == None
			else:	
				print("    * Helical phasing set to %1.2f base-pairs per turn" % self.paramdict['helicalphase'])

		# Set custom values for base-pair and base-par step parameters supplied by user
		if not self.paramdict['bpstep'] == None:
			self.paramdict['bpstep'] = self.paramdict['bpstep'].split(',')
			if len(self.paramdict['bpstep']) == 6:
				for param in xrange(len(self.paramdict['bpstep'])):
					try:
						if len(self.paramdict['bpstep'][param]) > 0:
							self.paramdict[BASEPAIR_STEPS[param]] = float(self.paramdict['bpstep'][param])	
							print("    * Custom value for base-pair step parameter %s set to %1.2f" % (BASEPAIR_STEPS[param],self.paramdict[BASEPAIR_STEPS[param]]))
					except:
						pass
			else:
				print "    * WARNING: No (correct) sequence of custom base-pair step values supplied. Skipping"
			del self.paramdict['bpstep']
			
		if not self.paramdict['bp'] == None:
			self.paramdict['bp'] = self.paramdict['bp'].split(',')
			if len(self.paramdict['bp']) == 6:
				for param in xrange(len(self.paramdict['bp'])):
					try:
						if len(self.paramdict['bp'][param]) > 0:
							self.paramdict[BASEPAIRS[param]] = float(self.paramdict['bp'][param])
							print("    * Custom value for base-pair step parameter %s set to %1.2f" % (BASEPAIRS[param],self.paramdict[BASEPAIRS[param]]))	
					except:
						pass
			else:
				print "    * WARNING: No (correct) sequence of custom base-pair values supplied. Skipping"
			del self.paramdict['bp']
	
	def Automodel(self):
	
		
		"""Module for automatic model generation. Modeling process is devided in three parts: 
		
		   Step 1: A default set of base-pair(step) parameters is generated either directly from the sequence or as a user supplied one.
		   Step 2: The global conformation of the ensemble average in terms of bend is modeld in according to the information provided
		           in the multibend.stat file (bend-angle and orientation). Global tolerance allows the user to specify
		           the range of variance around the average values for bend and orientation to be sampled. A global tolerance of 0 only
		           generates a strucure from the average bend and orientation of the ensemble. A global tolerance of 1 samples the complete
		           range of average plus/min the standard deviation for bend and orientation. It is possible to oversample the global
		           orientation by choosing a global tolerance value above 1. In addition to automatic range selection the user can specify
		           a custom range for bend angles and orientation.
		   Step 3: Introduce variation in the base-pair step parameters. The average parameters as extracted from the multiout.stat file
		           are weighted according to there convergens (as indicated by the standard deviation). The weighting is scaled by a fixed
		           scaling factor (defined in the constants.py file). The difference between the base value and weighted value for every 
		           base-pair step parameter can be scaled by the global variance factor and added to the base parameter. 
		   Step 4: Introduce variation in the base-pair parameters. The average parameters as extracted from the multiout.stat file
				   are weighted according to there convergens (as indicated by the standard deviation). The weighting is scaled by a fixed
				   scaling factor (defined in the constants.py file). The difference between the base value and weighted value for every 
				   base-pair step parameter can be scaled by the global variance factor and added to the base parameter. 
		"""
		
		# Check if base parameter file is present, if not try to generate it
		if self.paramdict['base'] == None:
			print "    * WARNING: No base parameter file supplied, try to aquire from sequence in base statistics file"
			seq = ConvertSeq(self.baseref['sequence'])
			FiberModule(''.join(seq.Export('base1')[0]),1,'BDNA','base')
			if os.path.isfile('base.pdb'):
				cmd = "X3DNAanalyze.py -f base.pdb"
				os.system(cmd)
				self.ReadPar('base.par')
			else:
				print "    * ERROR: failed to aquire base parameter file from sequence in base statistics file. Stopping"
				sys.exit(0)	
		
		# Check and Set parameters
		self.CheckParamdict()
		
		# Initiate summery file with data of all modeld structures
		MakeBackup("modelsummery.txt")
		outfile = file("modelsummery.txt","w")
		outfile.write("***************************************************************************************\n")
		outfile.write("Summery file generated nucleic acid models parameter files\n")
		outfile.write("Date/time: %s\n" % ctime())
		outfile.write("***************************************************************************************\n")
		outfile.write("\nUsed parameters\n")
		seq = ""
		for n in self.baseparameter['sequence']:
			seq = seq+n[0]+" "
		outfile.write(" - Nucleic acid sequence: %s\n" % seq)
		outfile.write(" - Number of models: %i\n" % self.paramdict['number'])
		outfile.write(" - Global reference frame: %i\n" % self.paramdict['refbp'])
		outfile.write(" - Zone over witch to bend: start %i, end %i\n" % (self.paramdict['startbp'],self.paramdict['endbp']))
		outfile.write(" - Global nucleic acid bending range: %1.1f to %1.1f degrees\n" % (self.paramdict['minangle'][0],self.paramdict['maxangle'][0]))
		outfile.write(" - Global nucleic acid orientation range: %1.1f to %1.1f degrees\n" % (self.paramdict['minorient'][0],self.paramdict['maxorient'][0]))
		outfile.write(" - Global tolerance set to: %1.1f\n" % self.paramdict['gltolerance'])
		outfile.write(" - Global smoothing set to: %1.1f\n" % self.paramdict['glsmoothing'])
		outfile.write("\n - Twist variance set to: %1.1f     scale: %1.1f\n" % (self.paramdict['glvariance'],TWISTSCALE))
		outfile.write(" - Roll variance set to: %1.1f      scale: %1.1f\n" % (self.paramdict['glvariance'],ROLLSCALE))
		outfile.write(" - Tilt variance set to: %1.1f      scale: %1.1f\n" % (self.paramdict['glvariance'],TILTSCALE))
		outfile.write(" - Rise variance set to: %1.1f      scale: %1.1f\n" % (self.paramdict['glvariance'],RISESCALE))
		outfile.write(" - Slide variance set to: %1.1f     scale: %1.1f\n" % (self.paramdict['glvariance'],SLIDESCALE))
		outfile.write(" - Shift variance set to: %1.1f     scale: %1.1f\n" % (self.paramdict['glvariance'],SHIFTSCALE))
		outfile.write("\n - Shear variance set to: %1.1f     scale: %1.1f\n" % (self.paramdict['lcvariance'],SHEARSCALE))
		outfile.write(" - Stretch variance set to: %1.1f   scale: %1.1f\n" % (self.paramdict['lcvariance'],STRETCHSCALE))
		outfile.write(" - Stagger variance set to: %1.1f   scale: %1.1f\n" % (self.paramdict['lcvariance'],STAGGERSCALE))
		outfile.write(" - Buckle variance set to: %1.1f    scale: %1.1f\n" % (self.paramdict['lcvariance'],BUCKLESCALE))
		outfile.write(" - Prop-Tw variance set to: %1.1f   scale: %1.1f\n" % (self.paramdict['lcvariance'],PROPTWSCALE))
		outfile.write(" - Opening variance set to: %1.1f   scale: %1.1f\n" % (self.paramdict['lcvariance'],OPENINGSCALE))
		outfile.write("\nStructure         GL bend         GL orient\n")
				
		# Step 1. Make Swap Database. All calculations will be done on this database and exported as a new parameter file
		self._MakeSwapDatabase(self.baseparameter['sequence'])
		
		# Start modeling
		struc = 1
		
		# Step 2. Generate smoothly bend global starting conformations
		anglematrix = self._MakeArray(minl=self.paramdict['minangle'],maxl=self.paramdict['maxangle'],step=self.paramdict['anglestep'])
		orientmatrix = self._MakeArray(minl=self.paramdict['minorient'],maxl=self.paramdict['maxorient'],step=self.paramdict['orientstep'],nondiv=True)
		
		print("--> Start generation of parameter files")
		models = self._LimitModels(anglematrix,orientmatrix)
		
		for model in models:
			# Always regenerate the Swap Database
			for parameters in BASEPAIRS+BASEPAIR_STEPS:
				self.swapdatabase.Update(parameters, deepcopy(self.baseparameter[parameters]))
		
			# Add custom values for base-pair and base-pair step parameters to swapdatabase
			for param in BASEPAIRS+BASEPAIR_STEPS:
				if self.paramdict.has_key(param):
					arr = [self.paramdict[param]]*len(self.swapdatabase['sequence'])
					self.swapdatabase.Update(param,arr)
			
			acctw = AccTwist(self.paramdict['refbp'],self.swapdatabase['twist'])
			directions = TwistCorrect(acctw,model[1])

			vector = []
			for direction in range(len(directions)):	
				frac = (DegreeToUnitvec(directions[direction]))
				vector.append(AngleVector(model[0][direction],frac[0],frac[1]))

			locroll = []
			loctilt = []
			for n in vector:
				locroll.append(n[1])
				loctilt.append(n[0])
			
			self.swapdatabase.Update('tilt',loctilt)
			self.swapdatabase.Update('roll',locroll)
			
			# Step 3. Introducing variation in base-pair step parameters
			vartwist = self._GlobalVariation(self.swapdatabase['twist'],self.baseref['twist'],self.baseref['twistsd'],self.paramdict['glvariance'],TWISTSCALE,smooth=self.paramdict['glsmoothing'])
			varroll = self._GlobalVariation(self.swapdatabase['roll'],self.baseref['roll'],self.baseref['rollsd'],self.paramdict['glvariance'],ROLLSCALE,smooth=self.paramdict['glsmoothing'])
			vartilt = self._GlobalVariation(self.swapdatabase['tilt'],self.baseref['tilt'],self.baseref['tiltsd'],self.paramdict['glvariance'],TILTSCALE,smooth=self.paramdict['glsmoothing'])
			varrise = self._GlobalVariation(self.swapdatabase['rise'],self.baseref['rise'],self.baseref['risesd'],self.paramdict['glvariance'],RISESCALE)
			varshift = self._GlobalVariation(self.swapdatabase['shift'],self.baseref['shift'],self.baseref['shiftsd'],self.paramdict['glvariance'],SHIFTSCALE)
			varslide = self._GlobalVariation(self.swapdatabase['slide'],self.baseref['slide'],self.baseref['slidesd'],self.paramdict['glvariance'],SLIDESCALE)
			
			self.swapdatabase.Update('twist',vartwist)
			self.swapdatabase.Update('roll',varroll)
			self.swapdatabase.Update('tilt',vartilt)
			self.swapdatabase.Update('rise',varrise)
			self.swapdatabase.Update('shift',varshift)
			self.swapdatabase.Update('slide',varslide)
			
			# Step 4. Introduce variation in base-pair parameters
			varshear = self._GlobalVariation(self.swapdatabase['shear'],self.baseref['shear'],self.baseref['shearsd'],self.paramdict['lcvariance'],SHEARSCALE)
			varstretch = self._GlobalVariation(self.swapdatabase['stretch'],self.baseref['stretch'],self.baseref['stretchsd'],self.paramdict['lcvariance'],STRETCHSCALE)
			varstagger = self._GlobalVariation(self.swapdatabase['stagger'],self.baseref['stagger'],self.baseref['staggersd'],self.paramdict['lcvariance'],STAGGERSCALE)
			varbuckle = self._GlobalVariation(self.swapdatabase['buckle'],self.baseref['buckle'],self.baseref['bucklesd'],self.paramdict['lcvariance'],BUCKLESCALE)
			varproptw = self._GlobalVariation(self.swapdatabase['proptw'],self.baseref['proptw'],self.baseref['proptwsd'],self.paramdict['lcvariance'],PROPTWSCALE)
			varopening = self._GlobalVariation(self.swapdatabase['opening'],self.baseref['opening'],self.baseref['openingsd'],self.paramdict['lcvariance'],OPENINGSCALE)
			
			self.swapdatabase.Update('shear',varshear)
			self.swapdatabase.Update('stretch',varstretch)
			self.swapdatabase.Update('stagger',varstagger)
			self.swapdatabase.Update('buckle',varbuckle)
			self.swapdatabase.Update('proptw',varproptw)
			self.swapdatabase.Update('opening',varopening)
			
			# Step 5. Write new parameter file
			outfile.write("%s      %1.1f      %1.1f\n" % (self.paramdict['name']+str(struc),sum(model[0]),mean(model[1])))
			WritePar(self.swapdatabase,self.paramdict['name']+str(struc),self.paramdict['verbose'])
			struc += 1
			
		outfile.close()		

	def Manualmodel(self):
	
		"""Manual modeling module. Modeling process is devided into three parts: using the validated input parameters a list or arrays is made
		   for angles and the orientation of the bend in space for all models that need to be made based on the parameters "anglestep" and/or 
		   "orientstep". Second a swapdatabase is made all calculations will be done on this database and exported as a new parameter file. Afterwards 
		   the actual modeling process starts the bend angles will be introduced, the swapdatabase will be updated with the resulting values for 
		   roll, tilt and twist and a new parameter file will be written.
		"""
		
		if self.paramdict['name'] == None:
			self.paramdict['name'] = (os.path.splitext(os.path.basename(self.paramdict['base']))[0])
			
		# 1. Make list of arrays   
		if self.paramdict['minangle'][0] == None and self.paramdict['maxangle'][0] == None:
			anglematrix = [self.paramdict['minangle']]
		else:	
			anglematrix = self._MakeArray(minl=self.paramdict['minangle'],maxl=self.paramdict['maxangle'],step=self.paramdict['anglestep'])
		
		orientmatrix = self._MakeArray(minl=self.paramdict['minorient'],maxl=self.paramdict['maxorient'],step=self.paramdict['orientstep'],nondiv=True)
	
		# 2. Make Swap Database
		self._MakeSwapDatabase(self.baseparameter['sequence'])
		basepars = BASEPAIRS+BASEPAIR_STEPS
		
		for parameters in basepars:
			self.swapdatabase.Update(parameters, deepcopy(self.baseparameter[parameters]))

		# 3. Set helical phasing (number of basepairs per turn)
		if not self.paramdict['helicalphase'] == None:
			twistperbp = 360/self.paramdict['helicalphase']
			arr = [twistperbp]*len(self.swapdatabase['sequence'])
			arr[0] = 0.0
			self.swapdatabase.Update('twist',arr)
	
		# 4. Add custom values for base-pair and base-pair step parameters to swapdatabase
		for param in BASEPAIRS+BASEPAIR_STEPS:
			if self.paramdict.has_key(param):
				arr = [self.paramdict[param]]*len(self.swapdatabase['sequence'])
				self.swapdatabase.Update(param,arr)
		
		# 5. Start modeling
		struc = 1
		
		print("--> Start generation of parameter files")
		models = []
		if self.paramdict.has_key('slidematrix'):
			for angle in anglematrix:
				for orient in orientmatrix:
					for slidefactor in self.paramdict['slidematrix']:
						models.append((angle,orient,slidefactor))
		else:
			for angle in anglematrix:
				for orient in orientmatrix:
					models.append((angle,orient))				

		models = self._LimitManualModels(models)
	
		# 4. Model in features
		for model in models:
			
			if not model[0][0] == None:
				#calculated accumulated twist values
				acctw = AccTwist(self.paramdict['refbp'],self.baseparameter['twist'])

				#correct global bend angles for twist values
				directions = TwistCorrect(acctw,model[1])
			
				#Calculate local roll and twist values
				vector = []
				for direction in range(len(directions)):	
					frac = (DegreeToUnitvec(directions[direction]))
					vector.append(AngleVector(model[0][direction],frac[0],frac[1]))
			
				#Construct new list of local roll and tilt values
				locroll = []
				loctilt = []
				for n in vector:
					locroll.append(n[1])
					loctilt.append(n[0])

				#Adjust new roll and tilt parameters for zone and add to database
				self.swapdatabase.Update('roll',self._AdjustParams(newparam=locroll,oldparam=self.swapdatabase['roll']))
				self.swapdatabase.Update('tilt',self._AdjustParams(newparam=loctilt,oldparam=self.swapdatabase['tilt']))
			
			#Construct new slide array to accomodate user defined major groove size changes.
			if self.paramdict.has_key('slidematrix'):
				slide = self._AdjustGroove(self.swapdatabase['slide'],self.swapdatabase['roll'],model[2])
				self.swapdatabase.Update('slide',slide)
				
			#Write new parameter file
			WritePar(self.swapdatabase,self.paramdict['name']+str(struc),self.paramdict['verbose'])
			struc += 1
				
if __name__ == '__main__':

	"""Running this plugin from the command line"""
	
	"""Import command line specific modules"""
	from optparse import *
	
	"""Setting up parameter dictionary"""
	paramdict = CommandlineOptionParser().option_dict
	
	"""Check for input"""
	if paramdict['input'] == None:
		print "    * Please supply .par and/or multibend.stat, multiout.stat file using option -b or use option -h/--help for usage"
		sys.exit(0)
	else:
		inputlist = paramdict['input']
	
	"""Running plugin main functions"""
	PluginCore(paramdict, inputlist)
	sys.exit(0)
