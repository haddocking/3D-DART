#!/usr/bin/env python2.7

USAGE= """
==========================================================================================

Author:				Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center for 
					Biomolecular Research, Utrecht university, The Netherlands.
Copyright (C):		2006 (DART project)
DART version:		1.2 (25-11-2008)
DART plugin: 		NABendAnalyze.py
Input:				3DNA PAR files or files containing a list of files. Any of the 
					allowed options, any order or combined.
Output:				.bend and/or a multibend.stat file.			
Plugin excecution:	Either command line driven (use -h/--help for the option) or as 
					part of a DART batch sequence.
Plugin function:	This plugin calculates the global bend angle of nucleic-acid
					structures. Calculation relies on the twist, roll and tilt
					values extracted from the .out files of a 3DNA analysis run.
					You can define the zone over which to calculate the bend and/or
					the reference base-pair to be used. Additionaly the script can 
					calculate statistical values for all bend parameters for an 
					ensamble of structures.
Dependencies:		Standard python2.3 or higher modules
Examples:			NABendAnalyze -f *.out
					NABendAnalyze -mf selection.list	

==========================================================================================
"""

"""Import modules"""
import os, sys, math, re, operator
from copy import deepcopy
from numpy import *
from time import ctime

"""Setting pythonpath variables if run from the command line"""
base, dirs = os.path.split(os.path.dirname(os.path.join(os.getcwd(), __file__)))

if base in sys.path: pass
else: sys.path.append(base)

"""Import DART specific modules"""
from QueryPDB import GetSequence, NAsummery
from PDBeditor import PDBeditor
from system.NAfunctionLib import ConvertSeq, UnitvecToDegree, AccTwist, Angle
from system.IOlib import InputOutputControl
from system.Constants import *

def PluginXML():
	PluginXML = """ 
<metadata>
 <name>Nucleic Acid global bend analysis</name>
 <input type="filetype">.out</input>
 <output type="filetype">.bend,multibend.stat</output>
</metadata>
<parameters>
 <option type="useplugin" form="hidden" text="None">True</option>
 <option type="inputfrom" form="hidden" text="None">2</option>
 <option type="refbp" form="text" text="Define base-pair used as global origin"></option>
 <option type="zone" form="text" text="Define zone over which to calculate"></option>
 <option type="multiana" form="checkbox" text="Perform multibend analysis">True</option>
 <option type="verbose" form="checkbox" text="Verbose output">False</option>
</parameters>"""
	
	return PluginXML

def PluginCore(paramdict, inputlist):
	
	"""Checking inputlist"""
	checked = InputOutputControl()
	checked.CheckInput(inputlist)
	
	print "--> Performing Global nucleic acid bend angle analysis"
	if paramdict['refbp']:
		print "    * Using predefined reference base-pair:", paramdict['refbp']
	if paramdict['zone']:
		print("    * Using predefined zone over whitch to calculate the bend angle: %s" % paramdict['zone'])
	
	bend = MeasureBend(refbp=paramdict['refbp'], zone=paramdict['zone'], verbose=paramdict['verbose'])	

	try:
		bend.ReadOutfiles(files=checked.checkedinput['.out'])
		if paramdict['multiana'] == True:
			bend.GetBaseseq()
			bend.CalcGlobalBend(multiana=paramdict['multiana'])
		else:
			bend.CalcGlobalBend(multiana=False)
	except:
		try:
			bend.ReadParfiles(files=checked.checkedinput['.par'])
			bend.CalcGlobalBend(multiana=False)
		except e:
			print "    * ERROR: No valid input found: {}".format(e)
			sys.exit(0)
	
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
		parser = OptionParser(usage)

		parser.add_option( "-f", "--file", action="callback", callback=self.varargs, dest="inputfile", type="string", help="Supply .par or .out inputfile(s). In case of .par no multianalysis")
		parser.add_option( "-r", "--refbp", action="store", dest="refbp", type="float", help="Define the reference base-pair over which to calculate the bend angle")
		parser.add_option( "-z", "--zone", action="store", dest="zone", help="Define the zone over which to calculate the bend angle")
		parser.add_option( "-m", "--multianalyze", action="store_true", dest="multiana", default=False, help="Perform multi-bend analysis over supplied files ")
		parser.add_option( "-v", "--verbose", action="store_true", dest="verbose", default=False, help="Output to standard out")
		
		(options, args) = parser.parse_args()
		
		self.option_dict['input'] = options.inputfile
		self.option_dict['refbp'] = options.refbp
		self.option_dict['zone'] = options.zone
		self.option_dict['multiana'] = options.multiana
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

class MeasureBend:

	"""Measure the global bend angle of irregular nucleic acids"""
	
	def __init__(self, refbp=None, zone=None, verbose=False):
		
		self.refbp = refbp		#user supplied or auto-set
		self.zone = zone		#user supplied or auto-set
		self.verbose = verbose
		
		self.bpstep = {}		#extracted from .out file
		self.pairs = {}			#extracted from .out file
		self.chainid = ''		#extracted from .out file
		
		self.tilt = {}
		self.roll = {}
		self.twist = {}
		self.freq = {}
	
		self.basenr = 0
		self.basemoltype = {}
		self.basechainlib = {}
		self.basesequence = {}
	
	def _RefBp(self, sequence):
		
		"""Determine reference base-pair if none is supplied. By default choose middel of zone"""
		
		if self.refbp == None:
			
			corr_zone = self._Zone(sequence)
			nrbp = self._NrBp(sequence)
			
			half = float(nrbp)/2
			print "    * No reference base-pair for bend calculation supplied. Set to middle of zone:",  (corr_zone[0]-1) + int(half)
				
			return (corr_zone[0]-1) + int(half)
		else:
			corr_zone = self._Zone(sequence)
			if self.refbp > corr_zone[1]:
				print "    * WARNING: supplied reference bp larger that sequence length, refbp set to middle of sequence"
				self.refbp = None
				self._RefBp(sequence)
			else:
				return self.refbp-1
		
	def _NrBp(self, sequence):
		
		"""Determine number of base-pairs in sequence. If zone=None then full sequence is used"""
		
		corr_zone = self._Zone(sequence)
		return ((corr_zone[1])+1) - (corr_zone[0])
	
	def _Zone(self, sequence):
	
		"""Check zone, if zone=None then full sequence is returned"""
		
		if self.zone == None:
			return [1, sequence]
		else:
			try:	
				corr_zone = self.zone.split(',')
				if int(corr_zone[0]) > sequence or int(corr_zone[1]) > sequence:
					print "    * WARNING: one or more base-pairs in supplied zone are outside the length of the sequence, zone set to full sequence length"
					return [1, sequence]
				else:	
					return [int(corr_zone[0]), int(corr_zone[1])]
			except:
				corr_zone = self.zone.split('-')
				if int(corr_zone[0]) > sequence or int(corr_zone[1]) > sequence:
					print "    * WARNING: one or more base-pairs in supplied zone are outside the length of the sequence, zone set to full sequence length"
					return [1, sequence]
				else:	
					return [int(corr_zone[0]), int(corr_zone[1])]
	
	def _GlobRollTilt(self, tilt, roll):		
	
		"""Calculate global roll and tilt"""
		
		base = int(self.corr_refbp)					
		frup = self.corr_refbp - base
		frdw = 1-frup
		p = math.pi/180
		
		if frup == 0:
			tiltarr = array(tilt)					
    			rollarr = array(roll)	
   			ctwistarr = array(self.ctwist)

			gltiltarr =  (tiltarr*cos(ctwistarr*p))+(rollarr*sin(ctwistarr*p))
    			glrollarr = (tiltarr*-sin(ctwistarr*p))+(rollarr*cos(ctwistarr*p))

			return gltiltarr.tolist(), glrollarr.tolist()		

		else:
			tilt.insert(base,(tilt[base]))
			tilt[base] = tilt[base]/2
			tilt[base+1] = tilt[base+1]/2
			
			roll.insert(base,(roll[base]))
			roll[base] = roll[base]/2
			roll[base+1] = roll[base+1]/2
			
    			tiltarr = array(tilt)					
    			rollarr = array(roll)	
   			ctwistarr = array(self.ctwist)
			
			gltiltarr =  (tiltarr*cos(ctwistarr*p))+(rollarr*sin(ctwistarr*p))
    			glrollarr = (tiltarr*-sin(ctwistarr*p))+(rollarr*cos(ctwistarr*p))
			
			self.corr_zone[1] = self.corr_zone[1]+1 #correction for duplicated base-pair
			
			return gltiltarr.tolist(), glrollarr.tolist()	
	
	def _ListFraction(self, v1,v2,v3):			#Calculation of angle fractions
		a = array(v1)
		b = array(v2)
		c = array(v3)
		d = (power(v2,2))/(power(v1,2))
		e = (power(v3,2))/(power(v1,2))

		f = d.tolist()
		g = e.tolist()

		return f,g
	
	def _ZoneCorrect(self, inlist):
		nrbp = len(inlist)
		nul_list = [0.0] * nrbp
		
		down = self.corr_zone[0]-1
		up = self.corr_zone[1]
		
		while down < up:
			nul_list[down] = inlist[down]
			down = down + 1
		
		return nul_list
	
	def _ConvertSeq2(self,inlist):
		
		dupl = deepcopy(inlist)
		full = len(self.basechainlib[self.chainid][0])
		count = 1
		while count < full:
			if self.basechainlib[self.chainid][0][count] == dupl[3][count]:
				pass
			else:
				dupl[0].insert(count,"X")
				dupl[1].insert(count,"X")
				dupl[2].insert(count,"X")
				dupl[3].insert(count,"X")
				dupl[4].insert(count,"X")
			count = count +1	
		
		newlist = []
		for na in range(len(dupl[1])):
			try:	
				pair = []
				if dupl[1][na] in baselist3 and dupl[1][na+1] in baselist3:
					pair.append(pairlist1_t[baselist3.index(dupl[1][na])])	
					pair.append(pairlist1_t[baselist3.index(dupl[1][na+1])])
				pair.append('/')
				pair.append(pairlist1_c[baselist3.index(dupl[1][na+1])])
				pair.append(pairlist1_c[baselist3.index(dupl[1][na])])
				pair = ''.join(pair)	
				newlist.append(pair)
			except:
				pass		
		
		newlist2 = []
		count = 1
		for na in range(len(dupl[3])):
			try:
				false = 0
				if dupl[3][na] == "X":
					false = 1
				if dupl[3][na+1] == "X":
					if false == 1:
						pass
					else:
						false = 1
				if false == 0:
					newlist2.append(count)
				count += 1
			except:
				pass
	
		return [newlist, newlist2]
	
	def _CalculateCore(self, files):
		
		"""Get reference basepair"""
		self.corr_refbp = self._RefBp(len(self.bpstep[files][0]))
		self.corr_zone = self._Zone(len(self.bpstep[files][0]))
		self.corr_nrbp = self._NrBp(len(self.bpstep[files][0]))
	
		"""Calculate accumulated twist, global roll and tilt, base-pair angle at each step and fraction roll/tilt"""
		self.ctwist = AccTwist(self.corr_refbp,self.bpstep[files][4])
		(self.gltilt, self.glroll) = self._GlobRollTilt(self.bpstep[files][2],self.bpstep[files][3])
		self.bpangle = Angle(self.gltilt,self.glroll) 
		(self.ftilt, self.froll) = self._ListFraction(self.bpangle,self.gltilt,self.glroll) 
	
		"""Correct above calculated values for the selected zone. return value 0 for resids outside zone"""
		self.ctwist = self._ZoneCorrect(self.ctwist)
		self.gltilt = self._ZoneCorrect(self.gltilt)
		self.glroll = self._ZoneCorrect(self.glroll)
		self.bpangle = self._ZoneCorrect(self.bpangle)
		self.ftilt = self._ZoneCorrect(self.ftilt)
		self.froll = self._ZoneCorrect(self.froll)
	
		"""Sum up all parameters"""	
		self.sgroll = sum(self.glroll)
		self.sgtilt = sum(self.gltilt)
		self.sftilt = sum(self.ftilt)
		self.sfroll = sum(self.froll)
		self.sbpangle = sum(self.bpangle)
		self.sctwist = sum(self.ctwist)

		"""Calculate average of all summed parameters"""	
		self.avgroll = self.sgroll/(self.corr_nrbp)
		self.avgtilt = self.sgtilt/(self.corr_nrbp)
		self.avftilt = self.sftilt/(self.corr_nrbp)
		self.avfroll = self.sfroll/(self.corr_nrbp)
		self.avbpangle = self.sbpangle/(self.corr_nrbp)
		self.avctwist = self.sctwist/(self.corr_nrbp)

		"""Calculate standard deviation of all summed parameters"""
		self.sdgroll = std(self.glroll[self.corr_zone[0]:self.corr_zone[1]])
		self.sdgtilt = std(self.gltilt[self.corr_zone[0]:self.corr_zone[1]])
		self.sdftilt = std(self.ftilt[self.corr_zone[0]:self.corr_zone[1]])
		self.sdfroll = std(self.froll[self.corr_zone[0]:self.corr_zone[1]])
		self.sdbpangle = std(self.bpangle[self.corr_zone[0]:self.corr_zone[1]])
		self.sdctwist = std(self.ctwist[self.corr_zone[0]:self.corr_zone[1]])

		"""Calculate global bend"""
		self.gbend = Angle(self.sgroll,self.sgtilt)
		
		"""Calculate orientation of global bend angle in Euler-Angle space"""
		self.orient = []
		for n in range(len(self.glroll)):
			self.orient.append(UnitvecToDegree(self.gltilt[n],self.glroll[n],self.bpangle[n]))
		self.sorient = sum(self.orient)
		self.avorient = self.sorient/(self.corr_nrbp)
		self.sdorient = std(self.orient[self.corr_zone[0]:self.corr_zone[1]])
				
		"""Report position of reference base-pair"""
		base = int(self.corr_refbp)					
		frup = self.corr_refbp-base
		frdw = 1-frup

		if frup == 0:
			self.emtylist = ['']*len(self.bpstep[files][0])
			self.emtylist[base] = 'x'
		else:	
			self.bpstep[files][1].insert(base,(self.bpstep[files][1][base]))

			self.emtylist = ['']*len(self.bpstep[files][1])
			self.emtylist[base] = 'x'
			self.emtylist[base+1] = 'x' 
	
	def _MultiCalculate(self):
		
		"""Get reference basepair"""
		self.corr_refbp = self._RefBp(len(self.basesequence[self.chainid]))
		self.corr_zone = self._Zone(len(self.basesequence[self.chainid]))
		self.corr_nrbp = self._NrBp(len(self.basesequence[self.chainid]))
		
		"""Calculate averages"""
		roll = []
		tilt = []
		twist = []
		for n in range(len(self.roll)):
			if len(self.roll[n]) > 1:
				roll.append(mean(self.roll[n]))
				tilt.append(mean(self.tilt[n]))
				twist.append(mean(self.twist[n]))
			else:
				roll.append(float(0))	
				tilt.append(float(0))	
				twist.append(float(0))	
		
		"""Calculate standard deviations"""
		stroll = []
		sttilt = []
		sttwist = []
		for n in range(len(self.roll)):
			if len(self.roll[n]) > 1:
				stroll.append(std(self.roll[n]))
				sttilt.append(std(self.tilt[n]))
				sttwist.append(std(self.twist[n]))
			else:
				stroll.append(float(0))	
				sttilt.append(float(0))	
				sttwist.append(float(0))	
		
		rolllist = []
		rolllist.append(array(roll))
		rolllist.append(array(roll)+array(stroll))
		rolllist.append(array(roll)-array(stroll))
		
		tiltlist = []
		tiltlist.append(array(tilt))
		tiltlist.append(array(tilt)+array(sttilt))
		tiltlist.append(array(tilt)-array(sttilt))
		
		twistlist = []
		twistlist.append(array(twist))
		twistlist.append(array(twist)+array(sttwist))
		twistlist.append(array(twist)-array(sttwist))
		
		"""Calculate accumulated twist, global roll and tilt, base-pair angle and angle orientation at each step and 
		   fraction roll/tilt"""
		
		self.ctwist = []
		for n in twistlist:
			self.ctwist.append(AccTwist(self.corr_refbp,n.tolist()))
		
		self.gltilt = []
		self.glroll = []
		for n in range(len(rolllist)):
			(gltilt,glroll) = self._GlobRollTilt(tiltlist[n].tolist(),rolllist[n].tolist())
			self.gltilt.append(gltilt[0])
			self.glroll.append(glroll[0])
		
		self.bpangle = []
		for n in range(len(self.gltilt)):
			self.bpangle.append(Angle(self.gltilt[n],self.glroll[n]))
		
		self.ftilt = []
		self.froll = []
		for n in range(len(self.bpangle)):	
			(ftilt, froll) = self._ListFraction(self.bpangle[n],self.gltilt[n],self.glroll[n])
			self.ftilt.append(ftilt)
			self.froll.append(froll) 
		
		self.orient = []
		for n in range(len(self.glroll)):
			tmp = []
			for k in range(len(self.glroll[n])):	
				tmp.append(UnitvecToDegree(self.gltilt[n][k],self.glroll[n][k],self.bpangle[n][k]))
			self.orient.append(tmp)	
		
		"""Supstract upper and lower limits to obtain interval"""
		self.ctwist = self.ctwist[0],(abs(array(self.ctwist[1])-array(self.ctwist[2]))).tolist()
		self.gltilt = self.gltilt[0],(abs(array(self.gltilt[1])-array(self.gltilt[2]))).tolist()
		self.glroll = self.glroll[0],(abs(array(self.glroll[1])-array(self.glroll[2]))).tolist()
		self.bpangle = self.bpangle[0],(abs(array(self.bpangle[1])-array(self.bpangle[2]))).tolist()
		self.orient = self.orient[0],(abs(array(self.orient[1])-array(self.orient[2]))).tolist()
		self.ftilt = self.ftilt[0],(abs(array(self.ftilt[1])-array(self.ftilt[2]))).tolist()
		self.froll = self.froll[0],(abs(array(self.froll[1])-array(self.froll[2]))).tolist()
		
		"""Correct above calculated values for the selected zone. return value 0 for resids outside zone"""
		self.ctwist = self._ZoneCorrect(self.ctwist[0]),self._ZoneCorrect(self.ctwist[1])  
		self.gltilt = self._ZoneCorrect(self.gltilt[0]),self._ZoneCorrect(self.gltilt[1])  
		self.glroll = self._ZoneCorrect(self.glroll[0]),self._ZoneCorrect(self.glroll[1])  
		self.bpangle = self._ZoneCorrect(self.bpangle[0]),self._ZoneCorrect(self.bpangle[1])
		self.orient = self._ZoneCorrect(self.orient[0]),self._ZoneCorrect(self.orient[1])
		self.ftilt = self._ZoneCorrect(self.ftilt[0]),self._ZoneCorrect(self.ftilt[1])	
		self.froll = self._ZoneCorrect(self.froll[0]),self._ZoneCorrect(self.froll[1])	
		
		"""Sum up all parameters"""	
		self.sgroll = sum(self.glroll[0])
		self.sgtilt = sum(self.gltilt[0])
		self.sftilt = sum(self.ftilt[0])
		self.sfroll = sum(self.froll[0])
		self.sbpangle = sum(self.bpangle[0])
		self.sorient = sum(self.orient[0])
		self.sctwist = sum(self.ctwist[0])

		"""Calculate average of all summed parameters"""	
		self.avgroll = self.sgroll/(self.corr_nrbp)
		self.avgtilt = self.sgtilt/(self.corr_nrbp)
		self.avftilt = self.sftilt/(self.corr_nrbp)
		self.avfroll = self.sfroll/(self.corr_nrbp)
		self.avbpangle = self.sbpangle/(self.corr_nrbp)
		self.avorient = self.sorient/(self.corr_nrbp)
		self.avctwist = self.sctwist/(self.corr_nrbp)

		"""Calculate standard deviation of all summed parameters"""
		self.sdgroll = std(self.glroll[0][self.corr_zone[0]:self.corr_zone[1]])
		self.sdgtilt = std(self.gltilt[0][self.corr_zone[0]:self.corr_zone[1]])
		self.sdftilt = std(self.ftilt[0][self.corr_zone[0]:self.corr_zone[1]])
		self.sdfroll = std(self.froll[0][self.corr_zone[0]:self.corr_zone[1]])
		self.sdbpangle = std(self.bpangle[0][self.corr_zone[0]:self.corr_zone[1]])
		self.sdorient = std(self.orient[0][self.corr_zone[0]:self.corr_zone[1]])
		self.sdctwist = std(self.ctwist[0][self.corr_zone[0]:self.corr_zone[1]])

		"""Calculate global bend"""
		self.gbend = Angle(self.sgroll,self.sgtilt)
		self.stgbend = Angle(self.sdgroll,self.sdgtilt)

		"""Report position of reference base-pair"""
		base = int(self.corr_refbp)					
		frup = self.corr_refbp-base
		frdw = 1-frup

		if frup == 0:
			self.emtylist = ['']*len(self.basesequence[self.chainid])
			self.emtylist[base] = 'x'
		else:	
			self.basesequence[self.chainid].insert(base,(self.basesequence[self.chainid][base]))

			self.emtylist = ['']*len(self.basesequence[self.chainid])
			self.emtylist[base] = 'x'
			self.emtylist[base+1] = 'x' 
		
		"""Correct frequency for duplicated base-pairs"""	
		freq = []
		for n in range(len(self.freq)):
			freq.append(self.freq[n])
		
		if frup == 0:
			self.freq = freq
		else:
			freq.insert(base,(self.freq[base]))
			self.freq = freq			
		
	def _WriteMultiBend(self):
		
		"""Print all data to file"""
		
		print("    * Writing Multi-bend analysis data to file: multibend.stat")
		
		if self.verbose == True:
			outfile = sys.stdout 
		else:
			outfile = file('multibend.stat','w')	
		
		outfile.write("**********************************************************************************************************************************************************************\n")
		outfile.write("Statistical info for the global bend analysis of %i files\n" % len(self.pairs))	
		outfile.write("Date/time: %s\n" % ctime())
		outfile.write("Reference base-pair: %.2f\n" % (self.corr_refbp+1))
		outfile.write("Zone: start %i end %i\n" % (self.corr_zone[0],self.corr_zone[1]))
		outfile.write("**********************************************************************************************************************************************************************\n")	
		outfile.write("Global bend angles of individual base-pair steps and of all identified base-pair steps together. Base-pair steps of unpaired\n")
		outfile.write("or mispaired base-pairs steps are not reported.\n")
		outfile.write("freq.         = Frequency in which the given base-pair step occurs in the ensemble of structures\n")
		outfile.write("frac. tilt    = contribution of tilt to the base-pair step bending angle (mean and standard deviation)\n")
		outfile.write("frac. roll    = contribution of roll to the base-pair step bending angle (mean and standard deviation)\n")
		outfile.write("bp-step angle = bending angle (deg.) at each base-pair step (mean and standard deviation)\n") 
		outfile.write("orient angle  = orientation of the bend angle relative to the reference base-pair in Euler-angle\n")
		outfile.write("                space (from 0 to 360 deg.)\n")  
		outfile.write("Acc. twist    = accumulated twist values above and below the reference base-pair (mean and standard deviation)\n")
		outfile.write("x	      = If the reference plane is located between two basepairs, one base-pair is\n") 
		outfile.write("		duplicated for the calculation (two x's). If the reference plane coincides with\n") 
		outfile.write("		one base-pair no duplication is necessary (one x)\n")
		outfile.write("\nindex  bp-step  freq.  frac.tilt std.      frac.roll std.   bp-step angle  std.  orient angle  std. global tilt   std.  global roll   std.   Acc.twist  std.   refbp\n")
		
		for p in range(len(self.basesequence[self.chainid])):
			outfile.write("%2d %9s %5.1d %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %8.3f %3s\n" % (p+1,
			          self.basesequence[self.chainid][p],self.freq[p],self.ftilt[0][p],self.ftilt[1][p],self.froll[0][p],self.froll[1][p],
			          self.bpangle[0][p],self.bpangle[1][p],self.orient[0][p],self.orient[1][p],self.gltilt[0][p],self.gltilt[1][p],self.glroll[0][p],self.glroll[1][p],
				  self.ctwist[0][p],self.ctwist[1][p],self.emtylist[p]))
							
		outfile.write("\n       sum." "%17.3f %19.3f %19.3f %19.3f %19.3f %19.3f %19.3f\n" % (self.sftilt,
		              self.sfroll,self.sbpangle,self.sorient,self.sgtilt,self.sgroll,self.sctwist))
		outfile.write("       ave." "%17.3f %19.3f %19.3f %19.3f %19.3f %19.3f %19.3f\n" % (self.avftilt,
			      self.avfroll,self.avbpangle,self.avorient,self.avgtilt,self.avgroll,self.avctwist))
		outfile.write("       S.D." "%17.3f %19.3f %19.3f %19.3f %19.3f %19.3f %19.3f\n" % (self.sdftilt,
			      self.sdfroll,self.sdbpangle,self.sdorient,self.sdgtilt,self.sdgroll,self.sdctwist))
		outfile.write("\nThe global bend-angle = %.2f +/- %.2f\n" % (self.gbend,self.stgbend))
		outfile.write("**********************************************************************************************************************************************************************\n")
		
		if self.verbose == False:
			outfile.close()
		
	def _WriteBendFile(self, files):	
		
		"""Print all data to file"""
		
		print("    * Writing global bend analysis data to file: %s" % os.path.splitext(os.path.basename(files))[0]+'.bend')
		
		if self.verbose == True:
			outfile = sys.stdout 
		else:
			outfile = file(os.path.splitext(os.path.basename(files))[0]+'.bend','w')
		
		outfile.write("**************************************************************************************************************\n")	
		outfile.write("Filename: %s\n" % os.path.basename(os.path.splitext(files)[0]+'.pdb'))
		outfile.write("Date/time: %s\n" % ctime())
		outfile.write("Reference base-pair: %.2f\n" % (self.corr_refbp+1))
		outfile.write("Zone: start %i end %i\n" % (self.corr_zone[0],self.corr_zone[1]))
		outfile.write("**************************************************************************************************************\n")	
		outfile.write("Global bend angles of individual base-pair steps and of all identified base-pair steps together:\n")
		outfile.write("frac. tilt    = contribution of tilt to the base-pair step bending angle\n")
		outfile.write("frac. roll    = contribution of roll to the base-pair step bending angle\n")
		outfile.write("bp-step angle = bending angle (deg.) at each base-pair step.\n") 
		outfile.write("orient angle  = orientation of the bend angle relative to the reference base-pair in Euler-angle\n")
		outfile.write("                space (from 0 to 360 deg.)\n")
                outfile.write("global tilt   = the value of Tilt (deg.) in the global reference frame.\n")
                outfile.write("global roll   = the value of Roll (deg.) in the global reference frame.\n")  
		outfile.write("Acc. twist    = accumulated twist values (deg.) above and below the reference base-pair.\n")
		outfile.write("x	     = If the reference plane is located between two basepairs, one base-pair is\n") 
		outfile.write("		       duplicated for the calculation (two x's). If the reference plane coincides with\n") 
		outfile.write("		       one base-pair no duplication is necessary (one x)\n")
		outfile.write("\nindex  bp-step  frac.tilt   frac.roll   bp-step angle orient angle global tilt  global roll  Acc.twist  refbp\n")
	
		for p in range(len(self.bpstep[files][1])):
			outfile.write("%2d %9s %10.3f %11.3f %12.3f %12.3f %13.3f %12.3f %12.3f %3s\n" % (p+1,self.bpstep[files][1][p],self.ftilt[p],
			              self.froll[p],self.bpangle[p],self.orient[p],self.gltilt[p],self.glroll[p],self.ctwist[p],self.emtylist[p]))
							
		outfile.write("\nsum." "%19.3f %11.3f %12.3f %12.3f %13.3f %12.3f %12.3f\n" % (self.sftilt,self.sfroll,self.sbpangle,self.sorient,self.sgtilt,self.sgroll,self.sctwist))
		outfile.write("ave." "%19.3f %11.3f %12.3f %12.3f %13.3f %12.3f %12.3f\n" % (self.avftilt,self.avfroll,self.avbpangle,self.avorient,self.avgtilt,self.avgroll,self.avctwist))
		outfile.write("S.D." "%19.3f %11.3f %12.3f %12.3f %13.3f %12.3f %12.3f\n" % (self.sdftilt,self.sdfroll,self.sdbpangle,self.sdorient,self.sdgtilt,self.sdgroll,self.sdctwist))
		outfile.write("\nThe global bend-angle     = %.3f\n" % self.gbend)
		outfile.write("**************************************************************************************************************\n")
		
		if self.verbose == False:
			outfile.close()
	
	def _MatchPairs(self,tablein=None,tableout=None,infile=None):	
		
		intnr = []
		residnr1 = []
		residnr2 = []
		resid1 = []
		resid2 = []
		
		splitter1 = re.compile('[.]')
		splitter2 = re.compile('_')
		splitter3 = re.compile('\]|\[')
		for lines in tablein:
			intnr.append(lines[0])
			tmp1 = splitter1.split(lines[2])
			tmp2 = []
			for n in tmp1:	
				if n == '':
					pass
				else:
					tmp2.append(n)		
			residnr1.append(int(splitter2.split(tmp2[1])[0]))
			try:													#Three letter NA-code
				residnr2.append(int(splitter2.split(tmp2[2])[0]))
				resid1.append(splitter3.split(tmp2[1])[1])
				resid2.append(splitter3.split(tmp2[1])[3])
			except:													#One letter NA-code	
				residnr2.append(int(splitter2.split(tmp2[4])[0]))
				resid1.append(splitter3.split(tmp2[2])[0])
			 	resid2.append(splitter3.split(tmp2[3])[0])
		#TODO Distinction between OSX and Linux version (tmp2[0][1] for OSX)	
		self.chainid = tmp2[0][0]	#extract chainid once from .out file
		tableout[infile] = []
		tableout[infile].append(intnr)
		tableout[infile].append(resid1)
		tableout[infile].append(resid2)
		tableout[infile].append(residnr1)
		tableout[infile].append(residnr2)	
	
	def _StepDevide(self):
	
		"""Devide all values for roll, tilt and twist according to the base-pair step number in full sequence"""
		
		for resid in range(len(self.basesequence[self.chainid])):
			self.roll[resid] = []
			self.tilt[resid] = []
			self.twist[resid] = []
	
		for structures in self.pairs:
			for n in range(len(self.basesequence[self.chainid])):
				if n+1 in self.pairs[structures][1]:
					self.roll[n].append(self.bpstep[structures][3][self.pairs[structures][1].index(n+1)])
					self.tilt[n].append(self.bpstep[structures][2][self.pairs[structures][1].index(n+1)])
					self.twist[n].append(self.bpstep[structures][4][self.pairs[structures][1].index(n+1)])
	
		for pairs in self.roll:
				self.freq[pairs] = len(self.roll[pairs])
		
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
	
	def _FormatFloat(self,number):
		
		if float(number) == 0.0 or float(number) == -0.0:
			return float(0.0001)
		else:
			return float(number)	
	
	def _SortTable(self,tablein=None,tableout=None,infile=None):
		
		"""Sorting tablelines from Readtable and extract colums. Temporary store in columns that 
		   append to library (self) with filename as key"""
		
		columns = []
		
		for n in range(5):
			columns.append([])
			
		for lines in tablein:
			try:
				columns[2].append(self._FormatFloat(lines[5]))
				columns[3].append(self._FormatFloat(lines[6]))
				columns[4].append(self._FormatFloat(lines[7]))
				columns[0].append(int(lines[0]))
				columns[1].append(lines[1])
			except:
				pass
		
		tableout[infile] = columns		
	
	def ReadOutfiles(self,files):
	
		"""Parses 3DNA .out file and extracts roll, tilt and twist values as well as sequence information"""
	
		rmsd = re.compile("RMSD of the bases")
		bpstep = re.compile("step       Shift     Slide      Rise      Tilt      Roll     Twist")
		
		for infile in files:
			print("    * Running global bend analysis on file: %s" % os.path.basename(infile))
			readfile = file(infile,'r')
			lines = readfile.readlines()	
			linecount = 1
			countstart = []
			for line in lines:
				line = line.strip()
				result1 = rmsd.match(line)
				result2 = bpstep.match(line)
				if result1:
					tablelines = self._ReadTable(linenr=len(countstart)+2,outfile=lines)
					self._MatchPairs(tablein=tablelines,tableout=self.pairs,infile=infile)
				elif result2:
					tablelines = self._ReadTable(linenr=len(countstart),outfile=lines)
					self._SortTable(tablein=tablelines,tableout=self.bpstep,infile=infile)	
				linecount += 1
        			countstart.append(linecount)	
		
	def ReadParfiles(self,files):
		
		"""Parses 3DNA .par files and extracts roll, tilt and twist values as well as sequence information.
		   Values in the parameter file only reflects the paired base-paires. Partly unpaired structures 
		   will result in missing base-pairs in the parameter file. The native sequence cannot be extracted
		   and therefor the bend calculation using .par files cannot be run in multianalysis mode
		"""
		
		zero = float(0.0001)
		
		for infile in files:
			print("    * Running global bend analysis on file: %s" % os.path.basename(infile))
			self.bpstep[infile] = []
			self.pairs[infile] = []
			readfile = file(infile, 'r')
			lines = readfile.readlines()
			resnr = []
			sequence = []
			tilt = []
			roll = []
			twist = []
			rescount = 1
			resnr.append(rescount)
			sequence.append((lines[3].split())[0])
			tilt.append(zero)
			roll.append(zero)
			twist.append(zero)
			rescount += 1
			for line in lines[4:]:
				line = line.strip()
				value = line.split()
				resnr.append(rescount)
				sequence.append(value[0])
				tilt.append(self._FormatFloat(value[10]))
				roll.append(self._FormatFloat(value[11]))
				twist.append(self._FormatFloat(value[12]))
				rescount += 1
				
			self.bpstep[infile].append(resnr)
			self.bpstep[infile].append(sequence)
			self.bpstep[infile].append(tilt)
			self.bpstep[infile].append(roll)
			self.bpstep[infile].append(twist)
		
	def CalcGlobalBend(self,multiana):
	
		"""Control module"""
		
		for files in self.pairs:
			self._CalculateCore(files)	
			self._WriteBendFile(files)
		
		if multiana == True:
			for files in self.pairs:
				self.pairs[files] = self._ConvertSeq2(self.pairs[files])
			self._StepDevide()
			self._MultiCalculate()
			self._WriteMultiBend()
	
	def GetBaseseq(self):
			
		"""Get the base sequence and pairing as it should be in a fully paired structure (input)
		   Calls GetSequence and NAsummery class from QueryPDB plugin. The output must be of a
		   fixed format: the chains in chainlib and sequence in pairs must can only have 2 lists
		   more than two chains are not supported. The two chains must be written as 5'->3' for 
		   the template strand (first list) and 3'->5' for the complementary strand (second list)"""
		
		master = os.path.splitext(self.bpstep.keys()[0])[0]+'.pdb'
		
		pdb = PDBeditor()
		pdb.ReadPDB(master)	
		xml = pdb.PDB2XML().xml()
		
		sequence = GetSequence()
		sequence.GetSequence(pdbxml=xml)
		
		naeval = NAsummery(pdbxml=xml,sequence=sequence.seqlib)
		naeval.Evaluate()	

		self.basechainlib = naeval.chainlib
		seq = ConvertSeq(naeval.pairs[self.chainid][0],naeval.pairs[self.chainid][1])
		self.basesequence[self.chainid] = seq.Export('step1')
			
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
