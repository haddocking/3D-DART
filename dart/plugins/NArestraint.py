#!/usr/bin/env python2.7

USAGE= """
==========================================================================================

Author:			Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center for 
			Biomolecular Research, Utrecht university, The Netherlands.
Copyright (C):		2006 (DART project)
DART version:		1.0 (01-01-2007)
DART plugin: 		NArestraint.py
Input:			Structure coordinate file in PDB format
Output:			dna-rna_restraints.def file		
Plugin excecution:	Either command line driven (use -h/--help for the option) or as 
			part of a DART batch sequence.
Plugin function:	Generates a HADDOCK dna-rna_restraint.def file from a PDB file.
					it uses 3DNA to find the watson-Crick pairs in the structure. By
					default the most commen settings are used. You can change them if
					needed have a look at the options.
Plugin dependencies:	Standard python2.3 or higher modules, 3DNA
Examples:		NABrestraint -f *.pdb	

for further information, please contact:
			- DART website (http://www.nmr.chem.uu.nl/DART)
			- email: abonvin@chem.uu.nl
			- For 3DNA: rutchem.rutgers.edu/~xiangjun/3DNA 
			
==========================================================================================
"""

"""Import modules"""
import os, sys, re, string

"""Setting pythonpath variables if run from the command line"""
base, dirs = os.path.split(os.path.dirname(os.path.join(os.getcwd(), __file__)))

if base in sys.path:
	pass
else:
	sys.path.append(base)

def PluginXML():
	
	PluginXML = """ 
<metadata>
 <name>dna-rna_rstraints.def file generator</name>
 <input type="filetype">.pdb</input>
 <output type="filetype">dna-rna_restraints.def</output>
</metadata>
<parameters>
 <option type="useplugin" form="hidden" text="None">True</option>
 <option type="inputfrom" form="hidden" text="None">8</option>
 <option type="verbose" form="checkbox" text="Verbose output">False</option>
 <option type="bpplan" form="checkbox" text="Restrain base-pair planarity">False</option>
 <option type="bplan" form="checkbox" text="Restrain base planarity">True</option>
 <option type="c1pick" form="checkbox" text="Restrain C1'-C1' virtual bond distance">False</option>
 <option type="c1lower" form="text" text="C1'-restraint lower distance limit">0.05</option>
 <option type="c1upper" form="text" text="C1'-restraint upper distance limit">0.05</option>
 <option type="pickwc" form="checkbox" text="Restrain Watson-Crick hydrogen bond length">True</option>
 <option type="wc_up" form="text" text="WC-restraint upper distance limit">0.05</option>
 <option type="wc_low" form="text" text="WC-restraint lower distance limit">0.05</option>
 <option type="wc_uri_up" form="text" text="WC-restraint (uracil) upper distance limit">0.01</option>
 <option type="wc_uri_low" form="text" text="WC-restraint (uracil) lower distance limit">0.01</option>
 <option type="pickpuck" form="checkbox" text="Restrain sugar pucker dihedrals of input structure">True</option>
 <option type="pickbackdih" form="checkbox" text="Restrain phosphate backbone dihedrals of input structure">True</option>
 <option type="pform_1" form="list" default="Bform,Aform,Other" text="Conformation of group 1">Bform</option>
 <option type="puck_1_start" form="text" text="Start residue"></option>
 <option type="puck_1_end" form="text" text="End residue"></option>
 <option type="puck_1_nu2" form="text" text="dihedral C1'-C2'-C3'-C4'">-34.9</option>
 <option type="puck_1_nu2err" form="text" text="dihedral C1'-C2'-C3'-C4' error">0.0</option>
 <option type="puck_1_nu3" form="text" text="dihedral C5'-C4'-C3'-C2'">-86.4</option>
 <option type="puck_1_nu3err" form="text" text="dihedral C5'-C4'-C3'-C2' error">0.0</option>
 <option type="puck_1_nu4" form="text" text="dihedral C1'-O4'-C4'-C5'">106.4</option>
 <option type="puck_1_nu4err" form="text" text="dihedral C1'-O4'-C4'-C5' error">0.0</option>
 <option type="pform_2" form="list" default="Bform,Aform,Other" text="Conformation of group 2">Bform</option>
 <option type="puck_2_start" form="text" text="Start residue"></option>
 <option type="puck_2_end" form="text" text="End residue"></option>
 <option type="puck_2_nu2" form="text" text="dihedral C1'-C2'-C3'-C4'">-34.9</option>
 <option type="puck_2_nu2err" form="text" text="dihedral C1'-C2'-C3'-C4' error">0.0</option>
 <option type="puck_2_nu3" form="text" text="dihedral C5'-C4'-C3'-C2'">-86.4</option>
 <option type="puck_2_nu3err" form="text" text="dihedral C5'-C4'-C3'-C2' error">0.0</option>
 <option type="puck_2_nu4" form="text" text="dihedral C1'-O4'-C4'-C5'">106.4</option>
 <option type="puck_2_nu4err" form="text" text="dihedral C1'-O4'-C4'-C5' error">0.0</option>
 <option type="pform_3" form="list" default="Bform,Aform,Other" text="Conformation of group 3">Bform</option>
 <option type="puck_3_start" form="text" text="Start residue"></option>
 <option type="puck_3_end" form="text" text="End residue"></option>
 <option type="puck_3_nu2" form="text" text="dihedral C1'-C2'-C3'-C4'">-34.9</option>
 <option type="puck_3_nu2err" form="text" text="dihedral C1'-C2'-C3'-C4' error">0.0</option>
 <option type="puck_3_nu3" form="text" text="dihedral C5'-C4'-C3'-C2'">-86.4</option>
 <option type="puck_3_nu3err" form="text" text="dihedral C5'-C4'-C3'-C2' error">0.0</option>
 <option type="puck_3_nu4" form="text" text="dihedral C1'-O4'-C4'-C5'">106.4</option>
 <option type="puck_3_nu4err" form="text" text="dihedral C1'-O4'-C4'-C5' error">0.0</option>
 <option type="dform_1" form="list" default="Bform,Aform,Other" text="Conformation of group 1">Bform</option>
 <option type="dih_1_start" form="text" text="Start residue"></option>
 <option type="dih_1_end" form="text" text="End residue"></option>
 <option type="dih_1_alpha" form="text" text="alpha dihedral O3'-P-O5'-C5'">-10.0</option>
 <option type="dih_1_alphaerr" form="text" text="alpha dihedral O3'-P-O5'-C5' error">10.0</option>
 <option type="dih_1_beta" form="text" text="beta dihedral P-O5'-C5'-C4'">136.4</option>
 <option type="dih_1_betaerr" form="text" text="beta dihedral P-O5'-C5'-C4' error">40.0</option>
 <option type="dih_1_gamma" form="text" text="gamma dihedral O5'-C5'-C4'-C3'">31.1</option>
 <option type="dih_1_gammaerr" form="text" text="gamma dihedral O5'-C5'-C4'-C3' error">20.0</option>
 <option type="dih_1_delta" form="text" text="delta dihedral C5'-C4'-C3'-O3'">-165.0</option>
 <option type="dih_1_deltaerr" form="text" text="delta dihedral C5'-C4'-C3'-O3' error">50.0</option>
 <option type="dih_1_eps" form="text" text="epsilon dihedral C4'-C3'-O3'-P">-165.0</option>
 <option type="dih_1_epserr" form="text" text="epsilon dihedral C4'-C3'-O3'-P error">10.0</option>
 <option type="dih_1_zeta" form="text" text="zeta dihedral C3'-O3'-P-O5'">-150.8</option>
 <option type="dih_1_zetaerr" form="text" text="zeta dihedral C3'-O3'-P-O5' error">50.0</option>
 <option type="dform_2" form="list" default="Bform,Aform,Other" text="Conformation of group 2">Bform</option>
 <option type="dih_2_start" form="text" text="Start residue"></option>
 <option type="dih_2_end" form="text" text="End residue"></option>
 <option type="dih_2_alpha" form="text" text="alpha dihedral O3'-P-O5'-C5'">-10.0</option>
 <option type="dih_2_alphaerr" form="text" text="alpha dihedral O3'-P-O5'-C5' error">10.0</option>
 <option type="dih_2_beta" form="text" text="beta dihedral P-O5'-C5'-C4'">136.4</option>
 <option type="dih_2_betaerr" form="text" text="beta dihedral P-O5'-C5'-C4' error">40.0</option>
 <option type="dih_2_gamma" form="text" text="gamma dihedral O5'-C5'-C4'-C3'">31.1</option>
 <option type="dih_2_gammaerr" form="text" text="gamma dihedral O5'-C5'-C4'-C3' error">20.0</option>
 <option type="dih_2_delta" form="text" text="delta dihedral C5'-C4'-C3'-O3'">-165.0</option>
 <option type="dih_2_deltaerr" form="text" text="delta dihedral C5'-C4'-C3'-O3' error">50.0</option>
 <option type="dih_2_eps" form="text" text="epsilon dihedral C4'-C3'-O3'-P">-165.0</option>
 <option type="dih_2_epserr" form="text" text="epsilon dihedral C4'-C3'-O3'-P error">10.0</option>
 <option type="dih_2_zeta" form="text" text="zeta dihedral C3'-O3'-P-O5'">-150.8</option>
 <option type="dih_2_zetaerr" form="text" text="zeta dihedral C3'-O3'-P-O5' error">50.0</option>
 <option type="dform_3" form="list" default="Bform,Aform,Other" text="Conformation of group 3">Bform</option>
 <option type="dih_3_start" form="text" text="Start residue"></option>
 <option type="dih_3_end" form="text" text="End residue"></option>
 <option type="dih_3_alpha" form="text" text="alpha dihedral O3'-P-O5'-C5'">-10.0</option>
 <option type="dih_3_alphaerr" form="text" text="alpha dihedral O3'-P-O5'-C5' error">10.0</option>
 <option type="dih_3_beta" form="text" text="beta dihedral P-O5'-C5'-C4'">136.4</option>
 <option type="dih_3_betaerr" form="text" text="beta dihedral P-O5'-C5'-C4' error">40.0</option>
 <option type="dih_3_gamma" form="text" text="gamma dihedral O5'-C5'-C4'-C3'">31.1</option>
 <option type="dih_3_gammaerr" form="text" text="gamma dihedral O5'-C5'-C4'-C3' error">20.0</option>
 <option type="dih_3_delta" form="text" text="delta dihedral C5'-C4'-C3'-O3'">-165.0</option>
 <option type="dih_3_deltaerr" form="text" text="delta dihedral C5'-C4'-C3'-O3' error">50.0</option>
 <option type="dih_3_eps" form="text" text="epsilon dihedral C4'-C3'-O3'-P">-165.0</option>
 <option type="dih_3_epserr" form="text" text="epsilon dihedral C4'-C3'-O3'-P error">10.0</option>
 <option type="dih_3_zeta" form="text" text="zeta dihedral C3'-O3'-P-O5'">-150.8</option>
 <option type="dih_3_zetaerr" form="text" text="zeta dihedral C3'-O3'-P-O5' error">50.0</option>
</parameters>"""
	
	return PluginXML

def PluginCore(paramdict, inputlist):
	
	print "--> Generating dna-rna_restraints.def file for use in HADDOCK"
	print "    * Check input parameters"
	
	if not paramdict['pickpuck']:
		valid = False
		for puck in [1,2,3]:
			if paramdict['puck_%i_start' % puck] and paramdict['puck_%i_end' % puck]: valid  = True
		if not valid:
			print "    * WARNING: automatic restraint definition for sugar pucker dihedral angles turned off"
			print "               but no conformational groups defined. Turn automatic defintion on."
			paramdict['pickpuck'] = True
	
	if not paramdict['pickbackdih']:
		valid = False
		for dih in [1,2,3]:
			if paramdict['dih_%i_start' % dih] and paramdict['dih_%i_end' % dih]: valid  = True
		if not valid:
			print "    * WARNING: automatic restraint definition for phosphate backbone dihedral angles turned off"
			print "               but no conformational groups defined. Turn automatic defintion on."
			paramdict['pickbackdih'] = True		
	
	restraints = NArestraints(paramdict)
	restraints.readpdb(inputlist[0])
	restraints.writedef()

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

		parser.add_option( "-f", "--file", action="callback", callback=self.varargs, dest="inputfile", type="string", help="Supply pdb. Standard UNIX selection syntax accepted")
		parser.add_option( "-a", "--bpplan", action="store_true", dest="bpplan", default=False, help="Restrain base-pair planarity, default=False")
		parser.add_option( "-b", "--bplan", action="store_true", dest="bplan", default=True, help="Restrain base planarity, default=True")
		parser.add_option( "-c", "--pickpuc", action="store_true", dest="pickpuc", default=True, help="Restrain sugar pucker dihedrals of input structure, default=True")
		parser.add_option( "-d", "--pickbacdih", action="store_true", dest="pickbacdih", default=True, help="Restrain phosphate backbone dihedrals of input structure, default=True")
		parser.add_option( "-e", "--c1pick", action="store_true", dest="c1pick", default=False, help="Restrain C1'-C1' virtual bond distance, default=False")
		parser.add_option( "-g", "--c1lower", action="store", dest="c1lower", default=0.05, help="C1'-restraint upper distance limit, default=0.05")
		parser.add_option( "-i", "--c1upper", action="store", dest="c1upper", default=0.05, help="C1'-restraint lower distance limit, default=0.05")
		parser.add_option( "-j", "--wcpairing", action="store_true", dest="wcpairing", default=True, help="Restrain Watson-Crick hydrogen bond length, default=True")
		parser.add_option( "-k", "--wc_up", action="store", dest="wc_up", default=0.05, help="WC-restraint lower distance limit, default=0.05")
		parser.add_option( "-l", "--wc_low", action="store", dest="wc_low", default=0.05, help="WC-restraint lower distance limit, default=0.05")
		parser.add_option( "-m", "--wc_uri_up", action="store", dest="wc_uri_up", default=0.01, help="WC-restraint (uracil) upper distance limit, default=0.01")
		parser.add_option( "-n", "--wc_uri_low", action="store", dest="wc_uri_low", default=0.01, help="WC-restraint (uracil) lower distance limit, default=0.01")		
		parser.add_option( "-v", "--verbose", action="store_true", dest="verbose", default=False, help="All output to standard output")
		
		(options, args) = parser.parse_args()
		
		self.option_dict['input'] = options.inputfile
		self.option_dict['verbose'] = options.verbose
		self.option_dict['bpplan'] = options.bpplan
		self.option_dict['bplan'] = options.bplan
		self.option_dict['pickpuc'] = options.pickpuc
		self.option_dict['pickbacdih'] = options.pickbacdih
		self.option_dict['c1pick'] = options.c1pick
		self.option_dict['c1lower'] = options.c1lower
		self.option_dict['c1upper'] = options.c1upper	
		self.option_dict['wcpairing'] = options.wcpairing	
		self.option_dict['wc_up'] = options.wc_up
		self.option_dict['wc_low'] = options.wc_low
		self.option_dict['wc_uri_up'] = options.wc_uri_up
		self.option_dict['wc_uri_low'] = options.wc_uri_low	
			
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

class NArestraints:
	
	def __init__(self, paramdict=None):
		
		self.paramdict = paramdict
		
		self.segid1 = []
		self.segid2 = []
		self.resid1 = []
		self.resid2 = []
		self.resnr1 = []
		self.resnr2 = []
		self.seglib = {}
	
	def getzones(self):

		self.makelib()

		while not len(self.resnr1) == 0:
			tmp,new,segid = self.loop(inlist=self.resnr1,insegid=self.segid1)
			self.seglib[segid].append(tmp)
			for n in tmp:
				if n in self.resnr1:
					self.resnr1.remove(n)
			self.segid1 = new

		while not len(self.resnr2) == 0:
			tmp,new,segid = self.loop(inlist=self.resnr2,insegid=self.segid2)
			self.seglib[segid].append(tmp)
			for n in tmp:
				if n in self.resnr2:
					self.resnr2.remove(n)
			self.segid2 = new
			
	def loop(self,inlist,insegid):

		segid = insegid[0]
		cut = 0
		tmp = []
		if inlist[1] == inlist[0]+1:
			sign = 1
		elif inlist[1] == inlist[0]-1:
			sign = -1
		for n in range(len(inlist)):
			try:
				if inlist[n+1] == inlist[n]+sign and insegid[n+1] == segid:
					tmp.append(inlist[n])
					cut += 1
				else:
					tmp.append(inlist[n])
					cut +=1
					break
			except:
					tmp.append(inlist[n])
					cut +=1

		new = []
		for n in range(cut,len(insegid)):
			new.append(insegid[n])

		return tmp,new,segid				

	def makelib(self):

		for segid in self.segid1:
			if not self.seglib.has_key(segid): self.seglib[segid] = []
		for segid in self.segid2:
			if not self.seglib.has_key(segid): self.seglib[segid] = []
	
	def readpdb(self,pdb):
	
		"""Run 3DNA find_pair on pdb file. Only generate .inp file and pass to 'importinp'"""
	
		os.system("find_pair -t %s output.inp" % pdb)
	
		if os.path.isfile('output.inp'): self.importinp('output.inp')
		else:
			print "--> No 3DNA input file to import, stopping"
			sys.exit(0)
		
		clean = ['output.inp','col_helices.scr','hel_regions.pdb','col_chains.scr','bp_order.dat','bestpairs.pdb','ref_frames.dat']
		
		for files in clean:
			if os.path.isfile(files): os.remove(files)
	
	def importinp(self,inpfile):
		
		readfile = file(inpfile,'r')
		lines = readfile.readlines()
		for line in lines:
			if len(line) == 87:
				self.segid1.append(line[23])
				self.resnr1.append(int(line[25:29].strip('.')))
				self.resid1.append(line[32:35])
				self.resid2.append(line[44:47])
				self.resnr2.append(int(line[49:53].strip('.')))
				self.segid2.append(line[55])	
		
		self.paramdict['pairs'] = []
		for n in range(len(self.segid1)):
			self.paramdict['pairs'].append((self.segid1[n],self.resid1[n],self.resnr1[n],self.segid2[n],self.resid2[n],self.resnr2[n]))
	
		self.getzones()
		
	def writedef(self):
	
		if self.paramdict['verbose']: outfile = sys.stdout
		else: outfile = file('dna-rna_restraints.def', 'w')
			
		self.header(outfile)
		self.bpplanarity(outfile)
		self.baseplanarity(outfile)
		self.pucker(outfile)
		self.sfbackbone(outfile)
		self.c1c1restraint(outfile)
		self.wcpairing(outfile)
		self.footer(outfile)
		
		if not self.paramdict['verbose']: outfile.close()	
	
	def bpplanarity(self,outfile):
	
		outfile.write("{=========================================== base-pair planarity ===========================================}\n")
		outfile.write("{* Use planarity restraints for Watson-Crick base pairing *}\n")
		outfile.write("{+ choice: true false +}\n\n")
		
		outfile.write("{===>} basepair_planar=%s;\n\n" % string.lower(str(self.paramdict['bpplan'])))
		
	def baseplanarity(self, outfile):	
		
		outfile.write ("{============================================== base planarity =============================================}\n\n")
		outfile.write ("{* Restrain base planarity. This selection must only include nucleotide residues *}\n\n")
		
		zone = ""
		for segid in self.seglib:
			resids = ""
			for n in self.seglib[segid]:
				if len(resids) == 0:
					resids = "resid "+str(n[0])+":"+str(n[-1])
				else:
					resids = resids+" or resid "+str(n[0])+":"+str(n[-1])
			if len(zone) == 0:
				zone = "("+resids+") and segid "+segid
			else:
				zone = zone+" or ("+resids+") and segid "+segid
		
		if self.paramdict['bplan'] == True:	
			outfile.write ("{===>} bases_planar=(%s);\n\n" % zone)	
		else:
			outfile.write ("{* Base planarity not restraint *}\n\n")
		
	def pucker(self,outfile):
	
		outfile.write("{=================================== sugar-pucker dihedral angle restraints ================================}\n\n")

		outfile.write("{* Pick the dihedral angles of the sugar pucker from the input structure\n")
		outfile.write("   and restrain them within the given error range *}\n")
		outfile.write("{+ choice: true false +}\n\n")
		
		if self.paramdict['pickpuck']:
			outfile.write("{===>} dna_pick_pucdih=true;\n")
	
			pucker = []
			for segid in self.seglib:
				for n in self.seglib[segid]:
					zone = "resid "+str(n[0])+":"+str(n[-1])+" and segid "+segid
					pucker.append(zone)
	    	
			puckercount = 1
			for pucker_group in pucker:
			
				outfile.write("{* residues with sugar pucker restrained - group %i *}\n" % puckercount)
				outfile.write("{===>} pucker_%i=(%s);\n\n" % (puckercount,pucker_group))
			
				outfile.write("{* conformation of group %i *}\n" % puckercount)
				outfile.write('{+ choice: "a-form" "b-form" "other" +}\n')
				outfile.write('{===>} pform_%i="other";\n\n' % puckercount)

				outfile.write("{* user defined sugar pucker for group %i *}\n" % puckercount)
			
				outfile.write("{* dihedral C1'-C2'-C3'-C4' *}\n")
				outfile.write("{===>} dihedral_nu2_%i=-34.9;\n" % puckercount) 
				outfile.write("{* dihedral C1'-C2'-C3'-C4' error range *}\n")
				outfile.write("{===>} error_nu2_%i=0.0;\n" % puckercount)
				outfile.write("{* dihedral C5'-C4'-C3'-C2' *}\n")
				outfile.write("{===>} dihedral_nu3_%i=-86.4;\n" % puckercount) 
				outfile.write("{* dihedral C5'-C4'-C3'-C2' error range *}\n")
				outfile.write("{===>} error_nu3_%i=0.0;\n" % puckercount)
				outfile.write("{* dihedral C1'-O4'-C4'-C5' *}\n")
				outfile.write("{===>} dihedral_nu4_%i=106.4;\n" % puckercount) 
				outfile.write("{* dihedral C1'-O4'-C4'-C5' error range *}\n")
				outfile.write("{===>} error_nu4_%i=0.0;\n\n" % puckercount)
	
				puckercount += 1
		else:
			outfile.write("{===>} dna_pick_pucdih=false;\n")
			
			for puck in [1,2,3]:
				if self.paramdict['puck_%i_start' % puck] and self.paramdict['puck_%i_end' % puck]:
					outfile.write("{* residues with sugar pucker restrained - group %i *}\n" % puck)
					
					segid = None
					for n in self.seglib:
						if self.paramdict['puck_%i_start' % puck] in self.seglib[n][0] and self.paramdict['puck_%i_end' % puck] in self.seglib[n][0]: segid = n
					
					outfile.write("{===>} pucker_%i=(resid %i:%i and segid %s);\n\n" % 
					             (puck,int(self.paramdict['puck_%i_start' % puck]),int(self.paramdict['puck_%i_end' % puck]), segid))

					outfile.write("{* conformation of group %i *}\n" % puck)
					outfile.write('{+ choice: "a-form" "b-form" "other" +}\n')
					outfile.write("{===>} pform_%i=\"%s\";\n\n" % (puck,self.paramdict['pform_%i' % puck]))

					outfile.write("{* user defined sugar pucker for group %i *}\n" % puck)

					outfile.write("{* dihedral C1'-C2'-C3'-C4' *}\n")
					outfile.write("{===>} dihedral_nu2_%i=%s;\n" % (puck,self.paramdict['puck_%i_nu2' % puck])) 
					outfile.write("{* dihedral C1'-C2'-C3'-C4' error range *}\n")
					outfile.write("{===>} error_nu2_%i=%s;\n" % (puck,self.paramdict['puck_%i_nu2err' % puck])) 
					outfile.write("{* dihedral C5'-C4'-C3'-C2' *}\n")
					outfile.write("{===>} dihedral_nu3_%i=%s;\n" % (puck,self.paramdict['puck_%i_nu3' % puck])) 
					outfile.write("{* dihedral C5'-C4'-C3'-C2' error range *}\n")
					outfile.write("{===>} error_nu3_%i=%s;\n" % (puck,self.paramdict['puck_%i_nu3err' % puck])) 
					outfile.write("{* dihedral C1'-O4'-C4'-C5' *}\n")
					outfile.write("{===>} dihedral_nu4_%i=%s;\n" % (puck,self.paramdict['puck_%i_nu4' % puck])) 
					outfile.write("{* dihedral C1'-O4'-C4'-C5' error range *}\n")
					outfile.write("{===>} error_nu4_%i=%s;\n\n" % (puck,self.paramdict['puck_%i_nu4err' % puck])) 
	
	def sfbackbone(self,outfile):
		
		outfile.write("{================================ phosphate backbone dihedral angle restraints =============================}\n\n")
		
		outfile.write("{* Pick the dihedral angles of the phosphate backbone from the input structure and\n")
		outfile.write("   restrain them within the given error range *}\n")
		outfile.write("{+ choice: true false +}\n\n")
		
		if self.paramdict['pickbackdih']: 
	
			outfile.write("{===>} dna_pick_bacdih=true;\n\n")
		
			bacdih = []
			for segid in self.seglib:
				for n in self.seglib[segid]:
					zone = "resid "+str(n[0])+":"+str(n[-1])+" and segid "+segid
					bacdih.append(zone)
		
			bacdihcount = 1
			for bacdih_group in bacdih:
			
				outfile.write("{* residues with phosphate backbone restrained - group %i *}\n" % bacdihcount)
				outfile.write("{===>} dihedral_%i=(%s);\n\n" % (bacdihcount,bacdih_group))
			
				outfile.write("{* conformation of group %i *}\n" % bacdihcount)
				outfile.write('{+ choice: "a-form" "b-form" "other" +}\n')
				outfile.write('{===>} dform_%i="other";\n\n' % bacdihcount)

				outfile.write("{* user defined posphate backbone for group %i *}\n" % bacdihcount)
			
				outfile.write("{* alpha dihedral O3'-P-O5'-C5' *}\n")
				outfile.write("{===>} dihedral_alpha_%i=-10.0;\n" % bacdihcount) 
				outfile.write("{* alpha dihedral range *}\n")
				outfile.write("{===>} error_alpha_%i=10.0;\n" % bacdihcount) 
				outfile.write("{* beta dihedral P-O5'-C5'-C4' *}\n")
				outfile.write("{===>} dihedral_beta_%i=136.4;\n" % bacdihcount) 
				outfile.write("{* beta dihedral range *}\n")
				outfile.write("{===>} error_beta_%i=40.0;\n" % bacdihcount) 
				outfile.write("{* gamma dihedral O5'-C5'-C4'-C3' *}\n")
				outfile.write("{===>} dihedral_gamma_%i=31.1;\n" % bacdihcount) 
				outfile.write("{* gamma dihedral range *}\n")
				outfile.write("{===>} error_gamma_%i=20.0;\n" % bacdihcount) 
				outfile.write("{* delta dihedral C5'-C4'-C3'-O3' *}\n")
				outfile.write("{===>} dihedral_delta_%i=-165.0;\n" % bacdihcount) 
				outfile.write("{* delta dihedral range *}\n")
				outfile.write("{===>} error_delta_%i=50.0;\n" % bacdihcount) 
				outfile.write("{* epsilon dihedral C4'-C3'-O3'-P *}\n")
				outfile.write("{===>} dihedral_eps_%i=-165.0;\n" % bacdihcount) 
				outfile.write("{* epsilon dihedral range *}\n")
				outfile.write("{===>} error_eps_%i=10.0;\n" % bacdihcount) 
				outfile.write("{* zeta dihedral C3'-O3'-P-O5' *}\n")
				outfile.write("{===>} dihedral_zeta_%i=-150.8;\n" % bacdihcount) 
				outfile.write("{* zeta dihedral range *}\n")
				outfile.write("{===>} error_zeta_%i=50.0;\n\n" % bacdihcount)
		
				bacdihcount += 1
		else:
			outfile.write("{===>} dna_pick_bacdih=false;\n\n")

			for dih in [1,2,3]:
				if self.paramdict['dih_%i_start' % dih] and self.paramdict['dih_%i_end' % dih]:
					outfile.write("{* residues with phosphate backbone restrained - group %i *}\n" % dih)

					segid = None
					for n in self.seglib:
						if self.paramdict['dih_%i_start' % dih] in self.seglib[n][0] and self.paramdict['dih_%i_end' % dih] in self.seglib[n][0]: segid = n

					outfile.write("{===>} dihedral_%i=(resid %i:%i and segid %s);\n\n" % 
					             (dih,int(self.paramdict['dih_%i_start' % dih]),int(self.paramdict['dih_%i_end' % dih]), segid))
					
					outfile.write("{* conformation of group %i *}\n" % dih)
					outfile.write('{+ choice: "a-form" "b-form" "other" +}\n')
					outfile.write('{===>} dform_%i="%s";\n\n' % (dih,self.paramdict['dform_%i' % dih]))

					outfile.write("{* user defined posphate backbone for group %i *}\n" % dih)

					outfile.write("{* alpha dihedral O3'-P-O5'-C5' *}\n")
					outfile.write("{===>} dihedral_alpha_%i=%s;\n" % (dih,self.paramdict['dih_%i_alpha' % dih])) 
					outfile.write("{* alpha dihedral range *}\n")
					outfile.write("{===>} error_alpha_%i=%s;\n" % (dih,self.paramdict['dih_%i_alphaerr' % dih]))
					outfile.write("{* beta dihedral P-O5'-C5'-C4' *}\n")
					outfile.write("{===>} dihedral_beta_%i=%s;\n" % (dih,self.paramdict['dih_%i_beta' % dih])) 
					outfile.write("{* beta dihedral range *}\n")
					outfile.write("{===>} error_beta_%i=%s;\n" % (dih,self.paramdict['dih_%i_betaerr' % dih])) 
					outfile.write("{* gamma dihedral O5'-C5'-C4'-C3' *}\n")
					outfile.write("{===>} dihedral_gamma_%i=%s;\n" % (dih,self.paramdict['dih_%i_gamma' % dih])) 
					outfile.write("{* gamma dihedral range *}\n")
					outfile.write("{===>} error_gamma_%i=%s;\n" % (dih,self.paramdict['dih_%i_gammaerr' % dih])) 
					outfile.write("{* delta dihedral C5'-C4'-C3'-O3' *}\n")
					outfile.write("{===>} dihedral_delta_%i=%s;\n" % (dih,self.paramdict['dih_%i_delta' % dih])) 
					outfile.write("{* delta dihedral range *}\n")
					outfile.write("{===>} error_delta_%i=%s;\n" % (dih,self.paramdict['dih_%i_deltaerr' % dih])) 
					outfile.write("{* epsilon dihedral C4'-C3'-O3'-P *}\n")
					outfile.write("{===>} dihedral_eps_%i=%s;\n" % (dih,self.paramdict['dih_%i_eps' % dih])) 
					outfile.write("{* epsilon dihedral range *}\n")
					outfile.write("{===>} error_eps_%i=%s;\n" % (dih,self.paramdict['dih_%i_epserr' % dih])) 
					outfile.write("{* zeta dihedral C3'-O3'-P-O5' *}\n")
					outfile.write("{===>} dihedral_zeta_%i=%s;\n" % (dih,self.paramdict['dih_%i_zeta' % dih])) 
					outfile.write("{* zeta dihedral range *}\n")
					outfile.write("{===>} error_zeta_%i=%s;\n\n" % (dih,self.paramdict['dih_%i_zetaerr' % dih]))
				
	def c1c1restraint(self,outfile):
		
		outfile.write("{============================================= C1'-C1' restraints ==========================================}\n\n")
		outfile.write("{* Have the length of the C1'-C1' virtual bonds measured and restraints. *}\n")
		outfile.write("{+ choice: true false +}\n")
		
		outfile.write("{===>} dna_pick_c1=%s;\n\n" % string.lower(str(self.paramdict['c1pick'])))
		
		outfile.write("{* Error range used for C1'-C1' virtual bonds  *}\n")
		outfile.write("{===>} c1_low=%1.3f;\n" % self.paramdict['c1lower'])
		outfile.write("{===>} c1_up=%1.3f;\n\n" % self.paramdict['c1upper'])	
	
	def wcpairing(self,outfile):
		
		outfile.write("{=========================================== Watson-Crick base pairs =======================================}\n\n")

		outfile.write("{* pick Watson-Crick restraint values from structure *}\n")
		outfile.write("{+ choice: true false +}\n")
		
		outfile.write("{===>} dna_pick_wc=%s;\n" % string.lower(str(self.paramdict['pickwc'])))
		outfile.write("{* error range used for dna_pick_wc defined Watson-Crick restraints *}\n")
		outfile.write("{===>} wc_low=%1.3f;\n" % self.paramdict['wc_low'])
		outfile.write("{===>} wc_up=%1.3f;\n" % self.paramdict['wc_up'])
		outfile.write("{* for URI, for default much lower range... why?*}\n")
		outfile.write("{===>} wc_low_uri=%1.3f;\n" % self.paramdict['wc_uri_low'])
		outfile.write("{===>} wc_up_uri=%1.3f;\n\n" % self.paramdict['wc_uri_up'])

		outfile.write("{* residues which form Watson-Crick pairs *}\n\n")

		paircount = 1
		for basepair in self.paramdict['pairs']:
			
			outfile.write("{* selection for pair %i base A *}\n" % paircount)
			outfile.write("{===>} base_a_%i=(resid %i and segid %s);\n" % (paircount,basepair[2],basepair[0]))
			outfile.write("{* selection for pair %i base B *}\n" % paircount)
			outfile.write("{===>} base_b_%i=(resid %i and segid %s);\n\n" % (paircount,basepair[5],basepair[3]))
		
			paircount += 1
		
	def header(self,outfile):
		
		outfile.write("""{+ file: dna-rna_restraints.def       directory: protocols +}
{+ description: Creates restraints to maintain conformation of DNA/RNA +}
{+ comment:This file is to be read by refinement files that modify atom coordinates +}
{+ authors: Axel T. Brunger, and Paul D. Adams, <br>
            modified by Alexandre Bonvin and Marc van Dijk for HADDOCK use
            <br><br>
Additions and changes were made to allow for flexibility during docking<br><br>
Changes include: <br>
<ul>
<li> flags to turn all options on or off
<li> separation of sugar
<li> pucker restraints and phosphate backbone restraints
<li> option to have sugar-phosphate backbone dihedrals measured and restrained within a user defined error range
<li> option to have the length of the Watson-Crick hydrogen bonds measured from the structure measured and restrained within a user defined error range.
</ul>
+}

set message=normal echo=on end

{- begin block parameter definition -} define(\n\n""")
		
	def footer(self,outfile):
		
		outfile.write("""{=========================================================================================================}
{                        things below this line do not normally need to be changed                        }
{=========================================================================================================}

 ) {- end block parameter definition -}

{- the planarity restraints for Watson-Crick base pairing -}

if (&basepair_planar=true) then
 evaluate ($pair=1)
 evaluate ($done=false)
 while ( $done = false ) loop plan_paired
   if ( &exist_base_a_$pair = true ) then
     if ( &exist_base_b_$pair = true ) then
       show (segid) ( &base_a_$pair and name C1' ) 
       evaluate ($Asegid=$result)
       show (resid) ( &base_a_$pair and name C1' ) 
       evaluate ($Aresid=$result)
       show (segid) ( &base_b_$pair and name C1' ) 
       evaluate ($Bsegid=$result)
       show (resid) ( &base_b_$pair and name C1' ) 
       evaluate ($Bresid=$result)
       evaluate ($plweight = 20)			! Enforce planarity by increasing plweight value.

       restraints plane

         group
           selection=(((segid $Asegid and resid $Aresid) or (segid $Bsegid and resid $Bresid)) and
                      (resname THY or resname CYT or resname GUA or
                       resname ADE or resname URI) and
                       not (name c#' or name h#' or name h#'' or name o#p or
                            name h7# or name o#' or name p or name h#t or name o#t))
           weight=$plweight
         end
       end
     end if
   else
     evalute ($done = true)
   end if
     evaluate ($pair = $pair + 1)
 end loop plan_paired
else
end if
flag include plan end

{- the planarity restraints single bases -}

 for $id in id ( &bases_planar and tag ) loop plan
   show (segid) (id $id)
   evaluate ($segid=$result)
   show (resid) (id $id)
   evaluate ($resid=decode($result))
   evaluate ($plweight = 20)

   restraints plane

     group
       selection=( segid $segid and resid $resid and
                  (resname THY or resname CYT or resname GUA or
                   resname ADE or resname URI) and
                   not (name c#' or name h#' or name h#'' or name o#p or
                        name h7# or name o#' or name p or name h#t or name o#t))
       weight=$plweight
     end
   end
 end loop plan

{- Dihedral restraints for the sugar pucker -}

if (&dna_pick_pucdih=true) then
  evaluate ($group=1)
  evaluate ($done=false)
  while ( $done = false ) loop dihe
   if ( &exist_pucker_$group = true ) then
     show sum(1) ( &pucker_$group )
     if ( $result > 0 ) then
       evaluate ($min_resid_$group = 99999)
       evaluate ($max_resid_$group = -99999)
       evaluate ($error_nu2=&error_nu2_$group)
       evaluate ($error_nu3=&error_nu3_$group)
       evaluate ($error_nu4=&error_nu4_$group)	 
       for $id in id ( &pucker_$group and tag ) loop resid
         show (segid) (id $id)
         evaluate ($segid=$result)
         show (resid) ( id $id )
         evaluate ($resid=decode($result))
	 evaluate ($min_resid_$group = max($min_resid_$group,$resid))
	 evaluate ($max_resid_$group = max($max_resid_$group,$resid))
         pick dihedral
                   ( segid $segid and resid $resid and name c1' )
                   ( segid $segid and resid $resid and name c2' )
                   ( segid $segid and resid $resid and name c3' )
                   ( segid $segid and resid $resid and name c4' ) 
		  geometry
	 evaluatate ($dihedral_nu2=$result)
	 pick dihedral
                   ( segid $segid and resid $resid and name c5' )
                   ( segid $segid and resid $resid and name c4' )
                   ( segid $segid and resid $resid and name c3' )
                   ( segid $segid and resid $resid and name c2' ) 
		  geometry
	 evaluatate ($dihedral_nu3=$result)
	 pick dihedral
                   ( segid $segid and resid $resid and name c1' )
                   ( segid $segid and resid $resid and name o4' )
                   ( segid $segid and resid $resid and name c4' )
                   ( segid $segid and resid $resid and name c5' ) 
		  geometry
	 evaluatate ($dihedral_nu4=$result)

      restraints dihedral
           assign  ( segid $segid and resid $resid and name c1' )
                   ( segid $segid and resid $resid and name c2' )
                   ( segid $segid and resid $resid and name c3' )
                   ( segid $segid and resid $resid and name c4' ) 
                                                       20.0 $dihedral_nu2 $error_nu2 2
           assign  ( segid $segid and resid $resid and name c5' )
                   ( segid $segid and resid $resid and name c4' )
                   ( segid $segid and resid $resid and name c3' )
                   ( segid $segid and resid $resid and name c2' ) 
                                                       20.0 $dihedral_nu3 $error_nu3 2
           assign  ( segid $segid and resid $resid and name c1' )
                   ( segid $segid and resid $resid and name o4' )
                   ( segid $segid and resid $resid and name c4' )
                   ( segid $segid and resid $resid and name c5' ) 
                                                       20.0 $dihedral_nu4 $error_nu4 2
           scale=20.0
         end
       end loop resid
     end if
   else
     evaluate ($done=true)
   end if
   evaluate ($group=$group+1)
  end loop dihe

 else

 evaluate ($group=1)
 evaluate ($done=false)
 while ( $done = false ) loop dihe
   if ( &exist_pucker_$group = true ) then
     show sum(1) ( &pucker_$group )
     if ( $result > 0 ) then
       if ( &pform_$group = "a-form" ) then
         evaluate ($dihedral_nu2=37.053)
         evaluate ($dihedral_nu3=-155.59)
         evaluate ($dihedral_nu4=144.26)
       elseif ( &pform_$group = "b-form" ) then
         evaluate ($dihedral_nu2=-34.9)
         evaluate ($dihedral_nu3=-86.4)
         evaluate ($dihedral_nu4=106.4)
       elseif ( &pform_$group = "other" ) then
         evaluate ($dihedral_nu2=&dihedral_nu2_$group)
         evaluate ($dihedral_nu3=&dihedral_nu3_$group)
         evaluate ($dihedral_nu4=&dihedral_nu4_$group)
       end if

       evaluate ($min_resid_$group = 99999)
       evaluate ($max_resid_$group = -99999)

       for $id in id ( &pucker_$group and tag ) loop resid

         show (segid) (id $id)
         evaluate ($segid=$result)
         show (resid) ( id $id )
         evaluate ($resid=decode($result))
	 evaluate ($min_resid_$group = max($min_resid_$group,$resid))
	 evaluate ($max_resid_$group = max($max_resid_$group,$resid))

         restraints dihedral
           assign  ( segid $segid and resid $resid and name c1' )
                   ( segid $segid and resid $resid and name c2' )
                   ( segid $segid and resid $resid and name c3' )
                   ( segid $segid and resid $resid and name c4' ) 
                                                       20.0 $dihedral_nu2 0.0 2
           assign  ( segid $segid and resid $resid and name c5' )
                   ( segid $segid and resid $resid and name c4' )
                   ( segid $segid and resid $resid and name c3' )
                   ( segid $segid and resid $resid and name c2' ) 
                                                       20.0 $dihedral_nu3 0.0 2
           assign  ( segid $segid and resid $resid and name c1' )
                   ( segid $segid and resid $resid and name o4' )
                   ( segid $segid and resid $resid and name c4' )
                   ( segid $segid and resid $resid and name c5' ) 
                                                       20.0 $dihedral_nu4 0.0 2

           scale=20.0
         end
       end loop resid
     end if
   else
     evaluate ($done=true)
   end if
   evaluate ($group=$group+1)
 end loop dihe
end if
flags include cdih end

{- Dihedral restraints for the phosphate backbone -}

if (&dna_pick_bacdih=true) then
  evaluate ($group=1)
  evaluate ($done=false)
  while ( $done = false ) loop bdihe
   if ( &exist_dihedral_$group = true ) then
     show sum(1) ( &dihedral_$group )
     if ( $result > 0 ) then
       evaluate ($resid=$min_resid_$group)
       evaluate ($nres=$max_resid_$group - $min_resid_$group + 1)
       evaluate ($error_alpha=&error_alpha_$group)
       evaluate ($error_beta=&error_beta_$group)
       evaluate ($error_gamma=&error_gamma_$group)
       evaluate ($error_zeta=&error_zeta_$group)
       evaluate ($error_epsilon=&error_eps_$group)
       evaluate ($error_delta=&error_delta_$group)
       for $id in id ( &dihedral_$group and tag ) loop resid
         show (segid) (id $id)
         evaluate ($segid=$result)
         show (resid) ( id $id )
         evaluate ($resid=decode($result))
         if ($resid > $min_resid_$group) then
           evaluate ($rprec = $resid - 1)
	   pick dihedral
                     ( segid $segid and resid $rprec and name O3' )
                     ( segid $segid and resid $resid and name P )
                     ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' ) 
		  geometry
	   evaluatate ($dihedral_alpha=$result)
	   pick dihedral
                     ( segid $segid and resid $resid and name P )
                     ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' ) 
		  geometry
	   evaluatate ($dihedral_beta=$result)

           restraint dihedral
	    ! alpha
             assign  ( segid $segid and resid $rprec and name O3' )
                     ( segid $segid and resid $resid and name P )
                     ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' ) 
                                                       1.0 $dihedral_alpha $error_alpha 2
	    ! beta					       
             assign  ( segid $segid and resid $resid and name P )
                     ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' ) 
                                                       1.0 $dihedral_beta $error_beta 2
             scale 200.0
           end
         end if

	 pick dihedral
                     ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' ) 
	  	  geometry
	 evaluatate ($dihedral_gamma=$result)
         pick dihedral
                     ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' ) 
		  geometry
	   evaluatate ($dihedral_delta=$result)

	 restraints dihedral
	    ! gamma
             assign  ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' ) 
                                                       1.0 $dihedral_gamma $error_gamma 2
	    ! delta					       
             assign  ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' ) 
                                                       1.0 $dihedral_delta $error_delta 2		    
	      scale=200.0
            end

	  if ($resid < $max_resid_$group) then
           evaluate ($rfoll = $resid + 1)
	   pick dihedral
                     ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' )
                     ( segid $segid and resid $rfoll and name P ) 
		  geometry
	   evaluatate ($dihedral_epsilon=$result)
	   pick dihedral
                     ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' )
                     ( segid $segid and resid $rfoll and name P )
                     ( segid $segid and resid $rfoll and name O5' ) 
		  geometry
	   evaluatate ($dihedral_zeta=$result)
           restraint dihedral
             ! epsilon
	     assign  ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' )
                     ( segid $segid and resid $rfoll and name P ) 
                                                       1.0 $dihedral_epsilon $error_epsilon 2
             ! zeta
	     assign  ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' )
                     ( segid $segid and resid $rfoll and name P )
                     ( segid $segid and resid $rfoll and name O5' ) 
                                                       1.0 $dihedral_zeta $error_zeta 2
             scale 200.0
           end
         end if
       end loop resid
     end if
   else
     evaluate ($done=true)
   end if
     evaluate ($group=$group+1)
 end loop bdihe

 else

 evaluate ($group=1)
 evaluate ($done=false)
 while ( $done = false ) loop bdihe
 if ( &exist_dihedral_$group = true ) then
     show sum(1) ( &dihedral_$group )
     if ( $result > 0 ) then
       evaluate ($resid=$min_resid_$group)
       evaluate ($nres=$max_resid_$group - $min_resid_$group + 1)
       if ( &dform_$group = "a-form" ) then
         evaluate ($dihedral_alpha=-70)
	 evaluate ($error_alpha=50)
         evaluate ($dihedral_beta=180)
	 evaluate ($error_beta=50)
         evaluate ($dihedral_gamma=60)
	 evaluate ($error_gamma=35)
         evaluate ($dihedral_delta=81)
         evaluate ($error_delta=20)
         evaluate ($dihedral_zeta=-85)
	 evaluate ($error_zeta=50)
         evaluate ($dihedral_epsilon=180)
	 evaluate ($error_epsilon=35)
       elseif ( &dform_$group = "b-form" ) then
         evaluate ($dihedral_alpha=-63.6)
	 evaluate ($error_alpha=6)
         evaluate ($dihedral_beta=176)
	 evaluate ($error_beta=7)
         evaluate ($dihedral_gamma=51.4)
	 evaluate ($error_gamma=7)
         evaluate ($dihedral_delta=128)
         evaluate ($error_delta=13)
         evaluate ($dihedral_epsilon=-171.7)
	 evaluate ($error_epsilon=3.7)
         evaluate ($dihedral_zeta=-103.8)
	 evaluate ($error_zeta=10)
       elseif ( &dform_$group = "other" ) then
         evaluate ($dihedral_alpha=&dihedral_alpha_$group)
	 evaluate ($error_alpha=&error_alpha_$group)
         evaluate ($dihedral_beta=&dihedral_beta_$group)
	 evaluate ($error_beta=&error_beta_$group)
         evaluate ($dihedral_gamma=&dihedral_gamma_$group)
	 evaluate ($error_gamma=&error_gamma_$group)
         evaluate ($dihedral_delta=&dihedral_delta_$group)
	 evaluate ($error_delta=&error_delta_$group)
	 evaluate ($dihedral_zeta=&dihedral_zeta_$group)
	 evaluate ($error_zeta=&error_zeta_$group)
         evaluate ($dihedral_epsilon=&dihedral_eps_$group)
	 evaluate ($error_epsilon=&error_eps_$group)
       end if

       for $id in id ( &dihedral_$group and tag ) loop resid
         show (segid) (id $id)
         evaluate ($segid=$result)
         show (resid) ( id $id )
         evaluate ($resid=decode($result))
         if ($resid > $min_resid_$group) then
           evaluate ($rprec = $resid - 1)
           restraint dihedral
             ! alpha
	     assign  ( segid $segid and resid $rprec and name O3' )
                     ( segid $segid and resid $resid and name P )
                     ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' ) 
                                                       1.0 $dihedral_alpha $error_alpha 2
             ! beta
	     assign  ( segid $segid and resid $resid and name P )
                     ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' ) 
                                                       1.0 $dihedral_beta $error_beta 2
             scale 200.0
           end
         end if

         restraints dihedral
             ! gamma
	     assign  ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' ) 
                                                       1.0 $dihedral_gamma $error_gamma 2
             !delta
	     assign  ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' ) 
                                                       1.0 $dihedral_delta $error_delta 2	         
	     scale=200.0
           end

	 if ($resid < $max_resid_$group) then
           evaluate ($rfoll = $resid + 1)
           restraint dihedral
             ! epsilon
	     assign  ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' )
                     ( segid $segid and resid $rfoll and name P ) 
                                                       1.0 $dihedral_epsilon $error_epsilon 2
             ! zeta
	     assign  ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' )
                     ( segid $segid and resid $rfoll and name P )
                     ( segid $segid and resid $rfoll and name O5' ) 
                                                       1.0 $dihedral_zeta $error_zeta 2
             scale 200.0
           end
         end if
       end loop resid
     end if
   else
     evaluate ($done=true)
   end if
   evaluate ($group=$group+1)
  end loop bdihe
 end if
flags include cdih end

{- C1'-C1' virtual bond length restraints -}

noe
   class hres
   averaging hres cent
   potential hres square
   sqconstant hres 1.
   sqexponent hres 2
   scale hres 70.
 end           

if (&dna_pick_c1 = true) then
  evaluate ($pair=1)
  evaluate ($done=false)
  while ( $done = false ) loop noe
   if ( &exist_base_a_$pair = true ) then
     if ( &exist_base_b_$pair = true ) then
       show ( resname ) ( &base_a_$pair and name C1' ) 
       evaluate ($ares=$result)
       show ( resname ) ( &base_b_$pair and name C1' ) 
       evaluate ($bres=$result)
        pick bond
			(&base_a_$pair and name C1') 
			(&base_b_$pair and name C1')
	   geometry
        evaluate ($c1c1=$result)
        noe
        assign (&base_a_$pair and name C1') 
               (&base_b_$pair and name C1') $c1c1 &c1_low &c1_up 
        end
     end if
   else
     evaluate ($done=true)
   end if         
     evaluate ($pair=$pair+1)
  end loop noe
 else
end if
flags include noe end

{- Watson-Crick base pairing -}

 noe
   class hres
   averaging hres cent
   potential hres square
   sqconstant hres 1.
   sqexponent hres 2
   scale hres 70.
 end           

 if (&dna_pick_wc = true) then
  evaluate ($pair=1)
  evaluate ($done=false)
  while ( $done = false ) loop noe
   if ( &exist_base_a_$pair = true ) then
     if ( &exist_base_b_$pair = true ) then
       show ( resname ) ( &base_a_$pair and name C1' ) 
       evaluate ($ares=$result)
       show ( resname ) ( &base_b_$pair and name C1' ) 
       evaluate ($bres=$result)
       if ( $ares = THY ) then
        pick bond
                  (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n6)
		  geometry
        evaluate ($o4n6=$result)
        pick bond
                  (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1)
		  geometry
        evaluate ($n3n1=$result)
        pick bond
                  (&base_a_$pair and name h3) 
                  (&base_b_$pair and name n1)
		  geometry
        evaluate ($h3n1=$result)
        pick bond
                  (&base_a_$pair and name o2) 
                  (&base_b_$pair and name h2)
		  geometry
        evaluate ($o2h2=$result)
        pick bond
                  (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n1)
		  geometry
        evaluate ($o4n1=$result)
        pick bond
                  (&base_a_$pair and name o2) 
                  (&base_b_$pair and name n1)
		  geometry
        evaluate ($o2n1=$result)
       elseif ( $ares = URI ) then
        pick bond
                  (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n6)
		  geometry
        evaluate ($o4n6=$result)
        pick bond
                  (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1)
		  geometry
        evaluate ($n3n1=$result)
        pick bond
                  (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n1)
		  geometry
        evaluate ($o4n1=$result)
        pick bond
                  (&base_a_$pair and name o2) 
                  (&base_b_$pair and name n6)
		  geometry
        evaluate ($o2n6=$result)
       elseif ( $ares = ADE ) then
        pick bond
                  (&base_b_$pair and name o4) 
                  (&base_a_$pair and name n6)
		  geometry
        evaluate ($o4n6=$result)
        pick bond
                  (&base_b_$pair and name n3) 
                  (&base_a_$pair and name n1)
		  geometry
        evaluate ($n3n1=$result)
        pick bond
                  (&base_b_$pair and name h3) 
                  (&base_a_$pair and name n1)
		  geometry
        evaluate ($h3n1=$result)
        pick bond
                  (&base_b_$pair and name o2) 
                  (&base_a_$pair and name h2)
		  geometry
        evaluate ($o2h2=$result)
        pick bond
                  (&base_b_$pair and name o4) 
                  (&base_a_$pair and name n1)
		  geometry
        evaluate ($o4n1=$result)
        pick bond
                  (&base_b_$pair and name o2) 
                  (&base_a_$pair and name n1)
		  geometry
        evaluate ($o2n1=$result)
       elseif ( $ares = CYT ) then
        pick bond
                  (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1)
		  geometry
        evaluate ($n3n1=$result)
        pick bond
                  (&base_a_$pair and name n3) 
                  (&base_b_$pair and name h1)
		  geometry
        evaluate ($n3h1=$result)
        pick bond
                  (&base_a_$pair and name n4) 
                  (&base_b_$pair and name o6)
		  geometry
        evaluate ($n4o6=$result)
        pick bond
                  (&base_a_$pair and name o2) 
                  (&base_b_$pair and name n2)
		  geometry
        evaluate ($o2h2=$result)
        pick bond
                  (&base_a_$pair and name n3) 
                  (&base_b_$pair and name o6)
		  geometry
        evaluate ($n3o6=$result)
        pick bond
                  (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n2)
		  geometry
        evaluate ($n3n2=$result)
       elseif ( $ares = GUA ) then
        pick bond
                  (&base_b_$pair and name n3) 
                  (&base_a_$pair and name n1)
		  geometry
        evaluate ($n3n1=$result)
        pick bond
                  (&base_b_$pair and name n3) 
                  (&base_a_$pair and name h1)
		  geometry
        evaluate ($n3h1=$result)
        pick bond
                  (&base_b_$pair and name n4) 
                  (&base_a_$pair and name o6)
		  geometry
        evaluate ($n4o6=$result)
        pick bond
                  (&base_b_$pair and name o2) 
                  (&base_a_$pair and name n2)
		  geometry
        evaluate ($o2n2=$result)
        pick bond
                  (&base_b_$pair and name n3) 
                  (&base_a_$pair and name o6)
		  geometry
        evaluate ($n3o6=$result)
        pick bond
                  (&base_b_$pair and name n3) 
                  (&base_a_$pair and name n2)
		  geometry
        evaluate ($n3n2=$result)

       end if
       noe
         if ( $ares = THY ) then
           assign (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n6) $o4n6 &wc_low &wc_up 
           assign (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1) $n3n1 &wc_low &wc_up
           assign (&base_a_$pair and name h3) 
                  (&base_b_$pair and name n1) $h3n1 &wc_low &wc_up 
           assign (&base_a_$pair and name o2) 
                  (&base_b_$pair and name h2) $o2h2 &wc_low &wc_up
           assign (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n1) $o4n1 &wc_low &wc_up
           assign (&base_a_$pair and name o2) 
                  (&base_b_$pair and name n1) $o2n1 &wc_low &wc_up
         elseif ( $ares = URI ) then
           assign (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n6) $o4n6 &wc_low_uri &wc_up_uri 
           assign (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1) $n3n1 &wc_low_uri &wc_up_uri
           assign (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n1) $o4n1 &wc_low_uri &wc_up_uri 
           assign (&base_a_$pair and name o2) 
                  (&base_b_$pair and name n6) $o2n6 &wc_low_uri &wc_up_uri
         elseif ( $ares = ADE ) then
           assign (&base_b_$pair and name o4) 
                  (&base_a_$pair and name n6) $o4n6 &wc_low &wc_up 
           assign (&base_b_$pair and name n3) 
                  (&base_a_$pair and name n1) $n3n1 &wc_low &wc_up
           assign (&base_b_$pair and name h3) 
                  (&base_a_$pair and name n1) $h3n1 &wc_low &wc_up 
           assign (&base_b_$pair and name o2) 
                  (&base_a_$pair and name h2) $o2h2 &wc_low &wc_up
           assign (&base_b_$pair and name o4) 
                  (&base_a_$pair and name n1) $o4n1 &wc_low &wc_up
           assign (&base_b_$pair and name o2) 
                  (&base_a_$pair and name n1) $o2n1 &wc_low &wc_up
         elseif ( $ares = CYT ) then
           assign (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1) $n3n1 &wc_low &wc_up 
           assign (&base_a_$pair and name n3) 
                  (&base_b_$pair and name h1) $n3h1 &wc_low &wc_up
           assign (&base_a_$pair and name n4)
                  (&base_b_$pair and name o6) $n4o6 &wc_low &wc_up
           assign (&base_a_$pair and name o2)
                  (&base_b_$pair and name n2) $o2n2 &wc_low &wc_up 
           assign (&base_a_$pair and name n3)
                  (&base_b_$pair and name o6) $n3o6 &wc_low &wc_up
           assign (&base_a_$pair and name n3)
                  (&base_b_$pair and name n2) $n3n2 &wc_low &wc_up
         elseif ( $ares = GUA ) then
           assign (&base_b_$pair and name n3) 
                  (&base_a_$pair and name n1) $n3n1 &wc_low &wc_up 
           assign (&base_b_$pair and name n3) 
                  (&base_a_$pair and name h1) $n3h1 &wc_low &wc_up
           assign (&base_b_$pair and name n4)
                  (&base_a_$pair and name o6) $n4o6 &wc_low &wc_up
           assign (&base_b_$pair and name o2)
                  (&base_a_$pair and name n2) $o2n2 &wc_low &wc_up 
           assign (&base_b_$pair and name n3)
                  (&base_a_$pair and name o6) $n3o6 &wc_low &wc_up
           assign (&base_b_$pair and name n3)
                  (&base_a_$pair and name n2) $n3n2 &wc_low &wc_up
         end if
       end
     end if
   else
     evaluate ($done=true)
   end if         
   evaluate ($pair=$pair+1)
  end loop noe

 else

 evaluate ($pair=1)
 evaluate ($done=false)
 while ( $done = false ) loop noe
   if ( &exist_base_a_$pair = true ) then
     if ( &exist_base_b_$pair = true ) then
       show ( resname ) ( &base_a_$pair and name C1' ) 
       evaluate ($ares=$result)
       show ( resname ) ( &base_b_$pair and name C1' ) 
       evaluate ($bres=$result)
       noe
         if ( $ares = THY ) then
           assign (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n6) 2.89 0.2 0.2 
           assign (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1) 2.92 0.2 0.2
           assign (&base_a_$pair and name h3) 
                  (&base_b_$pair and name n1) 1.87 0.2 0.2 
           assign (&base_a_$pair and name o2) 
                  (&base_b_$pair and name h2) 2.94 0.2 0.2
           assign (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n1) 3.69 0.2 0.2
           assign (&base_a_$pair and name o2) 
                  (&base_b_$pair and name n1) 3.67 0.2 0.2	    
         elseif ( $ares = URI ) then
           assign (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n6) 2.95 0.01 0.01 
           assign (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1) 2.82 0.01 0.01
           assign (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n1) 3.63 0.01 0.01 
           assign (&base_a_$pair and name o2) 
                  (&base_b_$pair and name n6) 5.40 0.01 0.01
         elseif ( $ares = ADE ) then
           assign (&base_b_$pair and name o4) 
                  (&base_a_$pair and name n6) 2.89 0.2 0.2 
           assign (&base_b_$pair and name n3) 
                  (&base_a_$pair and name n1) 2.92 0.2 0.2
           assign (&base_b_$pair and name h3) 
                  (&base_a_$pair and name n1) 1.87 0.2 0.2 
           assign (&base_b_$pair and name o2) 
                  (&base_a_$pair and name h2) 2.94 0.2 0.2
           assign (&base_b_$pair and name o4) 
                  (&base_a_$pair and name n1) 3.69 0.2 0.2
           assign (&base_b_$pair and name o2) 
                  (&base_a_$pair and name n1) 3.67 0.2 0.2
         elseif ( $ares = CYT ) then
           assign (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1) 2.87 0.2 0.2 
           assign (&base_a_$pair and name n3) 
                  (&base_b_$pair and name h1) 1.86 0.2 0.2
           assign (&base_a_$pair and name n4)
                  (&base_b_$pair and name o6) 2.81 0.2 0.2
           assign (&base_a_$pair and name o2)
                  (&base_b_$pair and name n2) 2.81 0.2 0.2 
           assign (&base_a_$pair and name n3)
                  (&base_b_$pair and name o6) 3.58 0.2 0.2
           assign (&base_a_$pair and name n3)
                  (&base_b_$pair and name n2) 3.63 0.2 0.2
         elseif ( $ares = GUA ) then
           assign (&base_b_$pair and name n3) 
                  (&base_a_$pair and name n1) 2.87 0.2 0.2 
           assign (&base_b_$pair and name n3) 
                  (&base_a_$pair and name h1) 1.86 0.2 0.2
           assign (&base_b_$pair and name n4)
                  (&base_a_$pair and name o6) 2.81 0.2 0.2
           assign (&base_b_$pair and name o2)
                  (&base_a_$pair and name n2) 2.81 0.2 0.2 
           assign (&base_b_$pair and name n3)
                  (&base_a_$pair and name o6) 3.58 0.2 0.2
           assign (&base_b_$pair and name n3)
                  (&base_a_$pair and name n2) 3.63 0.2 0.2
         end if
       end
     end if
   else
     evaluate ($done=true)
   end if         
   evaluate ($pair=$pair+1)
 end loop noe
 end if
 flags include noe end

set message=off echo=off end""")		

if __name__ == '__main__':

	"""Running from the command line"""
	from optparse import *

	"""Parse command line arguments"""
	metadict = {}
	paramdict = {'showonexec':'False','inputfrom':'self','autogenerateGui':'False'}
	option_dict = CommandlineOptionParser().option_dict

	for key in option_dict:
		paramdict[key] = option_dict[key]
	del option_dict

	"""Envoce main functions"""
	PluginCore(paramdict, metadict, inputlist=paramdict['input'])
	sys.exit(0)		
