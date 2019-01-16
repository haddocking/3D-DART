#!/usr/bin/env python2.7

USAGE= """
==========================================================================================

Author:				Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center
					for Biomolecular Research, Utrecht university, The Netherlands
Copyright (C):		2006 (DART project)
DART version:		1.2 (25-11-2008)
DART plugin: 		NAContacts.py
Plugin function:	Analyze contacts between protein and nucleic acids. Optional is
					a multicontact analysis. It uses nucplot for the actual 
					calculation of the contacts. Note that for nucplot to work the 
					chain ID has to be placed in the default location.
Dependencies:		Nucplot

==========================================================================================
"""

"""Import modules"""
import os, sys, re
from math import sqrt
from time import ctime

def PluginXML():
	PluginXML = """ 
<metadata>
 <name>Calculate protein/nucleic-acid contacts</name>
 <input type="Filetype">.pdb</input>
 <output type="Filetype">.bond</output>
</metadata>
<parameters>
 <option type="useplugin" form="hidden" text="None">True</option>
 <option type="inputfrom" form="hidden" text="None">1</option>
 <option type="nucplot" form="checkbox" text="NUCLOT driven contact analysis">False</option>
 <option type="contact" form="checkbox" text="Custom contact analysis">False</option>
 <option type="cutoff" form="text" text="Contacts distance cutoff">5</option>
</parameters>"""
	
	return PluginXML

def PluginGui():
	pass

def PluginCore(paramdict, metadict, inputlist):
	
	for files in inputlist:
		if os.path.basename(files) == 'nbcontacts.disp':
			print "--> Found HADDOCK nbcontats file in input. Calculating contact statistics"
			AnaNBcontacts()
			print "    * Statistics writen to contacts.stat"
		else:
			pass
	
	if paramdict['nucplot'] == 'True' or paramdict['nucplot'] == True:
		
		"""
		Perform nucplot analysis on all files. NOTE: files have to have chain ID in
		propper place otherwise nucplot will fail.
		"""
		
		for files in inputlist: 							
			NucplotAnalyze(files)									
						
	if paramdict['contact'] == 'True' or paramdict['contact'] == True:
		
		"""
		Custom calculation of interactomic distances
		"""
		
		if paramdict['cutoff'] == None:
			paramdict['cutoff'] = 5.0
		else:
			pass
			
		for files in inputlist:
			pdb = Structure(files)
			mol = []
			
			nrchain = 0
			for chain in pdb.nucleotide_chains:
				nrchain = nrchain+1
				mol.append(molecule_assign(nrchain,chain))
			for chain in pdb.peptide_chains:
				nrchain = nrchain+1
				mol.append(molecule_assign(nrchain,chain))
				
			contacts=CustomContact(mol[0],mol[1],mol[2],mol[3],paramdict['cutoff'])

 			for i in contacts.iterkeys():
  				print i[0],i[1],i[2],i[3],i[4],i[5],contacts[i][0],contacts[i][1],contacts[i][2]

								

#================================================================================================================================#
# 					PLUGIN SPECIFIC DEFINITIONS BELOW THIS LINE						 #
#================================================================================================================================#

def NucplotAnalyze(files):
	
	cmd = metadict['dependencies']+" "+(files)				
        os.system(cmd)  							    			   

        nucplot_files = ['hb2.log','hbadd.log','nb2.log','nucplot.ps']  	   			   
        extensions = ['.hb2log','.hbaddlog','.nb2log','.ps']			    			   
        for n in nucplot_files: 						    			   
        	if os.path.isfile(n):						    			   
        		basename,extension = os.path.splitext(files)		    			   
        		ext = extensions[nucplot_files.index(n)]		    			   
        		FileRootRename(n,ext,basename)  			    			   
        	else:								    			   
        		pass							    			   
												
	if os.path.isfile('nucplot.par'):						
		os.remove('nucplot.par')						
	else:										
		pass						

class Molecule:

	def __init__(self,id,chain):
   		self.id = id
   		self.chain=chain

def molecule_assign(nr,chain):
	return Molecule(nr,chain)

def FileRootRename(infile,extension,basename):
	outfile = basename+extension
	os.rename(infile,outfile)

def calc_dist(vect1,vect2):
 	d=0.0
 	for i in range(3):
  		d=d+(vect1[i]-vect2[i])**2

 	d=sqrt(d)
 	return d  

def CustomContact(molecule1,molecule2,molecule3,molecule4,cutoff):
    chain1=molecule1.chain
    chain2=molecule2.chain
    chain3=molecule3.chain
    chain4=molecule4.chain
    dict={}
    for i in chain1.residues:
        for atom1 in i:
            for p in chain3.residues:
                for atom2 in p:
		  if (not (atom1.name[0]=="H") and (not atom2.name[0]=="H")):
		     d=calc_dist(atom1.position,atom2.position)
       	             if (d<cutoff):
		      if not dict.has_key((molecule1.id,i.name,i.number,molecule2.id,p.name,p.number)):
		       dict.update({(molecule1.id,i.name,i.number,molecule2.id,p.name,p.number):(d,atom1.name,atom2.name)})
                       #print molecule1.id,i.name,i.number,atom1.name,molecule2.id,p.name,p.number,atom2.name,d
		      else:
		       tmp=dict[(molecule1.id,i.name,i.number,molecule2.id,p.name,p.number)]
		       if (d<tmp[0]):
		        newtmp=(d,atom1.name,atom2.name)
		        dict.update({(molecule1.id,i.name,i.number,molecule2.id,p.name,p.number):newtmp})
    return dict	

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
		parser.add_option( "-n", "--nucplot", action="store_true", dest="nucplot", default=False, help="Perform nucplot analysis on file")
		parser.add_option( "-c", "--contact", action="store_true", dest="contact", default=False, help="Determine contacts in structure within given cutoff range (-b option)")
		parser.add_option( "-b", "--cutoff", action="store", dest="cutoff", type="float", help="Define cutoff range for contact analysis, default = 5A")
		
		(options, args) = parser.parse_args()
		
		self.option_dict['input'] = options.inputfile
		self.option_dict['nucplot'] = options.nucplot
		self.option_dict['contact'] = options.contact
		self.option_dict['cutoff'] = options.cutoff
			
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

class AnaNBcontacts:

	"""Calculate statistics for intermolecular contatcs from the HADDOCK nbcontacts file"""
	
	def __init__(self):
		
		self.chainone = []
		self.residnrone = []
		self.residtypeone =[]
		self.atomtypeone = []
		self.chaintwo = []
		self.residnrtwo = []
		self.residtypetwo = []
		self.atomtypetwo = []
		self.distance = []
		self.files = []
		
		self.pairdict = {}
		
		self.ReadNBcontacts()
		self.MatchPaires()
		self.ReportStats()
	
	def _msort(self, liste, indice):
		
		""" This function sorts a list regarding values of the indice 'indice'. Indice start from 0"""
		
		tmp=[[tbl[indice]]+[tbl] for tbl in liste]
		tmp.sort(reverse=True)
		liste=[cl[1] for cl in tmp]
		del tmp
		
		return liste
	
	def ReadNBcontacts (self, debug=0):
		
		readfile = file('nbcontacts.disp','r')
		lines = readfile.readlines()
		
		self.ReadNBLines(lines, debug)
	
	def ReadNBLines(self, lines, debug=0):
		
		contacts = re.compile('( atoms)')
		files = re.compile('(PREVIT)')
		for line in lines:	
			if contacts.match(line):
				params = []
				for n in range(1,13):
					params.append(line.split()[n])
				self.chainone.append(params[0])
				self.residnrone.append(int(params[1]))
				self.residtypeone.append(params[2])
				self.atomtypeone.append(params[3])
				self.chaintwo.append(params[6])
				self.residnrtwo.append(int(params[7]))
				self.residtypetwo.append(params[8])
				self.atomtypetwo.append(params[9])
				self.distance.append(float(params[11]))
			if files.match(line):
				self.files.append(line)
	
	def MatchPaires(self):	
		
		pairlist = []
		
		count = 0
		while count < len(self.residnrone):
			pairlist.append([self.residnrone[count],self.chainone[count],self.residtypeone[count],self.residnrtwo[count],self.chaintwo[count],self.residtypetwo[count]])
			count = count+1	
		
		count = 0
		while count < len(pairlist):
			instance = pairlist[count]
			if self.pairdict.has_key((instance[0],instance[1],instance[2],instance[3],instance[4],instance[5])):
				pass
			else:
				index = []
				i = 0
				while i < len(pairlist):
					if pairlist[i] == instance:
						index.append(i)
					i = i+1
				self.pairdict[instance[0],instance[1],instance[2],instance[3],instance[4],instance[5]] = index	
			count = count+1	
	
	def ReportStats(self):
		
		"""Write contact statistics to file contact.stat"""
		
		lines = []
		for key, value in self.pairdict.items():
			lines.append([key[1],key[0],key[2],key[4],key[3],key[5],len(value),((float(len(value)))/(float(len(self.files)))*100)])
		sorted_lines = self._msort(lines, 6)
		
		outfile = file('contacts.stat','w')
		
		outfile.write('*************************************************************************************************************************\n')
		outfile.write('Contact analysis from HADDOCK nbcontacts.disp file for %i structures\n' % (len(self.files)))
		outfile.write('Time and date: %s\n' % ctime())
		outfile.write('*************************************************************************************************************************\n')
		outfile.write('\nStatistical for residue/residue contacts\n')
		outfile.write('chain 1  resid nr.  resid type  chain 2  resid nr.  resid type  count  frac(%)\n')
		
		for paires in sorted_lines:
			outfile.write('%1s%11i%11s%10s%10i%12s%10i%8.2f\n' %
			(paires[0],paires[1],paires[2],paires[3],paires[4],paires[5],paires[6],paires[7]))
				
		outfile.close()	
				 
		
if __name__ == '__main__':

	"""For testing purposes"""
	from optparse import *
	
	"""Parse command line arguments"""
	metadict = {'dependencies':'/Applications/Science/DART/third-party/nucplot/nucplot'}
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
	PluginCore(paramdict, metadict, inputlist)
	
	"""Say goodbye"""
	print "--> Thanks for using NAContacts, bye"
	sys.exit(0)
