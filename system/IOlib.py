#!/usr/bin/env python2.7

import os,sys,re
from Utils import *
from numpy import *
from Constants import *

def WritePar(database,filename,verbose=False):
	
	"""Write 3DNA base-pair and base-pair step parameter file (*.par)"""
	if filename == None:
		filename = 'parfile'

	if verbose == True:
		outfile = sys.stdout
	else:
		MakeBackup(os.path.splitext(filename)[0]+'.par')
		outfile = file(os.path.splitext(filename)[0]+'.par','w')
		print("    * Writing new parameter file with name %s" % filename)
		
	for param in BASEPAIR_STEPS:	#All first values for base-pair step parameters are always 0
		database.Update(param,float(0),0)
		
	outfile.write(" %s base-pairs\n" % len(database['sequence']))
   	outfile.write("  0  ***local base-pair & step parameters***\n")
        outfile.write("      Shear  Stretch  Stagger Buckle Prop-Tw Opening   Shift  Slide    Rise    Tilt    Roll   Twist\n")

	for n in range(len(database['sequence'])):
		outfile.write("%s %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n" % (database['sequence'][n],database['shear'][n],
		database['stretch'][n],database['stagger'][n],database['buckle'][n],database['proptw'][n],database['opening'][n],
		database['shift'][n],database['slide'][n],database['rise'][n],database['tilt'][n],database['roll'][n],database['twist'][n]))

	if verbose == False:
		outfile.close()

class InputOutputControl:

	"""Acts as a plugin input and output control. Input and output requirements are checked against extensions or files.
	   Generation and use of files within the plugin are controled by this class in a similar manner"""
	
	def __init__(self):
	
		self.checkedinput = {}
	
	def _CheckFile(self,files,requirements):
		
		if requirements == None or requirements == 'None':
			extensions = []
			for n in files:
				ext = os.path.splitext(n)[1]
				if ext in extensions:
					pass
				else:
					extensions.append(ext)
			for extension in extensions:
				self.checkedinput[extension] = []
		else:			
			splitter = re.compile(',')
			requirements = splitter.split(requirements)
			for requirement in requirements:
				self.checkedinput[requirement] = []
		
		for n in files:
			if os.path.isfile(n):
				extension = os.path.splitext(n)[1]
				if self.checkedinput.has_key(extension):
					self.checkedinput[extension].append(n)
				elif self.checkedinput.has_key(os.path.basename(n)):
					self.checkedinput[os.path.basename(n)].append(n)	
			else:
				print "    * InputCheck ERROR: file", n, "not found"   
	
	def CheckInput(self,files,requirements=None):

		filelist = []
		
		if type(files) == type(""):
			filelist.append(files)
		elif type(files) == type([]):
			filelist = files

		for files in filelist:
			if os.path.basename(files) == 'selection.list':
				readfile = file(files,'r')
	       	 		lines = readfile.readlines()
				for line in lines:
					filelist.append(line.strip())
		
		self._CheckFile(filelist,requirements)
	
	def InputUpdate(self,reference,required):
		
		"""Update the dictionary of files"""
		self.checkedinput[required] = []
		
		for n in self.checkedinput[reference]:
			expected = RenameFilepath(n,path=os.getcwd(),extension=required)
			if os.path.isfile(expected):
				self.checkedinput[required].append(expected)
			else:
				print "    * InputCheck ERROR: file", expected, "expected but not found"	
	
	def CheckOutput(self,files,requirements=None):
	
		"""Check if required output is generated"""
		
		if not requirements or requirements == 'self':
			requirements = []
		else:	
			splitter = re.compile(',')
			requirements = splitter.split(requirements)
		
		inputfiles = self.DictToList()
		output_expect = []
		for a in inputfiles:
			basename, extension = os.path.splitext(os.path.basename(a))
			for requirement in requirements:
				if requirement[0] in ['.','_']: output_expect.append(basename+requirement)
		for requirement in requirements:
			if not requirement[0] in ['.','_']: output_expect.append(requirement)	
		
		for a in output_expect:
			if not os.path.isfile(a): print "    * WARNING: file:", a, "is not present in the output"	
		
		"""True output"""
		output_true = []
		for a in files:
			output_true.append(os.path.join(os.getcwd(),a))

		return output_true
			
	def DictToList(self):
	
		"""Return dictionary as plain list of files"""
		
		filelist = []
		
		for n in self.checkedinput:
			filelist = filelist+self.checkedinput[n]
		
		return filelist				
				
class DatabaseDeamon:
	
	"""This module allows for the quick construction, maipualtion and export of databses.
      	   The database is constructed as aq library, data is stored depending on type as either
  	   an array of values (floats) or a list of strings"""
	    
	def __init__(self):
	
		self.database = {}
	
	def __getitem__(self,key):
		
		"""Retrieve data from database"""
		
		return self.database[key]	
		
	def _TypeCheck(self,data):		

		"""Check input on type 'float','integer','string','list',or single value. If None
		   of these return None""" 
		
		checked = []
		datatype = None
		
		if type(data) == type([]):
			for items in data:
				try:
					checked.append(float(items))
				except:
					checked.append(items)	
		elif isinstance(data, float) or isinstance(data, int):
			checked.append(float(data))	
		elif isinstance(data, str):
			checked.append(data)
		elif isinstance(data, type(array(0))):
			checked = data
		else:
			checked = None
		
		return checked	
	
	def Update(self,key,item,index=None):
		
		"""Update data in database. Update complete dataset for a given database entry or update a value
		   in the dataset."""
		
		if not index == None:
			if self.database.has_key(key) == True:
				self.database[key][index] = self._TypeCheck(item)[0]
				try:
					self.database[key][index] = self._TypeCheck(item)[0]
				except:
					print("    * ERROR: failed database update of dataset %s at index %i\n" % (key,index))
		else:	
			if self.database.has_key(key) == True:
				if len(item) == len(self.database[key]):
					del self.database[key]
					self.database[key] = self._TypeCheck(item)
				else:
					print("    * DATABASE-ERROR: new list of items does not match length of old list")	
		
	def Load(self,name,data):
		
		"""Load checked data in database"""
	
		self.database[name] = self._TypeCheck(data)

	
	
				 
