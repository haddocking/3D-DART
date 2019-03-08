#!/usr/bin/env python2.7

"""
==========================================================================================

Author:			    Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center for 
			        Biomolecular Research, Utrecht university, The Netherlands.
Copyright (C):		2006 (DART project)
DART version:		1.2 (25-11-2008)

Diverse but usefull functions
==========================================================================================
"""

"""Import modules"""
import os

def MakeBackup(infile,report=False):
	
	"""Make backup of files or directories by checking if file (or backup as _*) is there and rename"""
	
	if os.path.isfile(infile) or os.path.isdir(infile):
		i = 1
		while i < 100:
			if os.path.isfile(infile+'_'+str(i)) or os.path.isdir(infile+'_'+str(i)):
				i = i+1
			else:
				os.rename(infile, infile+'_'+str(i))
				break

	if report == True:
		return infile

def FileRootRename(infile,extension,basename):

	"""Rename a file but preserve extension"""

	outfile = basename+extension
	os.rename(infile,outfile)	

def TransformDash(n):
	
	"""Transform dashes often found in 3DNA tables to floats"""
	
	if n == '---' or n == '----':
		return float(0)
	else:
		return float(n)

def GetFullPath(inputfile):

	"""Return the full path of the file(s)"""
	
	try:
		filelist = []
		for files in inputfile:
			filelist.append(os.path.abspath(files))
		return filelist
	except:
		return os.path.abspath(files)

def RenameFilepath(inputfile,path=None,basename=None,extension=None):

	"""Rename a filepath in several ways: path, basename and/or extension"""
	
	orpath = os.path.dirname(inputfile)
	orbasename = os.path.splitext(os.path.basename(inputfile))[0] 
	orextension = os.path.splitext(os.path.basename(inputfile))[1] 
	
	newfile = ""
	
	if path == None:	# Set path
		newfile = orpath+"/"
	else:
		newfile = path+"/"
	if basename == None:	# Set file basename
		newfile = newfile+orbasename
	else:
		newfile = newfile+basename
	if extension == None:	# Set file extension
		newfile = newfile+orextension
	else:
		newfile = newfile+extension
		
	return newfile				

