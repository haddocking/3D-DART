#!/usr/bin/env python2.7

USAGE = """
==========================================================================================

Author:			Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center
			for Biomolecular Research, Utrecht university, The Netherlands
Copyright (C):		2006 (DART project)
DART version:		1.2  (25-11-2008)
DART module: 		RunDART.py
Input:			A premade workflow, a sequence of plugins to generate a workflow 
                    	from and a list of files as input for the workflow.
Output:			A directory structure containing the output of every plugin in the
 			workflow sequence. 
Module function:	Main script for running DART.
Examples:		dart -dp X3DNAanalyze NABendAnalyze
			dart -w NAanalyze
Module depenencies:	Standard python2.4 modules, Numeric
			
==========================================================================================
"""

DART_VERSION = "1.0"

"""Import modules"""
import sys, os
import time
from time import ctime

"""Setting variables"""
base, dirs = os.path.split(os.path.join(os.getcwd(), __file__))
split = base.split('/')
DARTbase = '/'.join(split[0:-1])

if os.uname()[0] == 'Darwin':
	os.environ["X3DNA"] = "%s/software/X3DNA-mac/" % DARTbase
	os.environ["PATH"] = os.getenv("PATH")+":%s/software/X3DNA-mac/bin" % DARTbase
else: 
	os.environ["X3DNA"] = "%s/software/X3DNA-linux/" % DARTbase
	os.environ["PATH"] = os.getenv("PATH")+":%s/software/X3DNA-linux/bin" % DARTbase
    
if base in sys.path: pass
else: sys.path.append(base)

sys.path.append(base+"/system/")

from system.CommandLineParser import CommandlineOptionParser 
from system import FrameWork

def SystemChecks():
	print "--> Performing System Checks:"
	
	"""Python version check"""
	version = float(sys.version[:3])
	if version < 2.3:
		print "   * Python version 2.3 or higher required"
		print "   * Current version is:", sys.version[:5]
		ExitMessage()
	else:
		print "   * Python version is:", sys.version[:5]

	"""Check Numeric/numarray package""" 
	try:
		import Numeric
		print "   * Importing Numeric package succesful"
	except:
		print "   * Could not import Numeric package try numpy"
		try:
			import numpy as Numeric
			print "   * Importing numpy package succesfull"
		except:
			print "   * Could not import Numeric or numpy package"
			ExitMessage()
	
	"""General messages"""
	print "--> Your current working directory is:", os.getcwd()

def ExitMessage():
	print "-"*110
	print "3D-DART workflow sequence is executed succesfully"
	print "-"*110
	sys.exit(0)

if __name__ == "__main__":
	SystemChecks()
	options = CommandlineOptionParser(DARTdir=base)
	opt_dict = options.option_dict
	FrameWork.PluginExecutor(opt_dict=opt_dict, DARTdir=base)
	ExitMessage()
