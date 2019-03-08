#!/usr/bin/env python2.7

from __future__ import print_function


USAGE = """
==========================================================================================

Author:			    Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center
                    for Biomolecular Research, Utrecht university, The Netherlands
Copyright (C):		2006 (DART project)
DART version:		1.3  (08-03-2019)
DART module: 		RunDART.py
Input:			    A premade workflow, a sequence of plugins to generate a workflow 
                    	from and a list of files as input for the workflow.
Output:			    A directory structure containing the output of every plugin in the
                        workflow sequence. 
Module function:	Main script for running DART.
Examples:		    dart -dp X3DNAanalyze NABendAnalyze
			        dart -w NAanalyze
Module depenencies:	Standard python2.x modules, Numeric/NumPy
			
==========================================================================================
"""

DART_VERSION = "1.3"

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
    
if base not in sys.path:
    sys.path.append(base)

sys.path.append(os.path.join(base, "/system/"))

from system.CommandLineParser import CommandlineOptionParser 
from system import FrameWork


def system_checks():
    """Checks the system setup"""
    print("--> Performing System Checks:")
	
    """Python version check"""
    version = float(sys.version[:3])
    if version < 2.3:
        print("   * Python version 2.3 or higher required")
        print("   * Current version is: %s" % sys.version[:5])
        exit_message()
    else:
        print("   * Python version is: %s" % sys.version[:5])

    """Check Numeric/numarray package""" 
    try:
        import Numeric
        print("   * Importing Numeric package succesful")
    except:
        print("   * Could not import Numeric package trying NumPy")
        try:
            import numpy as Numeric
            print("   * Importing NumPy package succesfull")
        except:
            print("   * Could not import Numeric or NumPy package")
            exit_message()
	
    """General messages"""
    print("--> Your current working directory is: %s" % os.getcwd())


def exit_message():
    """Exits gently"""
    print("-"*110)
    print("3D-DART workflow sequence is executed succesfully")
    print("-"*110)
    raise SystemExit


if __name__ == "__main__":
    system_checks()
    options = CommandlineOptionParser(DARTdir=base)
    opt_dict = options.option_dict
    FrameWork.PluginExecutor(opt_dict=opt_dict, DARTdir=base)
    exit_message()

