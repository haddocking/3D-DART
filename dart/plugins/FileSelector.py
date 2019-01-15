#!/usr/bin/env python2.7

USAGE = """
==========================================================================================

Author:				Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center
					for Biomolecular Research, Utrecht university, The Netherlands
Copyright (C):		2006 (DART project)
DART version:		1.2 (25-11-2008)
DART plugin: 		FileSelector.py
Input:				None
Output:				List of selected files
Plugin excecution:	As part of a DART batch sequence
Plugin function:	A basic popup in which you can select files as inut for follow
					up plugins
Dependencies:		Standard python2.3 or higher modules

==========================================================================================
"""

"""Import modules"""
import os, shutil

def PluginXML():
	
	PluginXML = """ 
<metadata>
 <name>Upload your files</name>
 <input type="Filetype"></input>
 <output type="Filetype">self</output>
</metadata>
<parameters>
 <option type="useplugin" form="hidden" text="None">True</option>
 <option type="inputfrom" form="hidden" text="None">1</option> 
 <option type="upload" form="file" text="Upload your file">None</option>
</parameters>"""
	
	return PluginXML
	
def PluginCore(paramdict, inputlist):
	
	"""FileSelector does not need any input so inputlist is not used"""
	
	if len(inputlist) > 0:
		print "--> Selecting a file(s) for the batch sequence"

		jobdir = os.getcwd()
		for files in inputlist:
			if os.path.isfile(files):
				print("    * Found: %s" % files)	
				fname = os.path.basename(files)
				destination = os.path.join(jobdir,fname)
				shutil.copyfile(files,destination)
			else:
				print("    * Not found: %s" % files)		

if __name__ == '__main__':

	"""For testing purposes"""
	print USAGE
	paramdict = {'inputfrom':1}
	inputlist = ['dummy.test'] #The dummy list as the plugin needs no input

	PluginCore(paramdict, inputlist)
	
