#!/usr/bin/env python3

"""A basic popup in which you can select files as input for follow up plugins"""

import os
import shutil

# Logging
import logging
logging.basicConfig(format='%(name)s [%(levelname)s] %(message)s', level=logging.INFO)
log = logging.getLogger("FileSelector")


def PluginXML():
	xml = """ 
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
	return xml


def PluginCore(paramdict, inputlist):
	"""FileSelector does not need any input so inputlist is not used"""
	if len(inputlist) > 0:
		log.info("Selecting a file(s) for the batch sequence")
		jobdir = os.getcwd()
		for files in inputlist:
			if os.path.isfile(files):
				log.info("    * Found: %s" % files)	
				fname = os.path.basename(files)
				destination = os.path.join(jobdir, fname)
				shutil.copyfile(files, destination)
			else:
				log.error("    * Not found: %s" % files)		


if __name__ == '__main__':
	"""For testing purposes"""
	paramdict = {'inputfrom':1}
	# The dummy list as the plugin needs no input
	inputlist = ['dummy.test']
	PluginCore(paramdict, inputlist)
