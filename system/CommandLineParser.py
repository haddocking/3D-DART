#!/usr/bin/env python2.7

USAGE = """
==========================================================================================

Author:				Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center
					for Biomolecular Research, Utrecht university, The Netherlands
Copyright (C):		2006 (DART project)
DART version:		1.2  (25-11-2008)
DART module:		CommandLineParser.py
Input:				The DART framework excecutes predefined workflows composed of DART 
					plugins with certain settings (option -w/--workflow). The workflows 
					can be generated using the -d/--dry option. Edit the plugin 
					parameters in the workflow XML file afterwards. By placing the 
					workflow file in the workflow directory of the DART package you can 
					always excecute this workflow from the commandline (it functions 
					as a sort of script repository). 
Output:				Output of every plugin defined in the workflow is placed in a sepperate
 					subdirectory of a parent workflow directory
Examples:			dart -dp X3DNAanalyze NABendAnalyze
					dart -w NAanalyze
Dependencies:		Standard python2.5 modules, Numeric

==========================================================================================
"""

"""import modules"""
import os, glob, sys
from Xpath import Xpath
from optparse import *
from time import ctime
from Utils import GetFullPath,RenameFilepath
from DARTserver import WebServer
from Constants import *

class CommandlineOptionParser:
	
	"""Parses command line arguments using optparse. Parser -h/--help for help, -p/--plugin
	   for custom assembly of workflow on the command line, -l/--list of listing all 
	   available plugins, -w/--workflow for the execution of a pre-supplied DART workflow and
	   -f/--files for suppling input files for the first plugin. Alle commandline options are
	   checked to ensure validtie.""" 
	
	def __init__(self, DARTdir):
		
		self.DARTdir = DARTdir
		self.option_dict = {} 
		
		self.CommandlineOptionParser()
		self.WorkflowInput()
		self.WorkflowXML()
		self.MakeWebForm()
	
	def _TryPluginImport(self, plugin):
		
		"""Try to import the key components of the plugin. This is only
		   for testing the validitie of the plugin. Modules do not stay
		   registered."""
		
		print "    * Try importing plugin:", plugin 
	
		try:
			exec "from plugins import "+plugin
			print "      - Plugin imported"
		except:
			print "      - ERROR: Could not import", plugin
			sys.exit(0)

		try:
			exec "from plugins."+plugin+" import PluginXML"
			print "      - Plugin XML file present"
		except:
			print "      - ERROR: Plugin has no XML data file, this is not valid"
			sys.exit(0)
		try:
			exec "from plugins."+plugin+" import PluginCore as PluginCore"
			print "      - Plugin Core present"
		except:
			print "      - ERROR: Plugin has no Core module"	
			sys.exit(0)
		
	def CommandlineOptionParser(self):
	
		"""Parsing command line arguments"""
	
		usage = "usage: %prog" + USAGE
		parser = OptionParser(usage)

		parser.add_option( "-l", "--list", action="store_true", dest="plugin_list", default=False, help="List available plugins and workflows")
		parser.add_option( "-d", "--dry", action="store_true", dest="dry", default=False, help="Only generate workflow file, do nothing more")		
		parser.add_option( "-w", "--workflow", action="store", dest="workflow", type="string", help="Execute predefined workflow")
		parser.add_option( "-p", "--plugin", action="callback", callback=self.varargs, dest="pluginseq", type="string", help="Execute custom workflow assembled on the command line. You can execute a single plugin by typing '-p pluginname' or a sequence of plugins by typing '-p plugin1 plugin2...'")
		parser.add_option( "-f", "--file", action="callback", callback=self.varargs, dest="filename", type="string", help="Supply one or a list of files as input for the plugin(sequence)")
		parser.add_option( "-s", "--server", action="store_true", dest="server", default=False, help="Generate html webform from workflow xml file to be used in DART web implementation")

		(options, args) = parser.parse_args()
		
		option_dict = {}

		self.option_dict['workflow'] = options.workflow
		self.option_dict['input'] = options.filename
		self.option_dict['pluginseq'] = options.pluginseq
		self.option_dict['dry'] = options.dry
		self.option_dict['server'] = options.server

		if not self.option_dict['input'] == None:
			parser.remove_option('-f')
			self.option_dict['input'].append(self.GetFirstArgument(parser, shorta='-f', longa='--file'))
			self.option_dict['input'] = GetFullPath(self.option_dict['input'])
			
		if not parser.has_option('-f'):
			parser.add_option( "-f", "--file", action="store", dest="dummy2", type="string") #only needs to be here to complete the argument list, not used!
		
		if not self.option_dict['pluginseq'] == None:
			parser.remove_option('-p')
			self.option_dict['pluginseq'].append(self.GetFirstArgument(parser, shorta='-p', longa='--plugin'))
		if options.plugin_list == True:
			self.GetPluginList()
		
		"""Report to the user"""
		self.Welcome()
		print "--> Parsing command-line arguments:"
	
		if not self.option_dict['workflow'] == None:
			print "    * The following 3D-DART batch configuration file will be executed:", self.option_dict['workflow']
		
		if not self.option_dict['pluginseq'] == None:
			print "    * The following command line plugin sequence will be excecuted:" 
			for plugin in self.option_dict['pluginseq']:
				print "      - ", plugin
			
		if not self.option_dict['input'] == None:
			print "    * The following files will be used as input for the batch sequence:"
			for files in self.option_dict['input']:
				print "      - ", files

	def GetFirstArgument(self, parser, shorta, longa):

		"""HACK, optparse has difficulties in variable argument lists. The varargs definition solves this but never reports the first 
		   argument of the list. This definition hacks this issue"""

		parser.add_option( shorta, longa, action="store", dest="temp", type="string", help="Execute custom workflow assembled on the command line. You can execute a single plugin by typing '-p pluginname' or a sequence of plugins by typing '-p plugin1,plugin2...'")

		(options, args) = parser.parse_args()
		parser.remove_option(shorta)
		
		return options.temp

	def Welcome(self):

		from time import ctime
		print "-" * 110
		print "Welcome to 3D-DART version %s  " % DART_VERSION, ctime()
		print "-" * 110
		
	def GetPluginList(self):

		"""List all available plugins in plugin directory and all workflows in workflow directory"""

		self.Welcome()
		print "--> List of available plugins. Excecute plugin with option -h/--help for more information"
		print "    about the function of the plugin"
		plugindir = self.DARTdir+'/plugins' 			     		     
		os.chdir(plugindir)				     		     
		plugins = glob.glob('*.py')			     		     
		plugins.remove('__init__.py')			     		     
		for plugin in plugins:				     		     
			basename,extension = os.path.splitext(plugin)		     
        		print "    * ",basename 		     		     
		
		print "--> List of available workflows:"
		workflowdir = self.DARTdir+'/workflows'
		os.chdir(workflowdir)
		workflows = glob.glob('*.xml')
		for workflow in glob.glob('*.xml'):
			basename,extension = os.path.splitext(workflow)
			print "    * ",basename
		
		sys.exit(0)					     		     

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

	def WorkflowInput(self):
		
		"""Check command line input file string and make list of it"""
	
		if not self.option_dict['input'] == None:
			print "    * Check if all files supplied as input on the command line are present"
			filelist = self.option_dict['input']
			for files in filelist:
				if os.path.isfile(files):
					pass
				else:
					print "      - ERROR: the file:", files, "cannot be found it will be removed from the list"
					filelist.remove(files)
			
			if len(filelist) == 0:
				print "      - ERROR: the -f or --files command line argument was used but no valid files are available"
			
			else:
				print "      - OK"
				self.option_dict['input'] = filelist
				
	def WorkflowXML(self):
		
		"""Generate workflow.xml file at startup. This file allways needs to be present. The source
		   can be a pregenerated workflow file or a sequence of plugins supplied on the command line 
		   with the -p or --plugin option."""
		
		if not self.option_dict['pluginseq'] == None:
			
			workflow_dict = self.WorkflowSequence()
			
			print "    * Generate workflow from command line constructed plugin batch sequence"
			print "    * Check if all plugins are present and of a valid type"
			
			plugins = []
			for n in workflow_dict.values(): 
				if not n in plugins: 
					plugins.append(n)
					self._TryPluginImport(n)
		
			print "    * Writing workflow XML file as workflow.xml"
			
			outfile = file('workflow.xml','w')
			outfile.write("""<?xml version="1.0" encoding="iso-8859-1"?>\n""")
			outfile.write("""<main id="DARTworkflow">\n""")
			outfile.write("<meta>\n")
			outfile.write("<name>workflow.xml</name>\n")
			outfile.write("<datetime>"+ctime()+"</datetime>\n")
			outfile.write("</meta>\n")
			for key in workflow_dict:
				outfile.write("<plugin id='"+workflow_dict[key]+"' job='"+str(key)+"'>")
				exec "from plugins import "+workflow_dict[key]
				exec "outfile.write(\n"+workflow_dict[key]+".PluginXML())"
				outfile.write("\n</plugin>\n")
			outfile.write("</main>\n")
			outfile.close()
		
			self.option_dict['workflow'] = 'workflow.xml'
			
			if self.option_dict['dry'] == True:
				print "    * Option 'dry' is True, only write workflow.xml file"
				sys.exit(0)
			
		elif not self.option_dict['workflow'] == None:
			
			print "    * Check if all plugins in pre-supplied 3D-DART workflow", self.option_dict['workflow'], "are present and of a valid type"
		
			if os.path.isfile(self.option_dict['workflow']):
				pass
			elif os.path.isfile(RenameFilepath(self.option_dict['workflow'],path=self.DARTdir+"/workflows/",extension=".xml")):
				self.option_dict['workflow'] = RenameFilepath(self.option_dict['workflow'],path=self.DARTdir+"/workflows",extension=".xml")
			else:
				print "     * ERROR: the workflow:", self.option_dict['workflow'], "cannot be found"
				sys.exit(0)
				 	
			query = Xpath(self.option_dict['workflow'])
			query.Evaluate(query={1:{'element':'plugin','attr':None}})
			for job in query.nodeselection[1]:
				query.getAttr(node=job,selection=['id'],export='string')					
			workflow = query.result
			query.ClearResult()

			plugins = []
			for n in workflow:
				if not n in plugins:
					plugins.append(n)
					self._TryPluginImport(n)

			print "    * Plugin sequence is valid"
			
	def WorkflowSequence(self):
		
		"""Construct squence of execution from command line string of plugins"""
		
		workflow_dict = {}
		counter = 1
		
		pluginlist =  self.option_dict['pluginseq']
		if len(pluginlist) > 1:
			pluginlist = [pluginlist[len(pluginlist)-1],]+pluginlist[:len(pluginlist)-1]
			if pluginlist[0] == "FileSelector":
				pass
			else:
				print "    * For program integrety reasons the FileSelector plugin needs to be present, including it"
				pluginlist[0:0] = ['FileSelector']
			for plugin in pluginlist:
				workflow_dict[str(counter)] = plugin
				counter += 1
		else:
			if pluginlist[0] == "FileSelector":
				pass
			else:
				print "    * For program integrety reasons the FileSelector plugin needs to be present, including it"
				pluginlist[0:0] = ['FileSelector']
			for plugin in pluginlist:
				workflow_dict[str(counter)] = plugin
				counter += 1
		
		return workflow_dict
	
	def MakeWebForm(self):
	
		"""Call DARTserver.py to make a webform from the workflow.xml file"""
		
		if self.option_dict['server'] == True and not self.option_dict['workflow'] == None:
			print "--> Generating DARTserver compatible HTML webform from the workflow xml file"
			server = WebServer()
			server.MakeWebForm(verbose=False,xml=self.option_dict['workflow'])
		
			print "    * Webform generated. Stopping DART"
			sys.exit(0)			

if __name__ == "__main__":
	
	"""For testing the script"""
	
	DARTdir = '/Applications/Science/DART/'
	options = CommandlineOptionParser(DARTdir)
	opt_dict = options.option_dict
	print opt_dict
