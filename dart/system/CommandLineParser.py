#!/usr/bin/env python3

import os
import glob
import sys
import argparse
from time import ctime
from dart.system.Xpath import Xpath
from dart.system.Utils import RenameFilepath
from dart.system.DARTserver import WebServer
from dart.system.version import __version__


usage = """
==========================================================================================

Authors:            Marc van Dijk, Brian Jimenez-Garcia, Mikael Trellet
                    Department of NMR spectroscopy, Bijvoet Center for 
                    Biomolecular Research, Utrecht university, The Netherlands
Copyright (C):      2006 (DART project)
DART version:       {0}

DART module:        3d_dart.py
Input:              A premade workflow, a sequence of plugins to generate a workflow 
                        from and a list of files as input for the workflow.
Output:             A directory structure containing the output of every plugin in the
                    workflow sequence. 
Module function:    Main script for running DART.

Examples:           3d_dart.py -dp X3DNAanalyze NABendAnalyze
                    3d_dart.py -w NAanalyze

Module depenencies: Python 3.x, NumPy
            
==========================================================================================
""".format(__version__)


# Logging
import logging
logging.basicConfig(format='%(name)s [%(levelname)s] %(message)s', level=logging.INFO)
log = logging.getLogger("3d_dart")


class CommandLineOptionParser:
	
	"""Parses command line arguments using optparse.

	-h/--help for help
	-p/--plugin for custom assembly of workflow on the command line
	-l/--list of listing all available plugins
	-w/--workflow for the execution of a pre-supplied DART workflow
	-f/--files for suppling input files for the first plugin

	""" 
	
	def __init__(self, dart_path='.'):
		
		self.darth_path = dart_path
		self.option_dict = {} 
		
		self.parse_command_line()
		
		#self.WorkflowXML()
		#self.MakeWebForm()
	
	def _test_plugin_import(self, plugin):
		"""Try to import the key components of the plugin. This is only
		   for testing the validity of the plugin. Modules do not stay
		   registered."""
		
		print("    * Try importing plugin:", plugin)
	
		try:
			exec("from plugins import "+plugin)
			print("      - Plugin imported")
		except:
			print("      - ERROR: Could not import", plugin)
			sys.exit(0)

		try:
			exec("from plugins."+plugin+" import PluginXML")
			print("      - Plugin XML file present")
		except:
			print("      - ERROR: Plugin has no XML data file, this is not valid")
			sys.exit(0)
		try:
			exec("from plugins."+plugin+" import PluginCore as PluginCore")
			print("      - Plugin Core present")
		except:
			print("      - ERROR: Plugin has no Core module")
			sys.exit(0)
		
	def parse_command_line(self):
		"""Parsing command line arguments"""
		log.info("Parsing command-line arguments:")

		parser = argparse.ArgumentParser()

		parser.add_argument( "-l", "--list", action="store_true", dest="plugin_list", default=False, help="List available plugins and workflows")
		parser.add_argument( "-d", "--dry", action="store_true", dest="dry", default=False, help="Only generate workflow file, do nothing more")		
		parser.add_argument( "-w", "--workflow", action="store", dest="workflow", help="Execute predefined workflow")
		parser.add_argument( "-p", "--plugins", action="store", dest="plugins", nargs='+', help="Execute custom workflow assembled on the command line. You can execute a single plugin by typing '-p pluginname' or a sequence of plugins by typing '-p plugin1 plugin2...'")
		parser.add_argument( "-f", "--files", action="store", dest="filenames", nargs='+', help="Supply one or a list of files as input for the plugin(sequence)")
		parser.add_argument( "-s", "--server", action="store_true", dest="server", default=False, help="Generate html webform from workflow xml file to be used in DART web implementation")

		args = parser.parse_args()

		# Show help if no argument is given
		if not len(sys.argv) > 1:
			parser.print_help(sys.stderr)
			log.error("Wrong number of arguments")
			raise SystemExit()

		self.option_dict = {}
		self.option_dict['workflow'] = args.workflow
		self.option_dict['input'] = None
		self.option_dict['pluginseq'] = args.plugins
		self.option_dict['dry'] = args.dry
		self.option_dict['server'] = args.server

		# Get full path for any file provided
		if args.filenames:
			self.option_dict['input'] = [os.path.abspath(file) for file in args.filenames]
			for file in self.option_dict['input']:
				if not os.path.exists(file):
					log.error("Provided file {} does not exist".format(file))
					raise SystemExit()
	
		if self.option_dict['workflow']:
			log.info("    * The following 3D-DART batch configuration file will be executed:", self.option_dict['workflow'])
		
		if self.option_dict['pluginseq']:
			log.info("    * The following command line plugin sequence will be excecuted:")
			for plugin in self.option_dict['pluginseq']:
				log.info("      - {}".format(plugin))
			
		if self.option_dict['input']:
			log.info("    * The following files will be used as input for the batch sequence:")
			for file in self.option_dict['input']:
				log.info("      - {}".format(file))
		
	def get_plugin_list(self):
		"""List all available plugins in plugin directory and all workflows in workflow directory"""
		log.info("List of available plugins. Execute plugin with option -h/--help for more information")
		log.info("    about the function of the plugin")
		plugindir = os.path.join(self.darth_path, 'plugins')
		os.chdir(plugindir)
		plugins = glob.glob('*.py')
		plugins.remove('__init__.py')
		for plugin in plugins:
			basename, extension = os.path.splitext(plugin)
			log.info("    * {}".format(basename))
		
		log.info("List of available workflows:")
		workflowdir = os.path.join(self.darth_path, 'workflows')
		os.chdir(workflowdir)
		workflows = glob.glob('*.xml')
		for workflow in glob.glob('*.xml'):
			basename, extension = os.path.splitext(workflow)
			print("    * {}".format(basename))
		raise SystemExit()
				
	def WorkflowXML(self):
		
		"""Generate workflow.xml file at startup. This file allways needs to be present. The source
		   can be a pregenerated workflow file or a sequence of plugins supplied on the command line 
		   with the -p or --plugin option."""
		
		if not self.option_dict['pluginseq'] == None:
			
			workflow_dict = self.WorkflowSequence()
			
			print("    * Generate workflow from command line constructed plugin batch sequence")
			print("    * Check if all plugins are present and of a valid type")
			
			plugins = []
			for n in workflow_dict.values(): 
				if not n in plugins: 
					plugins.append(n)
					self._TryPluginImport(n)
		
			print("    * Writing workflow XML file as workflow.xml")
			
			outfile = file('workflow.xml','w')
			outfile.write("""<?xml version="1.0" encoding="iso-8859-1"?>\n""")
			outfile.write("""<main id="DARTworkflow">\n""")
			outfile.write("<meta>\n")
			outfile.write("<name>workflow.xml</name>\n")
			outfile.write("<datetime>"+ctime()+"</datetime>\n")
			outfile.write("</meta>\n")
			for key in workflow_dict:
				outfile.write("<plugin id='"+workflow_dict[key]+"' job='"+str(key)+"'>")
				exec("from plugins import "+workflow_dict[key])
				exec("outfile.write(\n"+workflow_dict[key]+".PluginXML())")
				outfile.write("\n</plugin>\n")
			outfile.write("</main>\n")
			outfile.close()
		
			self.option_dict['workflow'] = 'workflow.xml'
			
			if self.option_dict['dry'] == True:
				print("    * Option 'dry' is True, only write workflow.xml file")
				sys.exit(0)
			
		elif not self.option_dict['workflow'] == None:
			
			print("    * Check if all plugins in pre-supplied 3D-DART workflow", self.option_dict['workflow'], "are present and of a valid type")
		
			if os.path.isfile(self.option_dict['workflow']):
				pass
			elif os.path.isfile(RenameFilepath(self.option_dict['workflow'],path=self.darth_path+"/workflows/",extension=".xml")):
				self.option_dict['workflow'] = RenameFilepath(self.option_dict['workflow'],path=self.darth_path+"/workflows",extension=".xml")
			else:
				print("     * ERROR: the workflow:", self.option_dict['workflow'], "cannot be found")
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

			print("    * Plugin sequence is valid")
			
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
				print("    * For program integrety reasons the FileSelector plugin needs to be present, including it")
				pluginlist[0:0] = ['FileSelector']
			for plugin in pluginlist:
				workflow_dict[str(counter)] = plugin
				counter += 1
		else:
			if pluginlist[0] == "FileSelector":
				pass
			else:
				print("    * For program integrety reasons the FileSelector plugin needs to be present, including it")
				pluginlist[0:0] = ['FileSelector']
			for plugin in pluginlist:
				workflow_dict[str(counter)] = plugin
				counter += 1
		
		return workflow_dict
	
	def MakeWebForm(self):
	
		"""Call DARTserver.py to make a webform from the workflow.xml file"""
		
		if self.option_dict['server'] == True and not self.option_dict['workflow'] == None:
			print("--> Generating DARTserver compatible HTML webform from the workflow xml file")
			server = WebServer()
			server.MakeWebForm(verbose=False,xml=self.option_dict['workflow'])
		
			print("    * Webform generated. Stopping DART")
			sys.exit(0)			
