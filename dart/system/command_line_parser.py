#!/usr/bin/env python3

import os
import glob
import sys
import argparse
from time import ctime
from dart.system.Xpath import Xpath
from dart.system.Utils import rename_file_path
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
log = logging.getLogger("3D-DART")


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

		self.geneterate_workflow_xml()

		#self.make_web_form()

	def parse_command_line(self):
		"""Parsing command line arguments"""
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
			log.info(usage)
			parser.print_help(sys.stderr)
			log.error("Wrong number of arguments")
			raise SystemExit

		self.option_dict = {}
		self.option_dict['workflow'] = args.workflow
		self.option_dict['input'] = None
		self.option_dict['pluginseq'] = args.plugins
		self.option_dict['dry'] = args.dry
		self.option_dict['server'] = args.server

		# List of available plugins
		if args.plugin_list:
			self.get_plugin_list()
			raise SystemExit

		# Get full path for any file provided
		if args.filenames:
			self.option_dict['input'] = [os.path.abspath(file) for file in args.filenames]
			for file in self.option_dict['input']:
				if not os.path.exists(file):
					log.error("Provided file {} does not exist".format(file))
					raise SystemExit

		if self.option_dict['workflow']:
			log.info("    * The following 3D-DART batch configuration file will be executed: {}".format(self.option_dict['workflow']))
			if self.option_dict['pluginseq']:
				log.warning("Incompatible plugin sequence option detected")
				raise SystemExit("Wrong command line")

		if self.option_dict['pluginseq']:
			log.info("    * The following command line plugin sequence will be excecuted:")
			for plugin in self.option_dict['pluginseq']:
				log.info("      - {}".format(plugin))
			if self.option_dict['workflow']:
				log.warning("Incompatible workflow option detected")
				raise SystemExit("Wrong command line")

		if self.option_dict['input']:
			log.info("    * The following files will be used as input for the batch sequence:")
			for file in self.option_dict['input']:
				log.info("      - {}".format(file))

	def get_plugin_list(self):
		"""List all available plugins in plugin directory and all workflows in workflow directory"""
		log.info("List of available plugins. Execute plugin with option -h/--help for more information")
		log.info("    about the function of the plugin:")
		# Plugins
		plugindir = os.path.join(self.darth_path, 'plugins')
		os.chdir(plugindir)
		plugins = glob.glob('*.py')
		plugins.remove('__init__.py')
		for plugin in plugins:
			basename, extension = os.path.splitext(plugin)
			log.info("    * {}".format(basename))
		# Workflow
		log.info("List of available workflows:")
		workflowdir = os.path.join(self.darth_path, 'workflows')
		os.chdir(workflowdir)
		workflows = glob.glob('*.xml')
		for workflow in glob.glob('*.xml'):
			basename, extension = os.path.splitext(workflow)
			log.info("    * {}".format(basename))

	def geneterate_workflow_xml(self):
		"""Generate workflow.xml file at startup.

		This file allways needs to be present. The source can be a pregenerated workflow
		file or a sequence of plugins supplied on the command line with the -p or --plugin option.
		"""
		if self.option_dict['pluginseq']:

			workflow_dict = self.build_workflow_sequence()

			log.info("    * Generate workflow from command line constructed plugin batch sequence")
			log.info("    * Check if all plugins are present and of a valid type")

			plugins = []
			for n in workflow_dict.values():
				if not n in plugins:
					plugins.append(n)

			log.info("    * Writing workflow XML file as workflow.xml")
			with open('workflow.xml', 'w') as outfile:
				outfile.write("""<?xml version="1.0" encoding="iso-8859-1"?>\n""")
				outfile.write("""<main id="DARTworkflow">\n""")
				outfile.write("<meta>\n")
				outfile.write("<name>workflow.xml</name>\n")
				outfile.write("<datetime>{}</datetime>\n".format(ctime()))
				outfile.write("</meta>\n")
				for key in workflow_dict:
					outfile.write("<plugin id='{0}' job='{1}'>".format(workflow_dict[key], str(key)))
					exec("from plugins import {}".format(workflow_dict[key]))
					exec("outfile.write(\n{}.PluginXML())".format(workflow_dict[key]))
					outfile.write("\n</plugin>\n")
				outfile.write("</main>\n")

			self.option_dict['workflow'] = 'workflow.xml'

			if self.option_dict['dry']:
				log.info("    * Option 'dry' is enabled, only write workflow.xml file")
				raise SystemExit

		elif self.option_dict['workflow']:

			log.info("    * Check if all plugins in pre-supplied 3D-DART workflow {} are present and of a valid type".format(self.option_dict['workflow']))

			if os.path.isfile(self.option_dict['workflow']):
				pass
			elif os.path.isfile(rename_file_path(self.option_dict['workflow'], path=self.darth_path+"/workflows/", extension=".xml")):
				self.option_dict['workflow'] = rename_file_path(self.option_dict['workflow'], path=self.darth_path+"/workflows", extension=".xml")
			else:
				log.error("     * ERROR: the workflow: {} cannot be found".format(self.option_dict['workflow']))
				raise SystemExit

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

			log.info("    * Plugin sequence is valid")

	def build_workflow_sequence(self):
		"""Build squence of execution from command line string of plugins"""
		workflow_dict = {}
		counter = 1

		pluginlist =  self.option_dict['pluginseq']
		if len(pluginlist) > 1:
			pluginlist = [pluginlist[len(pluginlist)-1],]+pluginlist[:len(pluginlist)-1]
			if pluginlist[0] == "FileSelector":
				pass
			else:
				log.info("    * For program integrety reasons the FileSelector plugin needs to be present, including it")
				pluginlist[0:0] = ['FileSelector']
			for plugin in pluginlist:
				workflow_dict[str(counter)] = plugin
				counter += 1
		else:
			if pluginlist[0] == "FileSelector":
				pass
			else:
				log.info("    * For program integrety reasons the FileSelector plugin needs to be present, including it")
				pluginlist[0:0] = ['FileSelector']
			for plugin in pluginlist:
				workflow_dict[str(counter)] = plugin
				counter += 1

		return workflow_dict

	def make_web_form(self):
		"""Call DARTserver.py to make a webform from the workflow.xml file"""
		if self.option_dict['server'] and self.option_dict['workflow']:
			log.info("Generating DARTserver compatible HTML webform from the workflow xml file")
			server = WebServer()
			server.MakeWebForm(verbose=False,xml=self.option_dict['workflow'])
			log.info("    * Webform generated. Stopping DART")
			raise SystemExit
