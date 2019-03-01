#!/usr/bin/env python3

import os
import sys
import glob
import shutil
import re
import copy
import importlib
from time import ctime
from dart.system.xpath import Xpath
from dart.system.XMLwriter import Node
from dart.system.utils import make_backup
from dart.system.IOlib import InputOutputControl

# Logging
import logging
logging.basicConfig(format='%(name)s [%(levelname)s] %(message)s', level=logging.INFO)
log = logging.getLogger("Plugins")


class PluginExecutor:

	"""This class takes care of the actual workflow execution"""
	def __init__(self, opt_dict=None, DARTdir=None):
		self.DARTdir = DARTdir
		self.rundir = ''
		self.opt_dict = opt_dict
		self.plugin_executor()

	def _main_xml_data_andler(self, mainxml):
		"""Gets the list of elements for all the first childs of the root, then
		makes a dictionary with element data, finally prints info and returns metadata dictionary"""
		metadata = {}

		# Get all elements
		mainxml.Evaluate(query={1:{'element':'main','attr':{'id':'DARTworkflow'}},2:{'element':'meta','attr':None},3:{'element':None,'attr':None}})
		for data in mainxml.nodeselection[3]:
			mainxml.getElem(node=data,export='string')
		elementlist = mainxml.result
		mainxml.ClearResult()

		# For element get data
		for element in mainxml.nodeselection[3]:
			mainxml.getData(node=element,export='string')
		elementdata = mainxml.result
		mainxml.ClearResult()

		for element in elementlist:
			metadata[element] = elementdata[elementlist.index(element)]

		"""Get workflowsequence"""
		mainxml.Evaluate(query={1:{'element':'plugin','attr':None}})
		for attribute in mainxml.nodeselection[1]:
			mainxml.getAttr(node=attribute,selection='id',export='string')
		pluginid = mainxml.result
		mainxml.ClearResult()
		for attribute in mainxml.nodeselection[1]:
			mainxml.getAttr(node=attribute,selection='job',export='string')
		jobnr = mainxml.result
		mainxml.ClearResult()

		metadata['workflowsequence'] = {}
		for plugin in range(len(pluginid)):
			metadata['workflowsequence'][jobnr[plugin]] = pluginid[plugin]

		"""Print data as info to user"""
		log.info("Metadata belonging to workflow: {}".format(metadata['name']))
		for tag in metadata:
			if type(metadata[tag]) == type({}):
				log.info("    * Workflow {}:".format(tag))
				flow = metadata[tag]
				for key in sorted(flow.keys()):
					log.info("      {0:.0f} = {1:s}".format(key, flow[key]))
			else:
				log.info("    * Workflow {}:".format(metadata[tag]))

		return metadata

	def _MakeRundir(self, name):

		"""Make a workflow directory were all data will be stored in. Make backup if neccesary"""

		self.rundir = os.path.join(os.getcwd(),name)
		self.rundir = make_backup(self.rundir, report=True)
		log.info("Create run directory: {}".format(os.path.basename(self.rundir)))
		os.mkdir(self.rundir)

	def _Executor(self, plugin, mainxml, step):

		metadict = self._MetadataHandler(mainxml,plugin,str(step))
		paramdict = self._ParamDictHandler(mainxml,plugin,str(step))

		if paramdict['useplugin']:
			log.info("Instantiate plugin: {}".format(plugin))
			if step == 1:
				if self.opt_dict['input'] is not None:
					inputlist = self.opt_dict['input']
				else:
					inputlist = self._GetInput(paramdict)
			else:
				inputlist = self._GetInput(paramdict)

			checked = InputOutputControl()
			checked.CheckInput(inputlist,metadict['input'])
			filelist = checked.DictToList()
			self._ChangeDir(plugin, step)
			#exec("from plugins.{} import PluginCore as PluginCore".format(plugin))

			module = importlib.import_module("dart.plugins.{}".format(plugin))
			PluginCore = getattr(module, "PluginCore")
			PluginCore(paramdict, filelist)
			outputlist = checked.CheckOutput(glob.glob('*.*'), metadict['output'])

			os.chdir(self.rundir)

			return outputlist
		else:
			return []

	def _MetadataHandler(self,mainxml,plugin,step):

		"""Get list of tags for all first childs of the root, then make dictionary with tag data, then
		print info and return metadata dictionary"""

		metadata = {}

		"""Get all Elements"""
		mainxml.Evaluate(query={1:{'element':'plugin','attr':{'job':step}},2:{'element':'metadata','attr':None},3:{'element':None,'attr':None}})
		for data in mainxml.nodeselection[3]:
			mainxml.getElem(node=data,export='string')
		elementlist = mainxml.result
		mainxml.ClearResult()

		"""For element get data"""
		for element in mainxml.nodeselection[3]:
			mainxml.getData(node=element,export='string')
		elementdata = mainxml.result
		mainxml.ClearResult()

		for element in elementlist:
			metadata[element] = elementdata[elementlist.index(element)]

		"""Print data as info to user"""
		log.info("Metadata belonging to plugin: {}".format(plugin))
		for data in metadata:
			log.info("    * Plugin {}: {}".format(data,metadata[data]))

		return metadata

	def _ParamDictHandler(self,mainxml,plugin,step):

		"""Returns dictionary with parameter data"""

		"""Get all option attributes"""
		mainxml.Evaluate(query={1:{'element':'plugin','attr':{'job':step}},2:{'element':'parameters','attr':None},3:{'element':'option','attr':None}})
		for data in mainxml.nodeselection[3]:
			mainxml.getAttr(node=data,selection='type',export='string')
		attributelist = mainxml.result
		mainxml.ClearResult()

		"""Get all values belonging to options"""
		paramdict = {}
		for data in mainxml.nodeselection[3]:
			query = Xpath(data.toxml())
			for attribute in attributelist:
				query.Evaluate(query={1:{'element':'option','attr':{'type':attribute}}})
				for node in query.nodeselection[1]:
					query.getData(node=node,export='list')
					if len(query.result) > 0:
						paramdict[attribute] = query.result[0][0]
					else:
						paramdict[attribute] = None
					query.ClearResult()

		"""Print data as info to user"""
		if paramdict['useplugin']:
			log.info("Variables defined for executing the plugin core:")
			for parameter in paramdict:
				log.info("    * Variable {}: {}".format(parameter,paramdict[parameter]))
		else:
			log.info("Use this plugin set to: False")

		return paramdict

	def _GetInput(self, paramdict):

		rawxml = self.xmlroot.rawxml()
		inputlist = []

		if paramdict['inputfrom'] == "self":
			xml = Xpath(rawxml)
			xml.Evaluate(query={1:{'element':'plugin','attr':{'nr':'1'}},2:{'element':None,'attr':None}})
			for elements in xml.nodeselection[2]:
				xml.getElem(node=elements,export='string')
			data = xml.result
			xml.ClearResult()
			return data
		elif paramdict['inputfrom'] == "None":
			return inputlist
		else:
			try:
				inputfrom = paramdict['inputfrom'].split(',')
			except:
				inputfrom = [(paramdict['inputfrom'])]

			data = []
			for n in inputfrom:
				xml = Xpath(rawxml)
				xml.Evaluate(query={1:{'element':'plugin','attr':{'nr':str(int(n))}},2:{'element':'file','attr':None}})
				for elements in xml.nodeselection[2]:
					xml.getData(node=elements,export='string')
				for result in xml.result:
					data.append(result)
			xml.ClearResult()

			return data

	def _ChangeDir(self, plugin, step):

		basedir = os.getcwd()
		dirname = "jobnr"+str(step)+"-"+plugin
		jobdir = os.path.join(basedir,dirname)

		log.info("Create job directory for plugin: {}".format(os.path.basename(jobdir)))
		os.mkdir(jobdir)
		os.chdir(jobdir)

	def _WriteOutput(self, outputlist, plugin, step):

		plugintag = Node("plugin", ID=plugin, nr=str(step))

		if len(outputlist) == 0:
			plugintag += Node("file", "None")
		else:
			for files in outputlist:
				plugintag += Node("file", files)

		self.xmlroot += plugintag

	def _OutputToFile(self):

		"""Write XML list of files for each plugin to file"""

		outfile = file('Filelist.xml','w')
		outfile.write(self.xmlroot.xml())
		outfile.close

	def plugin_executor(self):
		"""Executes the plugins specified in this workflow"""

		# Get meta data from workflow xml file
		mainxml = Xpath(self.opt_dict['workflow'])
		self.maindict = self._main_xml_data_andler(mainxml)

		# Make main workflow directory
		self._MakeRundir(os.path.basename(os.path.splitext(self.maindict['name'])[0]))
		shutil.copy(self.opt_dict['workflow'],self.rundir)
		os.chdir(self.rundir)

		# Write job output to xml file
		self.xmlroot = Node("container", ID="filelist")

		# Executing all plugins
		steps = len(self.maindict['workflowsequence'].keys())
		step = 1
		while step < (steps+1):
			plugin = self.maindict['workflowsequence'][float(step)]
			outputlist = self._Executor(plugin, mainxml, step)
			self._WriteOutput(outputlist, plugin, step)
			step = step+1

		self._OutputToFile()


if __name__ == "__main__":
	"""For testing the script"""
	opt_dict={'workflow':'workflow.xml'}
	PluginExecutor(opt_dict)
