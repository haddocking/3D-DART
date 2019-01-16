#!/usr/bin/env python3

"""
Parses all information from the workflow XML file as a webform to edit.
Retrieved webform information is saved as an updated XML file and executed as DART workflow.
"""

import cgi
import os
import sys
import shutil
import glob
import time
import subprocess
from dart.system.xpath import Xpath
from dart.system.constants import *


class WebServer:

	"""Make a HTML webform from the data in the workflow XML file. Retrieve data from
	   the webform and make a new (custom) workflow.xml file to be excecuted"""
	
	def __init__(self, DARTDIR=None, remote_env=None):
	
		self.DARTDIR = DARTDIR
		self.jobid = str(int(time.time()))
		self.filestring = "" 
		self.error = ""
	
		self.metadata = {}
		self.pluginmeta = {}
		self.pluginoptions = {}
		self.pluginform = {}
		self.plugindefault = {}
		self.plugintext = {}
		
		self.formdata = {}
		
		# Log webserver call to DARTserver.log
		log = open(self.DARTDIR+'/server-tmp/DARTserver.log','a')
		log.write("* DART server call created on %s with job ID %s with data:\n" % (time.ctime(),self.jobid))
		for n in remote_env:
			log.write("  - %s" % n)
		log.close()

		# Account for the server call in the WeNMR SSO accounting DB
		###sso = ssoxs_connect('3d_dart')
		###sso.accounting(ip=os.environ['REMOTE_ADDR'], status=0, jid=self.jobid)
			
	def _MainXMLdataHandler(self):
	
		"""Make a dictionary of workflow metadata"""
		
		#Get all elements
		self.xml.Evaluate(query={1:{'element':'main','attr':{'id':'DARTworkflow'}},2:{'element':'meta','attr':None},3:{'element':None,'attr':None}})
		for data in self.xml.nodeselection[3]:
			self.xml.getElem(node=data,export='string')					
		elementlist = self.xml.result
		self.xml.ClearResult()
		
		#For element get data
		for element in self.xml.nodeselection[3]:
			self.xml.getData(node=element,export='string')			
		elementdata = self.xml.result
		self.xml.ClearResult()

		for element in elementlist:
			self.metadata[element] = elementdata[elementlist.index(element)]

		#Get workflowsequence
		self.xml.Evaluate(query={1:{'element':'plugin','attr':None}})
		for attribute in self.xml.nodeselection[1]:
			self.xml.getAttr(node=attribute,selection='id',export='string')
		pluginid = self.xml.result
		self.xml.ClearResult()
		for attribute in self.xml.nodeselection[1]:
			self.xml.getAttr(node=attribute,selection='job',export='string')
		jobnr = self.xml.result
		self.xml.ClearResult()
		
		self.metadata['workflowsequence'] = {}
		for plugin in range(len(pluginid)):
			self.metadata['workflowsequence'][int(jobnr[plugin])] = pluginid[plugin]
		
	def _PluginXMLdataHandler(self):
		
		"""Make a dictionary of all plugin data"""	
	
		for job in self.metadata['workflowsequence']:
	
			plugindata = {}
			
			#Get all Elements
			self.xml.Evaluate(query={1:{'element':'plugin','attr':{'job':str(job)}},2:{'element':'metadata','attr':None},3:{'element':None,'attr':None}})
			for data in self.xml.nodeselection[3]:
				self.xml.getElem(node=data,export='string')					
			elementlist = self.xml.result
			self.xml.ClearResult()
			
			#For element get data
			for element in self.xml.nodeselection[3]:
				self.xml.getData(node=element,export='string')			
			elementdata = self.xml.result
			self.xml.ClearResult()
			
			for element in elementlist:
				plugindata[element] = elementdata[elementlist.index(element)]
	
			#Get all option attributes
			self.xml.Evaluate(query={1:{'element':'plugin','attr':{'job':str(job)}},2:{'element':'parameters','attr':None},3:{'element':'option','attr':None}})
			for data in self.xml.nodeselection[3]:
				self.xml.getAttr(node=data,selection='type',export='string')					
			optionlist = self.xml.result
			self.xml.ClearResult()
			
			options = {}
			formdata = {}
			default = {}
			text = {}
			#Get all values belonging to options
			for data in self.xml.nodeselection[3]:
				query = Xpath(data.toxml())
				for option in optionlist:
					query.Evaluate(query={1:{'element':'option','attr':{'type':option}}})
					for node in query.nodeselection[1]:
						query.getData(node=node,export='list')
						if len(query.result) > 0:
							if query.result[0][0] == None:
								options[option] = ''
							else:
								options[option] = query.result[0][0]
						else:
							options[option] = ''	
						query.ClearResult()
						query.getAttr(node=node,selection='form',export='string')
						if len(query.result) > 0:
							formdata[option] = query.result[0]
						else:
							formdata[option] = None	
						query.ClearResult()
						query.getAttr(node=node,selection='default',export='string')
						if len(query.result) > 0:
							default[option] = query.result[0]
						else:
							default[option] = None	
						query.ClearResult()
						query.getAttr(node=node,selection='text',export='string')
						if len(query.result) > 0:
							text[option] = query.result[0]
						else:
							text[option] = None	
						query.ClearResult()
					
			self.pluginmeta[job] = plugindata
			self.pluginoptions[job] = options
			self.pluginform[job] = formdata
			self.plugindefault[job] = default
			self.plugintext[job] = text
		
	def _CheckBox(self,options,form,jobnr):
		
		"""Format checkbox"""
		
		if self.pluginoptions[jobnr][options] == True:
			default = "checked=checked"
		else:
			default = ""	
		form.write('  <tr>\n')
		form.write('   <td width="200"><p>%s</p></td>\n' % self.plugintext[jobnr][options])
		form.write('   <td width="600"><input type=checkbox name=%s value=%s %s></td>\n' % (str(jobnr)+"."+options,self.pluginoptions[jobnr][options],default))
		form.write('  </tr>\n')
	
	def _TextField(self,options,form,jobnr):
	
		"""Format text field"""
		
		form.write('  <tr>\n')
		form.write('   <td width="200"><p>%s</p></td>\n' % self.plugintext[jobnr][options])
		form.write('   <td width="600"><input type=text name=%s value=%s></td>\n' % (str(jobnr)+"."+options,self.pluginoptions[jobnr][options]))
		form.write('  </tr>\n')
	
	def _DDmenu(self,options,form,jobnr):
	
		"""Format dropdown menu"""
		
		value = self.plugindefault[jobnr][options].split(',')
		
		form.write('  <tr>\n')
		form.write('   <td width="200"><p>%s</p></td>\n' % self.plugintext[jobnr][options])
		form.write('   <td width="600"><select name=%s>\n' % (str(jobnr)+"."+options))
		for n in value:
			form.write('    <option value=%s>%s</option>\n' % (n,n))
		form.write('       </select></td>\n')	
		form.write('  </tr>\n')
	
	def _FileField(self,options,form,jobnr):
	
		"""Format file field"""
		
		form.write('  <tr>\n')
		form.write('   <td width="200"><p>%s</p></td>\n' % self.plugintext[jobnr][options])
		form.write('   <td width="600"><input type=file name=%s value=%s></td>\n' % (str(jobnr)+"."+options,self.pluginoptions[jobnr][options]))
		form.write('  </tr>\n')
	
	def _HiddenField(self,options,form,jobnr):
	
		"""Format hidden fields"""
		
		form.write('  <tr>\n')
		form.write('   <td width="800"><input type=hidden name=%s value=%s></td>\n' % (str(jobnr)+"."+options,self.pluginoptions[jobnr][options]))
		form.write('  </tr>\n')
		
	def _WriteHTML(self):
	
		if self.verbose == True:
			print("Content-Type: text/html\n\n")
			form = sys.stdout
		else:
			form = file(os.path.splitext(self.metadata['name'])[0]+'.html','w')
		
		form.write("<html>")
		form.write("<body>")
		form.write("<form method=post action=http://localhost/cgi-bin/dart.py enctype=multipart/form-data>\n")
		
		for jobnr in self.metadata['workflowsequence']:
			form.write(" <br/>\n")
			form.write(" <p><strong>%i  %s</strong></p>\n" % (jobnr,self.pluginmeta[jobnr]['name']))
			form.write(" <br/>\n")
			
			form.write(' <table width="800">\n')
			for options in self.pluginoptions[jobnr]:
				if self.pluginform[jobnr][options] == "checkbox":
					self._CheckBox(options,form,jobnr)
				elif self.pluginform[jobnr][options] == "text":
					self._TextField(options,form,jobnr)
				elif self.pluginform[jobnr][options] == "list":
					self._DDmenu(options,form,jobnr)
				elif self.pluginform[jobnr][options] == "file":
					self._FileField(options,form,jobnr)
				elif self.pluginform[jobnr][options] == "hidden":
					self._HiddenField(options,form,jobnr)		
				else:
					pass		
			form.write(" </table>\n")
		form.write(" <br/>\n")		
		form.write(" <input type=reset name=reset>")
		form.write(" <input type=submit name=%s value=submit>" % self.metadata['name'])	
		form.write("</form>\n")
		form.write("</body>\n")
		form.write("</html>\n")
		
		if self.verbose == False:
			form.close()

	def _FormatFormData(self,webform):
		
		"""Make dictionary of webform data"""
		
		for keys in webform:								# Make a plugin centeric dictionary of all key/value pairs
			pl,op = keys.split('.')
			if not self.formdata.has_key(pl): self.formdata[pl] = {}
			self.formdata[pl][op] = webform[keys]
		
		for keys in self.formdata:							# get default values from workflow xml file
			for options in self.formdata[keys]: 
				if self.formdata[keys][options] == 'submit':
					if os.path.isfile("%s/workflows/%s.xml" % (self.DARTDIR, keys)):
						self.xml = Xpath("%s/workflows/%s.xml" % (self.DARTDIR, keys))
						self._MainXMLdataHandler()
						self._PluginXMLdataHandler()
		
		for plugin in self.pluginoptions:					# Check formdata against default workflow. Missing values
			for options in self.pluginoptions[plugin]:		# are set to there default ones.
				if self.formdata.has_key(str(plugin)):
					if not self.formdata[str(plugin)].has_key(options): 
						self.formdata[str(plugin)][options] = self.pluginoptions[plugin][options]
						#if self.pluginform[plugin][options] == 'checkbox':
						#	self.formdata[str(plugin)][options] = False
						#else:
						#	self.formdata[str(plugin)][options] = self.pluginoptions[plugin][options]
				else:
					self.formdata[str(plugin)] = self.pluginoptions[plugin]	# If complete plugin is missing in formdata than
																			# default values are parsed in			
	
	def _GetDirSize(self,pdb):
		
		"""Get the total size of the uploaded directory"""
		
		kb = 0
		for files in pdb:
			kb = kb + os.stat(files)[6]	
		return kb/1024.0
			
	def _ManageUploads(self):
	
		"""Save uploads to file and move to temporary directory"""
		
		for n in self.formdata['1']:
			if n == 'upload':
				if (self.formdata['1'][n]):	
					filename = self.formdata['1'][n]['name'] 					
					upload = open(filename,'w')
					upload.write(self.formdata['1'][n]['file'])
					upload.close()
					if os.path.splitext(filename)[1] == '.zip':
						cmd = 'unzip '+filename
						os.system(cmd)
						if os.path.isdir(os.path.splitext(filename)[0]):
							curdir = os.getcwd()
							os.chdir(curdir+os.path.splitext(filename)[0])
							pdb = glob.glob("*.pdb")
							os.chdir(curdir)
						else:
							pdb = glob.glob("*.pdb")
						if len(pdb) > 0:
							if self._GetDirSize(pdb) > MAXMB:
								self.error = self.error+("Total size of uploaded files exeeds limit of %1.0f MB" % MAXMB)
							else:
								for files in pdb:
									self.filestring = self.filestring+files+" "
								self.formdata['1']['upload'] = self.DARTDIR+'/server-tmp/'+filename
						else:
							self.error = self.error+"No valid upload found\n"	
					else:
						self.formdata['1']['upload'] = self.DARTDIR+'server-tmp/'+filename
						self.filestring = self.filestring+filename					
				else:
					self.formdata['1']['upload'] = None
			
	def _WriteNewXML(self):
	
		"""Write parsed data from webform to new workflow xml file"""
		outfile = open('workflow.xml','w')
		outfile.write('<?xml version="1.0" encoding="iso-8859-1"?>\n')
		outfile.write('<main id="DARTworkflow">\n')
		outfile.write('<meta>\n')
		
		for tag in self.metadata:
			if not tag == 'workflowsequence':
				outfile.write('<%s>%s</%s>\n' % (tag,self.metadata[tag],tag))
			
		outfile.write('</meta>\n')
		
		for job in self.pluginmeta:
			plugin = self.metadata['workflowsequence'][job]
			outfile.write('<plugin id="%s" job="%i">\n' % (plugin,job))
			outfile.write('<metadata>\n')
			for data in self.pluginmeta[job]:
				outfile.write('<%s>%s</%s>\n' % (data,self.pluginmeta[job][data],data))
			outfile.write('</metadata>\n')
			outfile.write('<parameters>\n')
			for options in self.formdata[str(job)]:
				if self.pluginform[job].has_key(options):
					outfile.write('<option type="%s" form="%s">%s</option>\n' % (options,self.pluginform[job][options],self.formdata[str(job)][options]))
			outfile.write('</parameters>\n')
			outfile.write('</plugin>\n')
		
		outfile.write('</main>\n')						

	def MakeWebForm(self,verbose=False,xml=None):
	
		"""Control module for generating a webform from XMLdata"""
		
		self.xml = Xpath(xml)
		self.verbose = verbose
		
		self._MainXMLdataHandler()
		self._PluginXMLdataHandler()
		self._WriteHTML()

	def CleanJobs(self):

		newlines = []
		curtime = time.time()
		readfile = file(self.DARTDIR+'/server-tmp/Joblist.txt','r')
		for line in readfile.readlines():
			line.strip()
			line = line.split()
			if 'CLEANED' in line: newlines.append(line)
			else:
				if curtime - float(line[1]) >= float(CLEANTIME):
					target = os.path.join(self.DARTDIR+'/results/',line[0])
					if os.path.isfile(target):
						os.remove(target)
						line.append('CLEANED')
						newlines.append(line)
					else:
						line.append('CLEANED')
						newlines.append(line)
				else:
					line.append(" ")
					newlines.append(line)		
		readfile.close()

		readfile = file(self.DARTDIR+'/server-tmp/Joblist.txt','w')
		for newline in newlines:
			readfile.write("%s    %s    %s\n" % (newline[0],newline[1],newline[2]))
		readfile.close()

	def RunDART(self, pythondict):
		
		"""This is the control definition of the Webserver class. It is responsible for the conversion
		   of the webform data to a workflow.xml file and first check of all parameters. It sets up a
		   temporary directory for the given job, excecutes the job, checks if the job finished, zip 
		   the job directory and put it on the FTP site, notify the user for the download location and
		   cleans up afterwards""" 
		
		"""Retrieve the data from the webform"""
		self._FormatFormData(pythondict)
	
		"""Prepaire temporary working directory"""
		os.chdir(self.DARTDIR+'/server-tmp/')		# Move to server temporary directory
		os.mkdir('job'+self.jobid)			# Create temporary working directory
		os.chdir('job'+self.jobid)			# Move to the temporary working directory

		"""Write all nessacary files to temporary working directory"""
		self._ManageUploads()
		self._WriteNewXML()
		
		dirname = os.path.splitext(self.metadata['name'])[0]
		
		"""Run DART in server mode"""
		if self.formdata['1']['upload']: 
		  cmd = "%s %s/RunDART.py -w workflow.xml -f %s > dart.out" % (PYTHON, self.DARTDIR, self.filestring)
		else:
		  cmd = "%s %s/RunDART.py -w workflow.xml > dart.out" % (PYTHON, self.DARTDIR)
		os.system(cmd)
		
		"""Rename project directory, compress and move to FTP directory"""
		if os.path.isdir(os.path.splitext(self.metadata['name'])[0]) and os.path.isfile('dart.out'):
			shutil.move('dart.out',dirname)							# Move dart.out file inside job directory
			os.rename(dirname,dirname+self.jobid)						# Rename job directory to include jobid
			shutil.copy(self.DARTDIR+'/server-tmp/readme.txt',dirname+self.jobid)		# Copy a version of the readme file to the job directory
			for n in ['/workflow.xml','/Filelist.xml']:					# Remove the Filelist.xml and workflow.xml file
				if os.path.isfile(dirname+self.jobid+n): 
					os.remove(dirname+self.jobid+n)
			
			package = dirname+self.jobid							# Compress job directory using zip
			cmd1 = ("zip %s -r %s" % (package,package))
			output = subprocess.call(cmd1 , shell=True)
			package = package +".zip"							
			shutil.move(package,self.DARTDIR+'/results/')					# Move compressed job directory to download location
			os.chdir(self.DARTDIR+'/server-tmp/')						# Move back to server temporary directory
			shutil.rmtree('job'+self.jobid)							# Remove temporary job directory
			
			joblist = open(self.DARTDIR+'/server-tmp/Joblist.txt','a')			# Write finished job info to joblist
			joblist.write("%s    %s\n" % (package,self.jobid))
			joblist.close()
			
			if os.path.isfile(SERVERCOUNTFILE):					
				joblist = open(self.DARTDIR+'/server-tmp/Joblist.txt','r')		# Update the number of served requests
				served_requests = open(SERVERCOUNTFILE,'w')		
				served_requests.write("%i" % len(joblist.readlines()))
				served_requests.close()
				joblist.close()
			
			self.CleanJobs()								# Clean the FTP site for jobs older then CLEANTIME
			
			downloadpath = os.path.join(FTP_LOCATION, package)
			##sso = ssoxs_connect('3d_dart')							# WeNMR SSO account, job done		
			##sso.accounting(status=5, jid=self.jobid, url=downloadpath)
			
			return downloadpath								# Report download location to user
		else:
			return ("An error orccured during processing: %s" % self.error) 
	
