#!/usr/bin/env python3

import os


def make_backup(infile, report=False):
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


def file_root_rename(infile, extension, basename):
	"""Renames a file but preserves the extension"""
	outfile = basename + extension
	os.rename(infile, outfile)	


def transform_dash(n):
	"""Transforms dashes often found in 3DNA tables to floats"""
	if n == '---' or n == '----':
		return float(0)
	else:
		return float(n)


def rename_file_path(input_file, path=None, basename=None, extension=None):
	"""Renames a filepath in several ways: path, basename and/or extension"""
	orpath = os.path.dirname(input_file)
	orbasename = os.path.splitext(os.path.basename(input_file))[0] 
	orextension = os.path.splitext(os.path.basename(input_file))[1] 
	
	newfile = ""
	if not path:
		# Set path
		newfile = orpath + "/"
	else:
		newfile = path + "/"
	if not basename:
		# Set file basename
		newfile = newfile + orbasename
	else:
		newfile = newfile + basename
	if not extension:
		# Set file extension
		newfile = newfile + orextension
	else:
		newfile = newfile + extension
	return newfile	
