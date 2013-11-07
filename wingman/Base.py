#!/usr/bin/env python3
# -*- coding: utf-8 -*-#
#  
#  Copyright 2012 Unknown <diogo@arch>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

class Base ():

	def autofinder (self, reference_file):
	""" Autodetects the type of file to be parsed. Based on headers """
		autofind = "unknown"
		sequence = ""
		file_handle = open(reference_file,'r')
	    
	    # Skips first empty lines, if any 
		header = file_handle.readline()
		while header.startswith("\n"):
			header = next(file_handle)
	     
	    # Recognition of NEXUS files is based on the existence of the string '#NEXUS' in the first non-empty line
		if header.upper().startswith("#NEXUS"):
			autofind = "nexus"
			for line in file_handle:
				if line.strip().lower() == "matrix":
					sequence = "".join(file_handle.readline().split()[1:]).strip()
					break

	    # Recognition of FASTA files is based on the existence of a ">" character as the first character of a non-empty line     
		elif header.startswith(">"):
			autofind = "fasta"
			for line in file_handle:
				if line.strip() != "" and line.strip()[0] != ">":
					sequence += line.strip()
				elif line.strip() != "" and line.strip()[0] == ">":
					break
	    
	    # Recognition of Phylip files is based on the existence of two integers separated by whitespace on the first non-empy line     
		elif len(header.strip().split()) == 2 and header.strip().split()[0].isdigit() and header.strip().split()[1].isdigit():
			autofind = "phylip"
			sequence = "".join(file_handle.readline().split()[1:]).strip()
		
		# Guessing the genetic code
		code = guess_code (sequence)
		
		return autofind, code

	def partition_format (self, partition_file):
		""" Tries to guess the format of the partition file (Whether it is Nexus of RAxML's) """
		file_handle = open(partition_file)
		
		# Skips first empty lines, if any 
		header = file_handle.readline()
		while header.startswith("\n"):
			header = next(file_handle)
			
		fields = header.split()
		if fields[0].lower() == "charset":
			p_format = "nexus"
		else:
			p_format = "phylip"
		
		return p_format

	def guess_code (self, sequence):
		""" Function that guesses the code of the molecular sequences (i.e., DNA or Protein) based on the first sequence of a reference file """
		sequence = sequence.upper().replace("-","") # Removes gaps from the sequence so that the frequences are not biased
		DNA_count = sequence.count("A") + sequence.count("T") + sequence.count("G") + sequence.count("C") + sequence.count("N")
		DNA_proportion = float(DNA_count)/float(len(sequence))
		if DNA_proportion > 0.9: # The 0.9 cut-off has been effective so far
			code = ("DNA","N")
		else:
			code = ("Protein","X")
		return code	