#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
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

from wingman.Base import Base
from collections import OrderedDict

class Alignment (Base):

	def __init__ (self, input_alignment):
		""" The basic Alignment class requires only an alignment file and returns an Alignment object """

		self.input_alignment = input_alignment

		# Get alignment format and code
		input_format, sequence_code = self.autofinder (input_alignment)

		# parsing the alignment and getting the basic class attributes
		self.read_alignment (input_alignment, input_format)
		
	def read_alignment (self, input_alignment, alignment_format, size_check=True):
		""" The read_alignment method is run when the class is initialized to parse an alignment an set all the basic attributes of the class.

		The 'alignment_storage' variable contains an ordered dictionary with the taxa names as keys and sequences as values
		The 'model' is an non essential variable that contains a string with a substitution model of the alignment. This only applies to Nexus input formats, as it is the only supported format that contains such information 
		The 'loci_lengths' variable contains a int value with the length of the current alignment """
		
		self.alignment_storage = OrderedDict() # Storage taxa names and corresponding sequences in an ordered Dictionary
		self.model = [] # Only applies for nexus format. It stores any potential substitution model at the end of the file
		
		file_handle = open(input_alignment)
		
		# PARSING PHYLIP FORMAT
		if alignment_format == "phylip":
			header = file_handle.readline().split() # Get the number of taxa and sequence length from the file header
			self.loci_lengths = int(header[1])
			for line in file_handle:
				if line != "":
					taxa = line.split()[0].replace(" ","")
					taxa = self.rm_illegal(taxa)
					taxa_order.append(taxa)
					sequence = line.split()[1].strip()
					self.alignment_storage[taxa] = sequence
					
					## TO DO: Read phylip interleave
			
		# PARSING FASTA FORMAT
		elif alignment_format == "fasta":
			for line in file_handle:
				if line.strip().startswith(">"):
					taxa = line[1:].strip().replace(" ","_")
					taxa = self.rm_illegal(taxa)
					taxa_order.append(taxa)
					self.alignment_storage[taxa] = ""
				elif line.strip() != "":
					self.alignment_storage[taxa] += line.strip()
			self.loci_lengths = len(list(self.alignment_storage.values())[0])
			
		# PARSING NEXUS FORMAT
		elif alignment_format == "nexus":
			counter = 0
			for line in file_handle:
				if line.strip().lower() == "matrix" and counter == 0: # Skips the nexus header
					counter = 1
				elif line.strip() == ";" and counter == 1: # Stop parser here
					counter = 2
				elif line.strip() != "" and counter == 1: # Start parsing here
					taxa = line.strip().split()[0].replace(" ","")
					taxa = self.rm_illegal(taxa)
					if taxa not in taxa_order: # Prevents duplications in the list when the input format is interleave
						taxa_order.append(taxa)
					if taxa in self.alignment_storage: # This accomodates for the interleave format
						self.alignment_storage[taxa] += "".join(line.strip().split()[1:])
					else:
						self.alignment_storage[taxa] = "".join(line.strip().split()[1:])
						
				# This bit of code will extract a potential substitution model from the file
				elif counter == 2 and line.lower().strip().startswith("lset"):
					self.model.append(line.strip())
				elif counter == 2 and line.lower().strip().startswith("prset"):
					self.model.append(line.strip())

			self.loci_lengths = len(list(self.alignment_storage.values())[0])
		
		# Checks the size consistency of the alignment
		if size_check == True:
			self.check_sizes (self.alignment_storage, input_alignment)
		
		# Checks for duplicate taxa
		if len(taxa_order) != len(set(taxa_order)):
			taxa = self.duplicate_taxa(taxa_order)
			print ("WARNING: Duplicated taxa have been found in file %s (%s). Please correct this problem and re-run the program\n" %(input_alignment,", ".join(taxa)))
			raise SystemExit
		
		return 0

	def iter_taxa (self):

		taxa = [sp for sp in self.alignment_storage]

		return taxa

	def iter_sequences (self):

		sequences = [seq for seq in self.alignment_storage]

		return sequences
