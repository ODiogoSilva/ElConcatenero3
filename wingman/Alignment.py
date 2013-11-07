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
		input_format, self.sequence_code = self.autofinder (input_alignment)

		# parsing the alignment and getting the basic class attributes
		# Three attributes will be assigned: alignment_storage, model and loci_lengths
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
		
	def iter_taxa (self):
		""" Returns a list with the taxa contained in the alignment """

		taxa = [sp for sp in self.alignment_storage]

		return taxa

	def iter_sequences (self):
		""" Returns a list with the sequences contained in the alignment """

		sequences = [seq for seq in self.alignment_storage]

		return sequences

	def collapse (self):
		""" Collapses equal sequences into haplotypes. This method changes the alignment_storage variable and only returns a dictionary with the correspondance between the haplotypes and the original taxa names  """

		collapsed_dic, correspondance_dic = orderedDict(), orderedDict()
		counter = 1

		for taxa, seq in self.alignment_storage.items():
			if seq in collapsed_dic:
				collapsed_dic[seq].append(taxa)
			else:
				collapsed_dic[seq] = [taxa]

		self.alignment_storage = orderedDict()
		for seq, taxa_list in collapsed_dic.items():
			haplotype = "Hap_%s" % (counter)
			self.alignment_storage[haplotype] = seq
			correspondance_dic[haplotype] = taxa_list
			counter += 1

		return correspondance_dic

	def write_to_file (self, output_format, output_file, seq_space_nex=40, seq_space_phy=30, seq_space_ima2=10, cut_space_nex=50, cut_space_phy=50, cut_space_ima2=8, conversion=None, form="leave", gap="-", missing="n"):
		""" Writes the alignment object into a specified output file, automatically adding the extension, according to the output format """

		# Writes file in phylip format
		if output_format == "phylip":

			out_file = open(output_file+".phy","w")
			out_file.write("%s %s\n" % (len(self.alignment_storage), self.loci_lengths))
			for key, seq in self.alignment_storage.items():
					out_file.write("%s %s\n" % (key[:cut_space_phy].ljust(seq_space_phy),seq))
			if self.loci_range = None:
				partition_file = open(output_file+"_part.File","a")
				for partition,lrange in self.loci_range:
					partition_file.write("%s, %s = %s\n" % (model,partition,lrange))

		# Writes file in nexus format
		if output_format == "nexus":

			out_file = open(output_file+".nex","w")
			
			# This writes the output in interleave format
			if form == "interleave":
				out_file.write("#NEXUS\n\nBegin data;\n\tdimensions ntax=%s nchar=%s ;\n\tformat datatype=%s interleave=yes gap=%s missing=%s ;\n\tmatrix\n" % (len(self.alignment_storage), self.loci_lengths, self.sequence_code, gap, missing))
				counter = 0
				for i in range (90,self.loci_lengths,90):
					for key, seq in self.alignment_storage.items():
						out_file.write("%s %s\n" % (key[:self.cut_space_nex].ljust(self.seq_space_nex),seq[counter:i]))
					else:
						out_file.write("\n")
						counter = i				
				else:
					for key, seq in self.alignment_storage.items():
						out_file.write("%s %s\n" % (key[:self.cut_space_nex].ljust(self.seq_space_nex),seq[i:self.loci_lengths]))
					else:
						out_file.write("\n")
				out_file.write(";\n\tend;")

			# This writes the output in leave format (default)
			else:
				out_file.write("#NEXUS\n\nBegin data;\n\tdimensions ntax=%s nchar=%s ;\n\tformat datatype=%s interleave=no gap=%s missing=%s ;\n\tmatrix\n" % (len(self.alignment_storage), self.loci_lengths, self.sequence_code, self.gap, self.missing))
				for key,seq in self.alignment_storage.items():
					out_file.write("%s %s\n" % (key[:self.cut_space_nex].ljust(self.seq_space_nex),seq))
				out_file.write(";\n\tend;")

		# Writes file in fasta format
		if output_format == "fasta":
			out_file = open(self.output_file+".fas","w")
			for key in self.taxa_order:
				out_file.write(">%s\n%s\n" % (key,alignment_dic[key]))				
				
		out_file.close()




