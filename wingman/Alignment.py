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
from wingman.ErrorHandling import *
from collections import OrderedDict

class Alignment (Base):

	def __init__ (self, input_alignment,input_format=None,model_list=None):
		""" The basic Alignment class requires only an alignment file and returns an Alignment object. In case the class is initialized with a dictionary object, the input_format and model_list arguments can be used to provide complementary information for the class. However, if the class is not initialized with specific values for these arguments, they can be latter set using the _set_format and _set_model functions """

		# TODO: I need to, somehow, incorporate the concatenated data set as an Alignment object in order to have access to its methods

		self.input_alignment = input_alignment

		# In case the class is initialized with an input file name
		if type(self.input_alignment) is str:

			# Get alignment format and code. Sequence code is a tuple of (DNA, N) or (Protein, X)
			self.input_format, self.sequence_code = self.autofinder (input_alignment)

			# In case the input format is specified, overwrite the attribute
			if input_format != None or input_format != "guess":
				self.input_format = input_format

			# parsing the alignment and getting the basic class attributes
			# Three attributes will be assigned: alignment, model and loci_lengths
			self.read_alignment (input_alignment, self.input_format)

		# In case the class is initialized with a dictionay object
		elif type(self.input_alignment) is dict:

			self._init_dicObj(self.input_alignment)
			self.input_format = input_format
			self.model = model_list

	def _set_format (self, input_format):
		""" Use this function to manually set the input format associated with the Alignment object """

		self.input_format = _set_format

	def _set_model (self, model_list):
		""" Use this function to manyally set the model associated with the Alignment object. Since this object supports concatenated alignments, the model specification must be in list format and the list size must be of the same size of the alignment partitions """

		self.model = model_list

	def _init_dicObj (self, dictionary_obj):
		""" In case the class is initialized with a dictionary as input, this function will retrieve the same information as the read_alignment function would do  """

		self.sequence_code = self.guess_code(list(dictionary_obj.values())[0])
		self.alignment = dictionary_obj
		self.loci_lengths = len(list(dictionary_obj.values())[0])

		
	def read_alignment (self, input_alignment, alignment_format, size_check=True):
		""" The read_alignment method is run when the class is initialized to parse an alignment an set all the basic attributes of the class.

		The 'alignment' variable contains an ordered dictionary with the taxa names as keys and sequences as values
		The 'model' is an non essential variable that contains a string with a substitution model of the alignment. This only applies to Nexus input formats, as it is the only supported format that contains such information 
		The 'loci_lengths' variable contains a int value with the length of the current alignment """
		
		self.alignment = OrderedDict() # Storage taxa names and corresponding sequences in an ordered Dictionary
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
					self.alignment[taxa] = sequence
					
					## TO DO: Read phylip interleave
			
		# PARSING FASTA FORMAT
		elif alignment_format == "fasta":
			for line in file_handle:
				if line.strip().startswith(">"):
					taxa = line[1:].strip().replace(" ","_")
					taxa = self.rm_illegal(taxa)
					taxa_order.append(taxa)
					self.alignment[taxa] = ""
				elif line.strip() != "":
					self.alignment[taxa] += line.strip()
			self.loci_lengths = len(list(self.alignment.values())[0])
			
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
					if taxa in self.alignment: # This accomodates for the interleave format
						self.alignment[taxa] += "".join(line.strip().split()[1:])
					else:
						self.alignment[taxa] = "".join(line.strip().split()[1:])
						
				# This bit of code will extract a potential substitution model from the file
				elif counter == 2 and line.lower().strip().startswith("lset"):
					self.model.append(line.strip())
				elif counter == 2 and line.lower().strip().startswith("prset"):
					self.model.append(line.strip())

			self.loci_lengths = len(list(self.alignment.values())[0])
		
		# Checks the size consistency of the alignment
		if size_check == True:
			self.check_sizes (self.alignment, input_alignment)
		
		# Checks for duplicate taxa
		if len(taxa_order) != len(set(taxa_order)):
			taxa = self.duplicate_taxa(taxa_order)
			print ("WARNING: Duplicated taxa have been found in file %s (%s). Please correct this problem and re-run the program\n" %(input_alignment,", ".join(taxa)))
			raise SystemExit
		
	def iter_taxa (self):
		""" Returns a list with the taxa contained in the alignment """

		taxa = [sp for sp in self.alignment]

		return taxa

	def iter_sequences (self):
		""" Returns a list with the sequences contained in the alignment """

		sequences = [seq for seq in self.alignment]

		return sequences

	def collapse (self):
		""" Collapses equal sequences into haplotypes. This method changes the alignment variable and only returns a dictionary with the correspondance between the haplotypes and the original taxa names  """

		collapsed_dic, correspondance_dic = orderedDict(), orderedDict()
		counter = 1

		for taxa, seq in self.alignment.items():
			if seq in collapsed_dic:
				collapsed_dic[seq].append(taxa)
			else:
				collapsed_dic[seq] = [taxa]

		self.alignment = orderedDict()
		for seq, taxa_list in collapsed_dic.items():
			haplotype = "Hap_%s" % (counter)
			self.alignment[haplotype] = seq
			correspondance_dic[haplotype] = taxa_list
			counter += 1

		return correspondance_dic

	def write_to_file (self, output_format, output_file, new_alignment = None, seq_space_nex=40, seq_space_phy=30, seq_space_ima2=10, cut_space_nex=50, cut_space_phy=50, cut_space_ima2=8, conversion=None, form="leave", gap="-", missing="n", loci_range=[], model_phylip="LG", model_list=[]):
		""" Writes the alignment object into a specified output file, automatically adding the extension, according to the output format """

		# If this function is called in the AlignmentList class, there may be a need to specify a new alignment dictionary, such as a concatenated one
		if new_alignment != None:
			alignment = new_alignment
		else:
			alignment = self.alignment

		# Writes file in phylip format
		if output_format == "phylip":

			out_file = open(output_file+".phy","w")
			out_file.write("%s %s\n" % (len(alignment), self.loci_lengths))
			for key, seq in alignment.items():
					out_file.write("%s %s\n" % (key[:cut_space_phy].ljust(seq_space_phy),seq))

			# In case there is a concatenated alignment being written
			if conversion == None:
				partition_file = open(output_file+"_part.File","a")
				for partition,lrange in loci_range:
					partition_file.write("%s, %s = %s\n" % (model_phylip,partition,lrange))

		# Writes file in nexus format
		if output_format == "nexus":

			out_file = open(output_file+".nex","w")
			
			# This writes the output in interleave format
			if form == "interleave":
				out_file.write("#NEXUS\n\nBegin data;\n\tdimensions ntax=%s nchar=%s ;\n\tformat datatype=%s interleave=yes gap=%s missing=%s ;\n\tmatrix\n" % (len(alignment), self.loci_lengths, self.sequence_code[0], gap, missing))
				counter = 0
				for i in range (90,self.loci_lengths,90):
					for key, seq in alignment.items():
						out_file.write("%s %s\n" % (key[:self.cut_space_nex].ljust(self.seq_space_nex),seq[counter:i]))
					else:
						out_file.write("\n")
						counter = i				
				else:
					for key, seq in alignment.items():
						out_file.write("%s %s\n" % (key[:self.cut_space_nex].ljust(self.seq_space_nex),seq[i:self.loci_lengths]))
					else:
						out_file.write("\n")
				out_file.write(";\n\tend;")

			if conversion == False:
				out_file.write("\nbegin mrbayes;\n")
				for partition,lrange in loci_range:
					out_file.write("\tcharset %s = %s;\n" % (partition,lrange))
				out_file.write("\tpartition part = %s: %s;\n\tset partition=part;\nend;\n" % (len(loci_range),", ".join([part[0] for part in loci_range])))
				
				# Concatenates the substitution models of the individual partitions
				if model_list != []:
					loci_number = 1
					out_file.write("begin mrbayes;\n")
					for model in self.model_list:
						m1 = model[0].split()
						m2 = model[1].split()
						m1_final = m1[0]+" applyto=("+str(loci_number)+") "+" ".join(m1[1:])
						m2_final = m2[0]+" applyto=("+str(loci_number)+") "+" ".join(m2[1:])
						out_file.write("\t%s\n\t%s\n" % (m1_final, m2_final))
						loci_number += 1
					out_file.write("end;\n")

			# This writes the output in leave format (default)
			else:
				out_file.write("#NEXUS\n\nBegin data;\n\tdimensions ntax=%s nchar=%s ;\n\tformat datatype=%s interleave=no gap=%s missing=%s ;\n\tmatrix\n" % (len(alignment), self.loci_lengths, self.sequence_code, self.gap, self.missing))
				for key,seq in alignment.items():
					out_file.write("%s %s\n" % (key[:self.cut_space_nex].ljust(self.seq_space_nex),seq))
				out_file.write(";\n\tend;")

		# Writes file in fasta format
		if output_format == "fasta":
			out_file = open(self.output_file+".fas","w")
			for key in self.taxa_order:
				out_file.write(">%s\n%s\n" % (key,alignment_dic[key]))				

		out_file.close()


class AlignmentList (Alignment, Base):
	""" At the most basic instance, this class contains a list of Alignment objects upon which several methods can be applied. It only requires a list of alignment files.

		It inherits methods from Base and Alignment classes for the write_to_file methods """

	def __init__ (self, alignment_list):

		self.alignment_object_list = []

		for alignment in alignment_list:

			alignment_object = Alignment(alignment)
			self.alignment_object_list.append(alignment_object)

	def _get_format (self):
		""" Gets the input format of the first alignment in the list """

		return self.alignment_object_list[0].input_format

	def concatenate (self, missing="n"):
		""" The concatenate method will concatenate the multiple sequence alignments and create several attributes 

		This method sets the first three variables below and the concatenation variable containing the dict object"""

		self.loci_lengths = [] # Saves the sequence lengths of the 
		self.loci_range = [] # Saves the loci names as keys and their range as values
		self.models = [] # Saves the substitution models for each one
		
		for alignment_object in self.alignment_object_list:
			
			# When set to True, this statement produces a progress status on the terminal
			if progress_stat == True: 
				print ("\rProcessing file %s out of %s" % (self.alignment_object_list.index(alignment_object)+1,len(self.alignment_object_list)),end="")

			# If input format is nexus, save the substution model, if any
			if alignment_object.input_format == "nexus" and alignment_object.model != []:
				self.models.append(alignment_object.model)

			# Algorithm that fills absent taxa with missing data
			if self.loci_lengths == []:
				self.concatenation = alignment_object.alignment # Create the main alignment dictionary from the first current alignment and never visit this statement again
				self.loci_lengths.append(alignment_object.loci_lengths)
				self.loci_range.append((alignment_object.input_alignment.split(".")[0],"1-%s" % (alignment_object.loci_lengths))) # Saving the range for the first loci
	
			else:
				for taxa, sequence in alignment_object.alignment.items(): 
					if taxa in self.concatenation: 
						self.concatenation[taxa] += sequence # Append the sequence from the current alignment to the respective taxa in the main alignment
					elif taxa not in self.concatenation:
						self.concatenation[taxa] = missing*sum(self.loci_lengths)+sequence # If the taxa does not yet exist in the main alignment, create the new entry with a sequence of 'n' characters of the same size as the length of the missed loci and the sequence from the current alignment
						main_taxa_order.append(taxa)

				# Saving the range for the subsequent loci
				self.loci_range.append((alignment_object.input_alignment.split(".")[0],"%s-%s" % (sum(self.loci_lengths)+1, sum(self.loci_lengths)+alignment_object.loci_lengths)))
				self.loci_lengths.append(alignment_object.loci_lengths)
				
				# Check if any taxa from the main alignment are missing from the current alignment. If yes, fill them with 'n'
				for taxa in self.concatenation:
					if taxa not in alignment_object.alignment: 
						self.concatenation[taxa] += missing*alignment_object.loci_lengths
		else:
			print ("\n")

		concatenated_alignment = Alignment(self.concatenation, input_format=self._get_format,model_list=self.models)
		return concatenated_alignment

	def iter_alignment_dic (self):

		return [alignment.alignment for alignment in self.alignment_object_list]

	def iter_alignment_obj (self):

		return [alignment for alignment in self.alignment_object_list]

	def write_to_file (self, output_format, output_file=None):
		""" This method writes a list of alignment objects or a concatenated alignment into a file """

		for alignment_obj in self.alignment_obj_list:
			output_file_name = alignment_obj.input_file.split(".")[0]
			alignment_obj.write_to_file(output_format, output_file=output_file_name)

