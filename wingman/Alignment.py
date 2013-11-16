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
import re

class Alignment (Base):

	def __init__ (self, input_alignment,input_format=None,model_list=None, alignment_name=None, loci_ranges=None):
		""" The basic Alignment class requires only an alignment file and returns an Alignment object. In case the class is initialized with a dictionary object, the input_format, model_list, alignment_name and loci_ranges arguments can be used to provide complementary information for the class. However, if the class is not initialized with specific values for these arguments, they can be latter set using the _set_format and _set_model functions 

			The loci_ranges argument is only relevant when an Alignment object is initialized from a concatenated data set, in which case it is relevant to incorporate this information in the object"""

		# In case the class is initialized with an input file name
		if type(input_alignment) is str:

			self.input_alignment = input_alignment
			# Get alignment format and code. Sequence code is a tuple of (DNA, N) or (Protein, X)
			self.input_format, self.sequence_code = self.autofinder (input_alignment)

			# In case the input format is specified, overwrite the attribute
			if input_format != None:
				self.input_format = input_format

			# parsing the alignment and getting the basic class attributes
			# Three attributes will be assigned: alignment, model and locus_length
			self.read_alignment (input_alignment, self.input_format)

		# In case the class is initialized with a dictionay object
		elif type(input_alignment) is OrderedDict:

			self.input_alignment = alignment_name # The name of the alignment (str)
			self._init_dicObj(input_alignment) # Gets several attributes from the dictionary alignment 
			self.input_format = input_format # The input format of the alignment (str)
			self.model = model_list # A list containing the alignment model(s) (list)
			self.loci_ranges = loci_ranges # A list containing the ranges of the alignment, in case it's a concatenation

	def _set_loci_ranges (self, loci_list):
		""" Use this function to mannyally set the list with the loci ranges """

		self.loci_ranges = loci_list

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
		self.locus_length = len(list(dictionary_obj.values())[0])

		
	def read_alignment (self, input_alignment, alignment_format, size_check=True):
		""" The read_alignment method is run when the class is initialized to parse an alignment an set all the basic attributes of the class.

		The 'alignment' variable contains an ordered dictionary with the taxa names as keys and sequences as values
		The 'model' is an non essential variable that contains a string with a substitution model of the alignment. This only applies to Nexus input formats, as it is the only supported format that contains such information 
		The 'locus_length' variable contains a int value with the length of the current alignment """
		
		self.alignment = OrderedDict() # Storage taxa names and corresponding sequences in an ordered Dictionary
		self.model = [] # Only applies for nexus format. It stores any potential substitution model at the end of the file
		
		file_handle = open(input_alignment)

		# PARSING PHYLIP FORMAT
		if alignment_format == "phylip":
			header = file_handle.readline().split() # Get the number of taxa and sequence length from the file header
			self.locus_length = int(header[1])
			for line in file_handle:
				if line != "":
					taxa = line.split()[0].replace(" ","")
					taxa = self.rm_illegal(taxa)
					sequence = line.split()[1].strip()
					self.alignment[taxa] = sequence
					
					## TO DO: Read phylip interleave
			
		# PARSING FASTA FORMAT
		elif alignment_format == "fasta":
			for line in file_handle:
				if line.strip().startswith(">"):
					taxa = line[1:].strip().replace(" ","_")
					taxa = self.rm_illegal(taxa)
					self.alignment[taxa] = ""
				elif line.strip() != "":
					self.alignment[taxa] += line.strip()
			self.locus_length = len(list(self.alignment.values())[0])
			
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
					if taxa in self.alignment: # This accomodates for the interleave format
						self.alignment[taxa] += "".join(line.strip().split()[1:])
					else:
						self.alignment[taxa] = "".join(line.strip().split()[1:])
						
				# This bit of code will extract a potential substitution model from the file
				elif counter == 2 and line.lower().strip().startswith("lset"):
					self.model.append(line.strip())
				elif counter == 2 and line.lower().strip().startswith("prset"):
					self.model.append(line.strip())

			self.locus_length = len(list(self.alignment.values())[0])
		
		# Checks the size consistency of the alignment
		if size_check == True:
			self.check_sizes (self.alignment, input_alignment)
		
		# Checks for duplicate taxa
		if len(list(self.alignment)) != len(set(list(self.alignment))):
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

	def remove_taxa (self, taxa_list):
		""" Removes specified taxa from the alignment. taxa_list must be a list """

		new_alignment = OrderedDict()

		for taxa,seq in self.alignment.items():
			if taxa not in taxa_list:
				new_alignment[taxa] = seq

		self.alignment = new_alignment

	def collapse (self, write_haplotypes=True, haplotypes_file=None):
		""" Collapses equal sequences into haplotypes. This method changes the alignment variable and only returns a dictionary with the correspondance between the haplotypes and the original taxa names  """

		collapsed_dic, correspondance_dic = OrderedDict(), OrderedDict()
		counter = 1

		for taxa, seq in self.alignment.items():
			if seq in collapsed_dic:
				collapsed_dic[seq].append(taxa)
			else:
				collapsed_dic[seq] = [taxa]

		self.alignment = OrderedDict()
		for seq, taxa_list in collapsed_dic.items():
			haplotype = "Hap_%s" % (counter)
			self.alignment[haplotype] = seq
			correspondance_dic[haplotype] = taxa_list
			counter += 1

		if write_haplotypes == True:
			# If no output file for the haplotype correspondance is provided, use the input alignment name as reference
			if haplotypes_file == None:
				haplotypes_file = self.input_alignment.split(".")[0]
			self._write_loci_correspondance(correspondance_dic, haplotypes_file)

	def _write_loci_correspondance (self, dic_obj, output_file):
		""" This function supports the collapse method by writing the correspondance between the unique haplotypes and the loci into a new file """

		output_handle = open(output_file+".haplotypes","w")

		for haplotype, taxa_list in dic_obj.items():
			output_handle.write("%s: %s\n" % (haplotype, "; ".join(taxa_list)))

		output_handle.close()

	def reverse_concatenate (self, partition_obj):
		""" This function divides a concatenated file according previously defined partitions and returns an AlignmentList object """

		alignment_list, models, names = [], [], []

		for model, name, part_range in partition_obj.partitions:
			partition_dic = OrderedDict()
			for taxon, seq in self.alignment.items():
				sub_seq = seq[int(part_range[0])-1:int(part_range[1])-1]
				if sub_seq.replace(self.sequence_code[1],"") != "":
					partition_dic[taxon] = sub_seq
			alignment_list.append(partition_dic)
			models.append(model)
			names.append(name)

		alignmentlist_obj = AlignmentList (alignment_list, model_list=models, name_list=names)

		return alignmentlist_obj

	def code_gaps (self):
		""" This method codes gaps present in the alignment in binary format, according to the method of Simmons and Ochoterena (2000), to be read by phylogenetic programs such as MrBayes. The resultant alignment, however, can only be outputed in the Nexus format """

		def gap_listing (sequence,gap_symbol="-"):
			""" Function that parses a sequence string and returns the position of indel events. The returned list is composed of tuples with the span of each indel """
			gap = "%s+" % (gap_symbol)
			span_regex = ""
			gap_list,seq_start = [],0
			while span_regex != None:
				span_regex = re.search(gap,sequence)
				if span_regex != None and seq_start == 0:
					gap_list.append(span_regex.span())
					sequence = sequence[span_regex.span()[1]+1:]
					seq_start = span_regex.span()[1]+1
				elif span_regex != None and seq_start != 0:
					gap_list.append((span_regex.span()[0]+seq_start,span_regex.span()[1]+seq_start))
					sequence = sequence[span_regex.span()[1]+1:]
					seq_start += span_regex.span()[1]+1
			return gap_list

		def gap_binary_generator (sequence,gap_list):
			""" This function contains the algorithm to construct the binary state block for the indel events """
			for cur_gap in gap_list:
				cur_gap_start,cur_gap_end = cur_gap
				if sequence[cur_gap_start:cur_gap_end] == "-"*(cur_gap_end - cur_gap_start) and sequence[cur_gap_start-1] != "-" and sequence[cur_gap_end] != "-":
					sequence += "1"
				elif sequence[cur_gap_start:cur_gap_end] == "-"*(cur_gap_end - cur_gap_start):
					if sequence[cur_gap_start-1] == "-" or sequence[cur_gap_end] == "-":
						sequence += "-"
				elif sequence[cur_gap_start:cur_gap_end] != "-"*(cur_gap_end - cur_gap_start):
					sequence += "0"
			return sequence

		complete_gap_list = []

		# Get the complete list of unique gap positions in the alignment
		for taxa, seq in self.alignment.items():

			current_list = gap_listing (seq)
			complete_gap_list += [gap for gap in current_list if gap not in complete_gap_list]

		# This will add the binary matrix of the unique gaps listed at the end of each alignment sequence
		for taxa, seq in self.alignment.items():
			self.alignment[taxa] = gap_binary_generator (seq, complete_gap_list)

		self.restriction_range = "%s-%s" % (int(self.locus_length), len(complete_gap_list) + int(self.locus_length) - 1)


	def write_to_file (self, output_format, output_file, new_alignment = None, seq_space_nex=40, seq_space_phy=30, seq_space_ima2=10, cut_space_nex=50, cut_space_phy=50, cut_space_ima2=8, form="leave", gap="-", missing="n", model_phylip="LG", model_list=[], outgroup_list=None):
		""" Writes the alignment object into a specified output file, automatically adding the extension, according to the output format 

		This function supports the writting of both converted (no partitions) and concatenated (partitioned files). The choice of this modes is determined by the presence or absence of the loci_range attribute of the object. If its None, there are no partitions and no partitions files will be created. If there are partitions, then the appropriate partitions will be written.

		The outgroup_list argument is used only for Nexus output format and consists in writing a line defining the outgroup. This may be usefull for analyses with MrBayes or other software that may require outgrups"""

		# If this function is called in the AlignmentList class, there may be a need to specify a new alignment dictionary, such as a concatenated one
		if new_alignment != None:
			alignment = new_alignment
		else:
			alignment = self.alignment

		# Checks if there is any other format besides Nexus if the alignment's gap have been coded
		try:
			self.restriction_range
			if output_format != ["nexus"]:
				print ("OutputFormatError: Alignments with gaps coded can only be written in Nexus format")
				return 0
		except:
			pass

		# Writes file in phylip format
		if "phylip" in output_format:

			out_file = open(output_file+".phy","w")
			out_file.write("%s %s\n" % (len(alignment), self.locus_length))
			for key, seq in alignment.items():
					out_file.write("%s %s\n" % (key[:cut_space_phy].ljust(seq_space_phy),seq))

			# In case there is a concatenated alignment being written
			try:
				self.loci_ranges
				partition_file = open(output_file+"_part.File","a")
				for partition,lrange in self.loci_ranges:
					partition_file.write("%s, %s = %s\n" % (model_phylip,partition,lrange))
			except: pass

			out_file.close()

		if "mcmctree" in output_format:

			out_file = open(output_file+".phy", "w")
			taxa_number = len(self.alignment)

			for element in self.loci_ranges:
				partition_range = [int(x) for x in element[1].split("-")]
				out_file.write("%s %s\n" % (taxa_number, (int(partition_range[1])-int(partition_range[0]))))
				for taxon, sequence in self.alignment.items():
					out_file.write("%s  %s\n" % (taxon[:cut_space_phy].ljust(seq_space_phy),sequence[(int(partition_range[0])-1):(int(partition_range[1])-1)]))

			out_file.close()

		# Writes file in nexus format
		if "nexus" in output_format:

			out_file = open(output_file+".nex","w")
			
			# This writes the output in interleave format
			if form == "interleave":
				try:
					self.restriction_range
					out_file.write("#NEXUS\n\nBegin data;\n\tdimensions ntax=%s nchar=%s ;\n\tformat datatype=mixed(%s:1-%s, restriction:%s) interleave=yes gap=%s missing=%s ;\n\tmatrix\n" % (len(alignment), self.locus_length, self.sequence_code[0], self.locus_length-1, self.restriction_range, gap, missing))
				except:
					out_file.write("#NEXUS\n\nBegin data;\n\tdimensions ntax=%s nchar=%s ;\n\tformat datatype=%s interleave=yes gap=%s missing=%s ;\n\tmatrix\n" % (len(alignment), self.locus_length, self.sequence_code[0], gap, missing))
				counter = 0
				for i in range (90,self.locus_length,90):
					for key, seq in alignment.items():
						out_file.write("%s %s\n" % (key[:cut_space_nex].ljust(seq_space_nex),seq[counter:i]))
					else:
						out_file.write("\n")
						counter = i				
				else:
					for key, seq in alignment.items():
						out_file.write("%s %s\n" % (key[:cut_space_nex].ljust(seq_space_nex),seq[i:self.locus_length]))
					else:
						out_file.write("\n")
				out_file.write(";\n\tend;")
				
			# This writes the output in leave format (default)
			else:
				try:
					self.restriction_range
					out_file.write("#NEXUS\n\nBegin data;\n\tdimensions ntax=%s nchar=%s ;\n\tformat datatype=mixed(%s:1-%s, restriction:%s) interleave=yes gap=%s missing=%s ;\n\tmatrix\n" % (len(alignment), self.locus_length, self.sequence_code[0], self.locus_length-1, self.restriction_range, gap, missing))
				except:
					out_file.write("#NEXUS\n\nBegin data;\n\tdimensions ntax=%s nchar=%s ;\n\tformat datatype=%s interleave=no gap=%s missing=%s ;\n\tmatrix\n" % (len(alignment), self.locus_length, self.sequence_code[0], gap, missing))
				for key,seq in alignment.items():
					out_file.write("%s %s\n" % (key[:cut_space_nex].ljust(seq_space_nex),seq))
				out_file.write(";\n\tend;")


			try:
				self.loci_ranges
				out_file.write("\nbegin mrbayes;\n")
				for partition,lrange in self.loci_ranges:
					out_file.write("\tcharset %s = %s;\n" % (partition,lrange))
				out_file.write("\tpartition part = %s: %s;\n\tset partition=part;\nend;\n" % (len(self.loci_ranges),", ".join([part[0] for part in self.loci_ranges])))
			except:pass

			# In case outgroup taxa are specified
			if outgroup_list != None:

				compliant_outgroups = [taxon for taxon in outgroup_list if taxon in self.iter_sequences()] # This assures that only the outgroups present in the current file are written
				if compliant_outgroups != []:
					out_file.write("\nbegin mrbayes;\n\toutgroup %s\nend;\n" % (" ".join(compliant_outgroups)))

				# Concatenates the substitution models of the individual partitions
				if self.model:
					loci_number = 1
					out_file.write("begin mrbayes;\n")
					for model in self.model:
						m1 = model[0].split()
						m2 = model[1].split()
						m1_final = m1[0]+" applyto=("+str(loci_number)+") "+" ".join(m1[1:])
						m2_final = m2[0]+" applyto=("+str(loci_number)+") "+" ".join(m2[1:])
						out_file.write("\t%s\n\t%s\n" % (m1_final, m2_final))
						loci_number += 1
					out_file.write("end;\n")

			out_file.close()

		# Writes file in fasta format
		if "fasta" in output_format:
			out_file = open(output_file+".fas","w")
			for key,seq in self.alignment.items():
				out_file.write(">%s\n%s\n" % (key,seq))				

			out_file.close()


class AlignmentList (Alignment, Base):
	""" At the most basic instance, this class contains a list of Alignment objects upon which several methods can be applied. It only requires either a list of alignment files or .

		It inherits methods from Base and Alignment classes for the write_to_file methods """

	def __init__ (self, alignment_list, model_list=None, name_list=None, verbose=True):

		self.alignment_object_list = []

		if type(alignment_list[0]) is str:

			for alignment in alignment_list:

				if verbose == True:
					print ("\rParsing file %s out of %s" % (alignment_list.index(alignment), len(alignment_list)), end="")

				alignment_object = Alignment(alignment)
				self.alignment_object_list.append(alignment_object)

		elif type(alignment_list[0]) is OrderedDict:

			for alignment, model, name in zip(alignment_list, model_list, name_list):

				alignment_object = Alignment(alignment, model_list=[model], alignment_name=name)
				self.alignment_object_list.append(alignment_object)


	def _get_format (self):
		""" Gets the input format of the first alignment in the list """

		return self.alignment_object_list[0].input_format

	def concatenate (self, progress_stat=True):
		""" The concatenate method will concatenate the multiple sequence alignments and create several attributes 

		This method sets the first three variables below and the concatenation variable containing the dict object"""

		self.loci_lengths = [] # Saves the sequence lengths of the 
		self.loci_range = [] # Saves the loci names as keys and their range as values
		self.models = [] # Saves the substitution models for each one
		
		for alignment_object in self.alignment_object_list:
			
			# When set to True, this statement produces a progress status on the terminal
			if progress_stat == True: 
				print ("\rConcatenating file %s out of %s" % (self.alignment_object_list.index(alignment_object)+1,len(self.alignment_object_list)),end="")

			# Setting the missing data symbol
			missing = alignment_object.sequence_code[1]

			# If input format is nexus, save the substution model, if any
			if alignment_object.input_format == "nexus" and alignment_object.model != []:
				self.models.append(alignment_object.model)

			# Algorithm that fills absent taxa with missing data
			if self.loci_lengths == []:
				self.concatenation = alignment_object.alignment # Create the main alignment dictionary from the first current alignment and never visit this statement again
				self.loci_lengths.append(alignment_object.locus_length)
				self.loci_range.append((alignment_object.input_alignment.split(".")[0],"1-%s" % (alignment_object.locus_length))) # Saving the range for the first loci
	
			else:
				for taxa, sequence in alignment_object.alignment.items(): 
					if taxa in self.concatenation: 
						self.concatenation[taxa] += sequence # Append the sequence from the current alignment to the respective taxa in the main alignment
					elif taxa not in self.concatenation:
						self.concatenation[taxa] = missing*sum(self.loci_lengths)+sequence # If the taxa does not yet exist in the main alignment, create the new entry with a sequence of 'n' characters of the same size as the length of the missed loci and the sequence from the current alignment

				# Saving the range for the subsequent loci
				self.loci_range.append((alignment_object.input_alignment.split(".")[0],"%s-%s" % (sum(self.loci_lengths)+1, sum(self.loci_lengths)+alignment_object.locus_length)))
				self.loci_lengths.append(alignment_object.locus_length)
				
				# Check if any taxa from the main alignment are missing from the current alignment. If yes, fill them with 'n'
				for taxa in self.concatenation:
					if taxa not in alignment_object.alignment: 
						self.concatenation[taxa] += missing*alignment_object.locus_length


		concatenated_alignment = Alignment(self.concatenation, input_format=self._get_format(),model_list=self.models, loci_ranges=self.loci_range)
		return concatenated_alignment

	def iter_alignment_dic (self):

		return [alignment.alignment for alignment in self.alignment_object_list]

	def iter_alignment_obj (self):

		return [alignment for alignment in self.alignment_object_list]

	def write_to_file (self, output_format, form="leave",outgroup_list=[]):
		""" This method writes a list of alignment objects or a concatenated alignment into a file """

		for alignment_obj in self.alignment_object_list:
			output_file_name = alignment_obj.input_alignment.split(".")[0]
			alignment_obj.write_to_file(output_format, output_file=output_file_name, form=form, outgroup_list=outgroup_list)

