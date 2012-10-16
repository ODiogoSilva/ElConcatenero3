#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ElParsito.py v3.0.0-3
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

class SeqUtils ():
	def __init__ (self, missing):
		self.missing = missing
	
	def check_format (self,input_alignment,alignment_format):
		""" This function performs some very basic checks to see if the format of the input file is in accordance to the input file format specified when the script is executed """
		input_handle = open(input_alignment)
		line = input_handle.readline()
		while line.strip() == "":
			line = next(input_handle)
		
		if alignment_format == "fasta":
			if line.strip()[0] != ">":
				print ("File not in Fasta format. First non-empty line of the input file %s does not start with '>'. Please verify the file, or the input format settings\nExiting..." % input_alignment)
				raise SystemExit
		elif alignment_format == "nexus":
			if line.strip().lower() != "#nexus":
				print ("File not in Nexus format. First non-empty line of the input file %s does not start with '#NEXUS'. Please verify the file, or the input format settings\nExiting..." % input_alignment)
				raise SystemExit
		elif alignment_format == "phylip":
			try:
				header = line.strip().split()
				float(header[0])
				float(header[1])
			except:
				print ("File not in correct Phylip format. First non-empty line of the input file %s does not start with two intergers separated by whitespace. Please verify the file, or the input format settings\nExiting..." % input_alignment)
				raise SystemExit
				
	def rm_taxa (self, alignment_dic, taxa_list):
		""" Function that removes specified taxa from the alignment """
		alignment_mod = {}

		for taxa, sequence in alignment_dic.items():
			if taxa not in taxa_list:
				alignment_mod[taxa] = sequence
				self.taxa_order.append(taxa)
		return alignment_mod, self.taxa_order
	
	def pickle_taxa (self, alignment_dic, mode):
		""" Function that exports the list of taxa from an alignment """
		import pickle
		self.taxa_list = []
		if mode == "dump":
			self.taxa_list = [taxa for taxa in alignment_dic.keys()]
			pickle.dump(self.taxa_list, open("taxa_list","wb"))
			print ("Taxa names have been saved in the pickle file 'taxa_list'\nExiting...")
			raise SystemExit
		elif mode == "load":
			self.taxa_list = pickle.load(open("taxa_list","rb"))
		return self.taxa_list
		
	def import_taxa (self, alignment_dic):
		""" Function that imports new taxa. It mainly exists to complete single locus aligments with taxa that are not present in the current alignment but occur in other alignments """
		alignment_len = self.loci_lengths[0]
		for taxa in self.taxa_list:
			if taxa not in alignment_dic:
				alignment_dic[taxa] = self.missing*alignment_len
		return alignment_dic, self.taxa_list

	def read_alignment (self, input_alignment, alignment_format):
		""" ONLY FOR SINGLE FILE/LOCI INPUT: Function that parses an input file alignment and returns a dictionary with the taxa as keys and sequences as values """ 
		
		self.check_format (input_alignment, alignment_format)
		
		def rm_illegal (string):
			""" Function that removes illegal characters from taxa names """
			illegal_chars = [":",",",")","(",";","[","]","'"]
			clean_name = "".join([char for char in string if char not in illegal_chars])
			if string != clean_name:
				print ("\nWARNING: Removed illegal characters from the taxa %s" % string)
			return clean_name

		alignment_storage = {} # Save the taxa and their respective sequences
		self.taxa_order = [] # Save taxa names to maintain initial order
		file_handle = open(input_alignment)
		
		# PARSING PHYLIP FORMAT
		if alignment_format == "phylip":
			header = file_handle.readline().split() # Get the number of taxa and sequence length from the file header
			self.loci_lengths = [int(header[1])]
			for line in file_handle:
				if line != "":
					taxa = line.split()[0].replace(" ","")
					taxa = rm_illegal(taxa)
					self.taxa_order.append(taxa)
					sequence = line.split()[1].strip()
					alignment_storage[taxa] = sequence
			
		# PARSING FASTA FORMAT
		elif alignment_format == "fasta":
			for line in file_handle:
				if line.strip().startswith(">"):
					taxa = line[1:].strip().replace(" ","_")
					taxa = rm_illegal(taxa)
					self.taxa_order.append(taxa)
					alignment_storage[taxa] = ""
				else:
					alignment_storage[taxa] += line.strip()
			self.loci_lengths = [len(list(alignment_storage.values())[0])]
			
		# PARSING NEXUS FORMAT
		elif alignment_format == "nexus":
			counter = 0
			for line in file_handle:
				if line.strip().lower() == "matrix" and counter == 0: # Skips the nexus header
					counter = 1
				elif line.strip() == ";" and counter == 1: # Stop parser here
					counter = 0
				elif line != "" and counter == 1: # Start parsing here
					taxa = line.strip().split()[0].replace(" ","")
					taxa = rm_illegal(taxa)
					self.taxa_order.append(taxa)
					if taxa in alignment_storage: # This accomodates for the interleave format
						alignment_storage[taxa] += "".join(line.strip().split()[1:])
					else:
						alignment_storage[taxa] = "".join(line.strip().split()[1:])
			self.loci_lengths = [len(list(alignment_storage.values())[0])]
			
		self.loci_range = self.loci_lengths
		self.check_sizes (alignment_storage, input_alignment)
		return (alignment_storage, self.taxa_order, self.loci_lengths, self.loci_range)
		
	def read_alignments (self, input_list, alignment_format, progress_stat=True):
		""" Function that parses multiple alignment/loci files and returns a dictionary with the taxa as keys and sequences as values as well as two integers corresponding to the number of taxa and sequence length """
		loci_lengths = [] # Saves the sequence lengths of the 
		loci_range = [] # Saves the loci names as keys and their range as values
		self.main_taxa_order = []
		for infile in input_list:
			if progress_stat == True: # When set to True, this statement produces a progress status on the terminal
				print ("\rProcessing file %s out of %s" % (input_list.index(infile)+1,len(input_list)),end="")
			current_alignment, taxa_order, current_sequence_len, loci_range_temp = self.read_alignment(infile,alignment_format) # Parse the current alignment
			current_sequence_len = "".join("%s" % x for x in current_sequence_len) # Due to technical reasons that increase code homogeneity, the "current_sequence_len" is a list by default. Therefore, for this particular function, I need to convert it to a string
			if loci_lengths == []:
				main_alignment = current_alignment # Create the main alignment dictionary from the first current alignment and never visit this statement again
				self.main_taxa_order = taxa_order
				loci_lengths.append(int(current_sequence_len))
				loci_range.append((infile.split(".")[0],"1-%s" % (int(current_sequence_len)-1))) # Saving the range for the first loci
			else:
				for taxa, sequence in current_alignment.items(): 
					if taxa in main_alignment: 
						main_alignment[taxa] += sequence # Append the sequence from the current alignment to the respective taxa in the main alignment
					elif taxa not in main_alignment:
						main_alignment[taxa] = self.missing*sum(loci_lengths)+sequence # If the taxa does not yet exist in the main alignment, create the new entry with a sequence of 'n' characters of the same size as the length of the missed loci and the sequence from the current alignment
						self.main_taxa_order.append(taxa)
				loci_range.append((infile.split(".")[0],"%s-%s" % (loci_lengths[-1], int(current_sequence_len)-1))) # Saving the range for the subsequent loci
				loci_lengths.append(int(current_sequence_len))
				for taxa in main_alignment.keys():
					if taxa not in current_alignment: # Check if any taxa from the main alignment are missing from the current alignment. If yes, fill them with 'n'
						main_alignment[taxa] += self.missing*int(current_sequence_len)
		self.main_taxa_order = self.taxa_order
		return (main_alignment, self.taxa_order, loci_lengths, loci_range)
		
	def check_sizes (self, alignment_dic, current_file):
		""" Function to test whether all sequences are of the same size and, if not, which are different """
		# Determine the most common length 
		commonSeq = max(set([v for v in alignment_dic.values()]),key=[v for v in alignment_dic.values()].count)
		# Creates a dictionary with the sequences, and respective length, of different length
		difLength = dict((key,value) for key, value in alignment_dic.items() if len(commonSeq) != len(value))
		if difLength != {}:
			print ("\nWARNING: Unequal sequence lenght detected in %s for the following taxa" % i_file)
			
	def zorro2rax (self, alignment_file_list, zorro_sufix="_zorro.out"):
		""" Function that converts the floating point numbers contained in the original zorro output files into intergers that can be interpreted by RAxML. If multiple alignment files are provided, it also concatenates them in the same order """
		weigths_storage = []
		for alignment_file in alignment_file_list:
			zorro_file = alignment_file.split(".")[0]+zorro_sufix # This assumes that the prefix of the alignment file is shared with the corresponding zorro file
			zorro_handle = open(zorro_file)
			weigths_storage += [round(float(weigth.strip())) for weigth in zorro_hanlde]
		return zorro_storage
	
	class writer ():
			
		def __init__ (self, output_file, taxa_order, coding, loci_lengths, loci_range, gap = "-", missing = "n", conversion = None):
			self.output_file = output_file
			self.taxa_order = taxa_order
			self.coding = coding
			self.loci_lengths = loci_lengths
			self.loci_range = loci_range
			self.gap = gap
			self.missing = missing
			# The space (in characters) available for the taxon name before the sequence begins
			self.seq_space_nex = 40
			self.seq_space_phy = 30
			self.seq_space_ima2 = 10
			# Cut the taxa names by the following character:
			self.cut_space_nex = 50
			self.cut_space_phy = 50
			self.cut_space_ima2 = 8
				
		def phylip (self, alignment_dic, conversion=None):
			""" Writes a pre-parsed alignment dictionary into a new phylip file """
			out_file = open(self.output_file+".phy","w")
			out_file.write("%s %s\n" % (len(alignment_dic), sum(self.loci_lengths)))
			for key in self.taxa_order:
					out_file.write("%s %s\n" % (key[:self.cut_space_phy].ljust(self.seq_space_phy),alignment_dic[key]))
			if conversion == None:
				partition_file = open(self.output_file+"_part.File","a")
				for partition,lrange in self.loci_range:
					partition_file.write("%s, %s = %s\n" % (self.coding,partition,lrange))
			out_file.close()
					
		def fasta (self, alignment_dic, conversion=None):
			""" Writes a pre-parsed alignment dictionary into a new fasta file """
			out_file = open(self.output_file+".fas","w")
			for key in self.taxa_order:
				out_file.write(">%s\n%s\n" % (key,alignment_dic[key]))
			out_file.close()
				
		def nexus (self, alignment_dic, conversion=None):
			""" Writes a pre-parsed alignment dictionary into a new nexus file """
			out_file = open(self.output_file+".nex","w")
			out_file.write("#NEXUS\n\nBegin data;\n\tdimensions ntax=%s nchar=%s ;\n\tformat datatype=%s interleave=no gap=%s missing=%s ;\n\tmatrix\n" % (len(alignment_dic), sum(self.loci_lengths), self.coding, self.gap, self.missing))
			for key in self.taxa_order:
				out_file.write("%s %s\n" % (key[:self.cut_space_nex].ljust(self.seq_space_nex),alignment_dic[key]))
			out_file.write(";\n\tend;")
			if conversion == None:
				out_file.write("\nbegin mrbayes;\n")
				for partition,lrange in self.loci_range:
					out_file.write("\tcharset %s = %s;\n" % (partition,lrange))
				out_file.write("\tpartition part = %s: %s;\nend;" % (len(self.loci_range),"".join([part[0] for part in self.loci_range])))
			out_file.close()
				
		def zorro (self, zorro_weigths):
			""" Creates a concatenated file with the zorro weigths for the corresponding alignment files """
			outfile = self.output_file+"_zorro.out"
			for weigth in zorro_weigths:
				outfile.write("%s\n" % weigth)
			outfile.close()
